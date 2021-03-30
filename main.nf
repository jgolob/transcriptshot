#!/usr/bin/env nextflow

/*
  transcriptshot: A pipeline to robustly identify which alleles (n.e.e peptide coding sequences)
  are transcribed in a microbial community.

  This is intended for microbial transcriptional analysis, including metatranscriptomes.
  Initial version is intended for defined communities, for which every microbe present has an available genome 
  (does not need to be circularized).

  Attempts to find co-transcribed genes as well.
  Related to geneshot.
*/

container__diamond = 'quay.io/fhcrc-microbiome/docker-diamond:v2.0.6-biopython'
container__pandas = "quay.io/fhcrc-microbiome/python-pandas:v1.0.3"
container__barcodecop = "golob/barcodecop:0.5__bc_1"
container__trimgalore = 'quay.io/biocontainers/trim-galore:0.6.6--0'


// Default parameters
params.output = './output'

params.dmnd_min_identity = 80 // DIAMOND
params.dmnd_min_coverage = 50 // DIAMOND
params.dmnd_top_pct = 1 // DIAMOND
params.dmnd_min_score = 20 // DIAMOND
params.gencode = 11 //DIAMOND
params.sd_mean_cutoff = 3.0 // FAMLI
params.famli_batchsize = 10000000 // FAMLI
params.famli_folder = false // Import FAMLI-filtered alignments



// Using DSL-2
nextflow.enable.dsl=2



workflow TranscriptToGenes {
    take: 
        paired_read_manifest_ch
        failed_read_ch

    main:

    // trimgalore to clean up / trim
    TrimGalore(
        paired_read_manifest_ch.map{
            [it.specimen, file(it.R1), file(it.R2)]
        }
    )
    // Make DB of allele catalog
    DiamondDB(
        file(params.genecat)
    )
    // Align with diamond
    DiamondALN(
        TrimGalore.out,
        DiamondDB.out,
    )

    // FAMLI to cleanup alignments
    Famli(
        DiamondALN.out
    )
    
    JoinFamli(
        Famli.out.collect{ it[0] },
        Famli.out.collect{ it[1] },
    )

    JoinFamli.out
    
    
    // Output: Transcript counts long
    // CTG with ANN / cosine distance
    // Output: gene <-> CTG link
    // Output: GTC <-> specimen <-> RC long
    // Function which prints help message text


// */
}

// Use trim_galore to handle adapters / etc
process TrimGalore {
    container "${container__trimgalore}"
    label 'io_limited'
    errorStrategy 'finish'

    input:
    tuple val(specimen), file(R1), file(R2)

    output:
    tuple val(specimen), file("${specimen}.R1.tg.fastq.gz"), file("${specimen}.R2.tg.fastq.gz")

    """
    set -e

    cp ${R1} R1.fastq.gz
    cp ${R2} R2.fastq.gz

    trim_galore \
    --gzip \
    --cores ${task.cpus} \
    --paired \
    --fastqc \
    R1.fastq.gz R2.fastq.gz

    rm R1.fastq.gz
    rm R2.fastq.gz
    mv R1_val_1.fq.gz "${specimen}.R1.tg.fastq.gz"
    mv R2_val_2.fq.gz "${specimen}.R2.tg.fastq.gz"
    """
}
process DiamondDB {
    tag "Make a DIAMOND database"
    container "quay.io/fhcrc-microbiome/docker-diamond:v2.0.6-biopython"
    label 'mem_veryhigh'
    publishDir "${params.output}/ref/", mode: "copy"
    
    input:
    file fasta

    output:
    file "genes.dmnd"

    """
    set -e
    diamond \
      makedb \
      --in ${fasta} \
      --db genes.dmnd \
      --threads ${task.cpus}
    """

}

// Align each sample against the reference database of genes using DIAMOND
process DiamondALN {
    tag "Align to the gene catalog"
    container "${container__diamond}"
    label 'mem_veryhigh'
    
    input:
    tuple val(sample_name), file(R1), file(R2)
    file refdb
    
    output:
    tuple val(sample_name), file("${sample_name}.aln.gz")

    """
    set -e

    for fp in ${R1} ${R2}; do
        cat \$fp | \
        gunzip -c | \
        awk "{if(NR % 4 == 1){print \\"@\$fp\\" NR }else{if(NR % 4 == 3){print \\"+\\"}else{print}}}" | \
        gzip -c \
        > query.fastq.gz

        diamond \
        blastx \
        --query query.fastq.gz \
        --out \$fp.aln.gz \
        --threads ${task.cpus} \
        --db ${refdb} \
        --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen \
        --min-score ${params.dmnd_min_score} \
        --query-cover ${params.dmnd_min_coverage} \
        --id ${params.dmnd_min_identity} \
        --top ${params.dmnd_top_pct} \
        --block-size ${task.memory.toMega() / (1024 * 6 * task.attempt)} \
        --query-gencode ${params.gencode} \
        --compress 1 \
        --unal 0
    done

    cat *aln.gz > ${sample_name}.aln.gz
    """

}

// Filter the alignments with the FAMLI algorithm
process Famli {
    tag "Deduplicate multi-mapping reads"
    container "quay.io/fhcrc-microbiome/famli:v1.5"
    label 'mem_veryhigh'
    publishDir "${params.output}/abund/details/", mode: "copy"
    
    input:
    tuple val(sample_name), file(input_aln)
    
    output:
    tuple val(sample_name), file("${sample_name}.json.gz")

    """
    set -e
    famli \
      filter \
      --input ${input_aln} \
      --output ${sample_name}.json \
      --threads ${task.cpus} \
      --batchsize ${params.famli_batchsize} \
      --sd-mean-cutoff ${params.sd_mean_cutoff}
    gzip ${sample_name}.json
    """

}

// Combine FAMLI outputs to long-format of specimen <-> gene <-> count
process JoinFamli {
    tag "Combine FAMLI outputs to one file"
    container "${container__pandas}"
    label 'io_limited'
    publishDir "${params.output}/", mode: "copy"

    input:
        val(sample_ids)
        file(famli_jsons)

    output:
    file "specimen_gene_count_long.csv.gz"

"""
#!/usr/bin/env python3
import pandas as pd

sample_ids = "${sample_ids.join(";;")}".split(';;')
famli_json_files = "${famli_jsons.join(";;")}".split(';;')
sp_gene_long_list = []
for sp, js_fn in zip(sample_ids, famli_json_files):
    sp_js = pd.read_json(js_fn)[[
        'id',
        'nreads'
    ]].rename({'id': 'gene_id'}, axis=1)
    sp_js['specimen'] = sp
    sp_gene_long_list.append(sp_js)
pd.concat(sp_gene_long_list)[['specimen', 'gene_id', 'nreads']].to_csv(
    "specimen_gene_count_long.csv.gz",
    index=None
)
"""
}

def ReadManifest(manifest_file){
    manifest_file.splitCsv(
        header: true, 
        sep: ","
    ).branch{
        //valid_paired_indexed:  (it.specimen != null) && (it.R1 != null ) && (it.R1 != "" ) && (!file(it.R1).isEmpty()) && (it.R2 != null ) && (it.R2 != "") && (!file(it.R2).isEmpty()) && (it.I1 != null ) && (it.I1 != "" ) && (!file(it.I1).isEmpty()) && (it.I2 != null ) && (it.I2 != "") && (!file(it.I2).isEmpty())
        valid_paired:  (it.specimen != null) && (it.R1 != null ) && (it.R1 != "" ) && (!file(it.R1).isEmpty()) && (it.R2 != null ) && (it.R2 != "") && (!file(it.R2).isEmpty())
        //valid_unpaired:  (it.specimen != null) && (it.R1 != null ) && (it.R1 != "" ) && (!file(it.R1).isEmpty())
        other: true
    }
}

def helpMessage() {
    log.info"""
    Workflow to align RNA reads to a catalog of possible protein coding genes

    Usage:
    --manifest      Path to a read-pairs to be analyzed (REQUIRED)
                    'specimen'                  MUST be unique for this data set.
                    'R1'                        First read, in fasta format and gzipped
                    'R2'                        Second / reverse read, in fasta format and gzipped
    --genecat       gzipped FASTP file of possible alleles / protein coding genes (REQUIRED)
    
    --output        Path where to place the output files:
    
    Parameters:

    """.stripIndent()
}


workflow {
    if ((params.manifest == null) | (params.genecat == null)) {
        helpMessage()
        exit 0
    }

    // Load manifest!
    manifest = ReadManifest(
        Channel.from(
            file(params.manifest)
        )
    )

    TranscriptToGenes(
        manifest
    )

}