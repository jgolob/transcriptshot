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

// containers in use
container__prodigal = 'quay.io/biocontainers/prodigal:2.6.3--h516909a_2'
container__diamond = 'quay.io/fhcrc-microbiome/docker-diamond:v2.0.6-biopython'
container__pandas = "quay.io/fhcrc-microbiome/python-pandas:v1.0.3"
container__biopython = 'quay.io/biocontainers/biopython:1.78'

// Using DSL-2
nextflow.enable.dsl=2

// Default parameters
params.cluster_min_identity = 100
params.output = './output'


workflow GenomeToAlleles {
    take: genome_manifest_ch

    main:
        // For each genome in the manifest extract protein coding genes with prodigal
        Prodigal(
            genome_manifest_ch
        )
        Extract_GB_CDS(
            genome_manifest_ch
        )
        Extract_GB_CDS.out.view()
        // Concatenate the per-genome catalogs into one catalog with each gene assigned a unique name
        CombineAndRenameAlleles(
            Prodigal.out[0].collect{ it[0] },
            Prodigal.out[0].collect{ it[1] },
            Prodigal.out[0].collect{ it[2] },
        )
        // Deduplicate the resulting catalog using diamond
        DiamondDedup(
            CombineAndRenameAlleles.out[0]
        )
    
        JoinMap(
            DiamondDedup.out[1],
            CombineAndRenameAlleles.out[1]
        )


// */
}


//
// Reference genomes
//
// Parse gb of genomes
// Prodigal of raw sequence for protein coding sequences
// Convert to diamond DB

def ReadGenomeManifest(genome_manifest_file){
    genome_manifest_file.splitCsv(
        header: true, 
        sep: ","
    ).branch{
        valid:  (it.organism != null) && (it.organism_shortname != null) && (it.path != null ) && (it.path != "" ) && (!file(it.path).isEmpty())
        other: true
    }
}

process Extract_GB_CDS {
    tag "Extract annotated CDS from genbank genomes"
    container "${container__biopython}"
    label 'io_limited'
    publishDir "${params.output}/gene_catalog/${shortname}/", mode: "copy"

    input:
        tuple val(organism), val(shortname), file(genbank_f)

    output: 
        tuple val(organism), val(shortname), path("${shortname}.fasta"), path("${shortname}.cds.csv")
"""
#!/usr/bin/env python

from Bio import SeqIO
import gzip
import csv

srs = SeqIO.parse(
    gzip.open(
        '${genbank_f}'
        , 'rt')
, 'gb')

locus_details = []
with open('${shortname}.fasta', 'wt') as out_h:
    for sr in srs:
        for feat in sr.features:
            if feat.type != 'CDS':
                continue
            # Implicit else
            gene = feat.qualifiers.get('gene',[''])[0]
            locus_tag = feat.qualifiers.get('locus_tag',[''])[0]
            if locus_tag == "":
                print("missing tag")
            out_h.write(">{} {}\\n{}\\n".format(
                locus_tag,
                "; ".join(
                    feat.qualifiers.get('gene', [])+
                    feat.qualifiers.get('product', [])
                ),
                feat.extract(sr.seq)
            ))
            locus_details.append({
                'locus_tag': locus_tag,
                'gene': '; '.join(feat.qualifiers.get('gene', [])),
                'product': '; '.join(feat.qualifiers.get('product', []))
            })

with open('${shortname}.cds.csv', 'wt') as out_h:
    out_w = csv.DictWriter(out_h, fieldnames=['locus_tag', 'gene', 'product'])
    out_w.writerows(locus_details)

"""
}

process Prodigal {
    tag "Identify protein-coding genes"
    container "${container__prodigal}"
    label 'io_limited'
    //publishDir "${params.output_folder}/assembly/${specimen}/", mode: "copy"

    input:
        tuple val(organism), val(shortname), file(contigs)
    
    output:
        tuple val(organism), val(shortname), file("${shortname}.faa.gz")
        file("${shortname}.gff.gz")
    
"""
set -Eeuo pipefail 

gunzip -c ${contigs} > ${shortname}.contigs.fna

prodigal \
    -a ${shortname}.faa \
    -i  ${shortname}.contigs.fna \
    -f gff \
    -o ${shortname}.gff \
    -p single

gzip ${shortname}.gff
gzip ${shortname}.faa

"""
}

process DiamondDedup {
    tag "Deduplicate sequences by alignment with DIAMOND"
    container "${container__diamond}"
    label 'multithread'
    errorStrategy 'finish'
    publishDir "${params.output}/", mode: "copy"
    
    input:
        file(input_genes)
    
    output:
        file("output.genes.fasta.gz")
        file("output_map.csv.gz")
    
"""
#!/bin/bash

set -e

# Make a DIAMOND database
diamond \
    makedb \
    --in ${input_genes} \
    --db input.genes.dmnd \
    --threads ${task.cpus}

# Align the genes against themselves and filter any which align
diamond \
    blastp \
    --query ${input_genes} \
    --threads ${task.cpus} \
    --db input.genes.dmnd \
    --out input.genes.aln \
    --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen \
    --query-cover 90 \
    --id ${params.cluster_min_identity} \
    --top 0 \
    --block-size ${task.memory.toMega() / (1024 * 6)} \
    --unal 0

python3 << END

from Bio.SeqIO.FastaIO import SimpleFastaParser
import gzip
import sys
import csv

# Keep track of the genes which have been filtered out
# Key is the duplicated gene id. Value is the gene into which it was combined
duplicate_genes = {}

# Iterate over the alignment
print("Reading alignments")
ix = 0
for line in open("input.genes.aln", "r"):

    ix += 1

    if ix % 100000 == 0:
        print("Read %d lines of alignment - found %d duplicated genes" % (ix, len(duplicate_genes)))

    qname, sname, _, _, _, _, _, _, _, _, _, _, qlen, slen = line.rstrip("\\n").split("\\t")

    # Skip self alignments
    if qname == sname:
        continue
    # If we have excluded either of the genes before, skip the line
    if qname in duplicate_genes or sname in duplicate_genes:
        continue
    # For non-self alignments, remove the smaller of the two
    if int(slen) < int(qlen):
        duplicate_genes[sname] = qname
    else:
        duplicate_genes[qname] = sname

assert ix > 0, "Didn't read any alignments"
print("Read %d lines of alignment - found %d duplicated genes" % (ix, len(duplicate_genes)))
print("Done reading alignments")

# Now let's make the filtered FASTA with all duplicated genes removed
n_found = 0
n = 0
with gzip.open('${input_genes}', "rt") as handle_in:
    with gzip.open("output.genes.fasta.gz", "wt") as handle_out, gzip.open("output_map.csv.gz", 'wt') as map_h:
        map_w = csv.writer(map_h)
        map_w.writerow(['seq_id', 'centroid_id'])
        for header, seq in SimpleFastaParser(handle_in):

            header = header.split(" ", 1)[0]

            n += 1
            if header in duplicate_genes:
                map_w.writerow([header, duplicate_genes[header]])
                n_found += 1
            else:
                handle_out.write(">%s\\n%s\\n" % (header, seq))
                map_w.writerow([header, header])

# Make sure that we encountered all of the duplicated genes
print("Read in %d sequences, filtered out %d, wrote out the rest" % (n, n_found))
assert n_found == len(duplicate_genes), "%d != %d" % (n_found, len(duplicate_genes))

END

echo "Done"
"""
}

// Combine and Assign a new, shorter name to a set of alleles
process CombineAndRenameAlleles {
    tag "Make concise unique gene names"
    container "${container__diamond}"
    label 'io_limited'

    input:
    val(organisms)
    val(shortnames)
    file(genomes_faa_gz)

    output:
    file "alleles.faa.gz"
    file "allele_id_map.csv.gz"

"""
#!/usr/bin/env python3

from Bio.SeqIO.FastaIO import SimpleFastaParser
import gzip
import uuid
import csv

def random_string(n=8):
    return str(uuid.uuid4())[:n]

used_strings = set([])

organisms = "${organisms.join(';;')}".split(';;')
shortnames = "${shortnames.join(';;')}".split(';;')
org_faa_fn =  "${genomes_faa_gz.join(';;')}".split(';;')

old_new_id_long = []
with gzip.open("alleles.faa.gz", "wt") as fo:
    for org, sn, in_fn in zip(organisms, shortnames, org_faa_fn):
        with gzip.open(in_fn, "rt") as fi:
            org_names = []
            for header, seq in SimpleFastaParser(fi):
                old_id = header.split(" ")[0]
                new_string = "{}__{}".format(sn, random_string())
                while new_string in used_strings:
                    new_string = "{}__{}".format(sn, random_string())
                used_strings.add(new_string)
                org_names.append({
                    'organism': org,
                    'organism_shortname': sn,
                    'old_id': old_id,
                    'seq_id': new_string
                })
                fo.write(">{} {}\\n{}\\n".format(
                    new_string,
                    org,
                    seq
                ))
            old_new_id_long += org_names
# Output mapping of organism, sn, old_id, new_id
with gzip.open("allele_id_map.csv.gz", 'wt') as map_h:
    map_w = csv.DictWriter(
        map_h,
        fieldnames=[
            'organism',
            'organism_shortname',
            'old_id',
            'seq_id'
    ])
    map_w.writeheader()
    for row in old_new_id_long:
        map_w.writerow(row)

"""
}

// Combine and Assign a new, shorter name to a set of alleles
process JoinMap {
    tag "Join source gene IDs to centroids"
    container "${container__pandas}"
    label 'io_limited'
    publishDir "${params.output}/", mode: "copy"

    input:
        file(centroid_map)
        file(org_map)

    output:
    file "genome_allele_gene_map.csv.gz"

"""
#!/usr/bin/env python3
import pandas as pd

centroid_map = pd.read_csv("${centroid_map}")
org_map = pd.read_csv("${org_map}")

full_map = pd.merge(
    org_map,
    centroid_map,
    on='seq_id',
    how='inner'
)

full_map.to_csv("genome_allele_gene_map.csv.gz", index=None)

"""
}



// Function which prints help message text
def helpMessage() {
    log.info"""
    Workflow to convert reference genomes to a gene catalog

    Usage:
    --manifest      Path to a CSV with the minimum of these columns, one row per genome to be analyzed (REQUIRED)
                    'organism'                  A person-friendly name for this organism 
                                                (e.g. Ruminococcus bromii ATCC 27255)
                    'organism_shortname'        A computer-friendly (i.e. no spaces or symbols) short name. 
                                                MUST be unique for this analysis (e.g. RB_ATCC27255)
                    'path'                      Path to a gzipped fasta (FNA) or genbank file of the genome and/or contigs.
    
    --output        Path where to place the output files:
                        output.genes.fasta.gz                       A catalog of genes (centroids)
                        genome_allele_gene_map.csv.gz               A mapping of organism - shortname - gene_id - centroid_id
    
    Parameters:
    --cluster_min_identity                      Min PERCENT identity for clustering stage. (Default 100)

    """.stripIndent()
}

workflow {
    if (params.manifest == null) {
        helpMessage()
        exit 0
    }

    // Load manifest!
    manifest = ReadGenomeManifest(
        Channel.from(
            file(params.manifest)
        )
    )
    GenomeToAlleles(
        manifest.valid.map{
            [
                it.organism,
                it.organism_shortname,
                file(it.path)
            ]        
        }
    )
    manifest.other.view()

}