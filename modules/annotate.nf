
container__eggnogmapper = 'quay.io/biocontainers/eggnog-mapper:2.0.1--py_1'

// Using DSL-2
nextflow.enable.dsl=2


workflow AnnotateAlleles {
    take: 
        alleles_fastp_gz
        eggnog_db
        eggnog_proteins_dmnd

    main:
        Eggnog(
            alleles_fastp_gz,
            eggnog_db,
            eggnog_proteins_dmnd
        )
        Eggnog.out.view()
// */
}

process Eggnog {
    tag "Annotate genes by predicted function"
    container "${container__eggnogmapper}"
    label 'multithread'
    publishDir "${params.output}/", mode: "copy"
    
    input:
    path query
    path eggnog_db
    path eggnog_dmnd

    output:
    path "genes.emapper.annotations.gz"

    
    """
set -Eeuo pipefail

mkdir data
mkdir TEMP
mkdir SCRATCH

cp ${eggnog_db} data/eggnog.db
cp ${eggnog_dmnd} data/eggnog_proteins.dmnd

emapper.py \
    -i ${query} \
    --output genes \
    -m "diamond" \
    --cpu ${task.cpus} \
    --data_dir data/ \
    --scratch_dir SCRATCH/ \
    --temp_dir TEMP/

gzip genes.emapper.annotations
    
    """
}

// Function which prints help message text
def helpMessage() {
    log.info"""
    Annotate an allele catalog

    Usage:
    --alleles       Path to alleles in fastp.gz format
    --eggnog_db     Eggnog DB
    -eggnog_dmnd    Eggnog Proteins in diamond      
    
    --output        Path where to place the output files:
                    
    

    """.stripIndent()
}

workflow {
    if (
        (params.alleles == null) |
        (params.eggnog_db == null) |
        (params.eggnog_dmnd == null) |
        (params.help != null)
    ) {
        helpMessage()
        exit 0
    }

    AnnotateAlleles(
        file(params.alleles),
        file(params.eggnog_db),
        file(params.eggnog_dmnd)
    )

    

}