process SEQKIT_CONCAT {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqkit:2.8.1--h9ee0642_0':
        'biocontainers/seqkit:2.8.1--h9ee0642_0' }"

    input:
    tuple val(meta), path(input, stageAs: 'in/*')

    output:
    tuple val(meta), path("${prefix}.*")    , emit: fastx
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args         ?: ""
    prefix      = task.ext.prefix       ?: "${meta.id}"
    def extension   = "fastq"
    if (input[0] ==~ /.+\.fasta|.+\.fasta.gz|.+\.fa|.+\.fa.gz|.+\.fas|.+\.fas.gz|.+\.fna|.+\.fna.gz|.+\.fsa|.+\.fsa.gz/ ) {
        extension   = "fasta"
    }
    extension       = input[0].toString().endsWith('.gz') ? "${extension}.gz" : extension
    """
    seqkit \\
        concat \\
        $args \\
        in/* > ${prefix}.${extension}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$(seqkit version | cut -d' ' -f2)
    END_VERSIONS
    """

    stub:
    prefix  = task.ext.prefix   ?: "${meta.id}"
    def extension   = "fastq"
    if (input[0] ==~ /.+\.fasta|.+\.fasta.gz|.+\.fa|.+\.fa.gz|.+\.fas|.+\.fas.gz|.+\.fna|.+\.fna.gz|.+\.fsa|.+\.fsa.gz/ ) {
        extension   = "fasta"
    }
    extension       = input[0].toString().endsWith('.gz') ? "${extension}.gz" : extension
    """
    touch ${prefix}.${extension}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$(seqkit version | cut -d' ' -f2)
    END_VERSIONS
    """
}
