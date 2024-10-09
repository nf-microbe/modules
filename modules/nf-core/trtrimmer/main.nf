process TRTRIMMER {
    tag "${meta.id}"
    label 'process_high'

    input:
    tuple val(meta), path(fasta)
    path tr_trimmer

    output:
    tuple val(meta), path("${prefix}.tr-trimmer.fasta.gz")  , emit: fasta
    path "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    ${tr_trimmer} \\
        ${fasta} \\
        ${args} \\
        > ${prefix}.tr-trimmer.fasta

    gzip ${prefix}.tr-trimmer.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        trtrimmer: \$(tr-trimmer --version | sed 's/tr-trimmer //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.tr-trimmer.fasta.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        trtrimmer: \$(tr-trimmer --version | sed 's/tr-trimmer //')
    END_VERSIONS
    """
}
