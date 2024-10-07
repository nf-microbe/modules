process SEQHASHER_SEQHASHER {
    tag "${meta.id}"
    label 'process_high'

    input:
    tuple val(meta), path(fasta)
    path seq_hasher

    output:
    tuple val(meta), path("${prefix}.seq-hasher.tsv")   , emit: tsv
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    ${seq_hasher} \\
        ${fasta} \\
        ${args} \\
        > ${prefix}.seq-hasher.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqhasher: \$(${seq_hasher} --version | sed 's/seq-hasher //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.seq-hasher.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqhasher: \$(${seq_hasher} --version | sed 's/seq-hasher //')
    END_VERSIONS
    """
}
