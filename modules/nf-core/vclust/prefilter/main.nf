process VCLUST_PREFILTER {
    tag "${meta.id}"
    label 'process_high'

    input:
    tuple val(meta), path(fasta)
    path vclust

    output:
    tuple val(meta), path("${prefix}.vclust_prefilter.txt") , emit: prefilter
    path "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    ${vclust}/vclust.py \\
        prefilter \\
        --in ${fasta} \\
        --out ${prefix}.vclust_prefilter.txt \\
        --threads ${task.cpus} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vclust: \$( ${vclust}/vclust.py --version | sed 's/v//' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.vclust_prefilter.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vclust: \$( ${vclust}/vclust.py --version | sed 's/v//' )
    END_VERSIONS
    """
}
