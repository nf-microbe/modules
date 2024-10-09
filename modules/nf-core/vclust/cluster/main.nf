process VCLUST_CLUSTER {
    tag "${meta.id}"
    label 'process_high'

    input:
    tuple val(meta) , path(ani)
    tuple val(meta2), path(ids)
    path vclust

    output:
    tuple val(meta), path("${prefix}.vclust_clusters.tsv")  , emit: clusters
    path "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    ${vclust}/vclust.py \\
        cluster \\
        --in ${ani} \\
        --ids ${ids} \\
        --out ${prefix}.vclust_clusters.tsv \\
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
    touch ${prefix}.vclust_clusters.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vclust: \$( ${vclust}/vclust.py --version | sed 's/v//' )
    END_VERSIONS
    """
}
