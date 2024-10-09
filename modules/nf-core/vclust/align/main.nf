process VCLUST_ALIGN {
    tag "${meta.id}"
    label 'process_high'

    input:
    tuple val(meta) , path(fasta)
    tuple val(meta2), path(prefilter)
    path vclust

    output:
    tuple val(meta), path("${prefix}.vclust_ani.tsv")       , emit: ani
    tuple val(meta), path("${prefix}.vclust_ani.ids.tsv")   , emit: ids
    tuple val(meta), path("${prefix}.vclust_align.tsv")     , emit: align
    path "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    ${vclust}/vclust.py \\
        align \\
        --in ${fasta} \\
        --out ${prefix}.vclust_ani.tsv \\
        --out-aln ${prefix}.vclust_align.tsv \\
        --filter ${prefilter} \\
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
    touch ${prefix}.vclust_ani.tsv
    touch ${prefix}.vclust_align.tsv
    touch ${prefix}.vclust_ani.ids.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vclust: \$( ${vclust}/vclust.py --version | sed 's/v//' )
    END_VERSIONS
    """
}
