process MGEFINDER_FORMATBAM {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://carsonjm/mgefinder:0.0.1' :
        'docker://carsonjm/mgefinder:0.0.1' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${prefix}.mgefinder.bam*")   , emit: bam
    path "versions.yml"                                 , emit: versions

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    mgefinder \\
        formatbam \\
        ${bam} \\
        ${prefix}.mgefinder.bam \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mgefinder: \$(echo \$( mgefinder --version | tail -n 1 | sed "s/MGEfinder version://" ))
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.mgefinder.bam
    touch ${prefix}.mgefinder.bam.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mgefinder: \$(echo \$( mgefinder --version | tail -n 1 | sed "s/MGEfinder version://" ))
    END_VERSIONS
    """
}
