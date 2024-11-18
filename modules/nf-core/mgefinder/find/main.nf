process MGEFINDER_FIND {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://carsonjm/mgefinder:0.0.1' :
        'docker://carsonjm/mgefinder:0.0.1' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${prefix}.mgefinder.find")   , emit: find
    path "versions.yml"                                 , emit: versions

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    mgefinder \\
        find \\
        ${bam[0]} \\
        -o ${prefix}.mgefinder.find \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mgefinder: \$(echo \$( mgefinder --version | tail -n 1 | sed "s/MGEfinder version://" ))
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.mgefinder.find

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mgefinder: \$(echo \$( mgefinder --version | tail -n 1 | sed "s/MGEfinder version://" ))
    END_VERSIONS
    """
}
