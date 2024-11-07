process MGEFINDER_PAIR {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://carsonjm/mgefinder:0.0.1' :
        'docker://carsonjm/mgefinder:0.0.1' }"

    input:
    tuple val(meta) , path(fasta)
    tuple val(meta2), path(bam)
    tuple val(meta3), path(find)

    output:
    tuple val(meta), path("${prefix}.mgefinder.pair")   , emit: pair
    path "versions.yml"                                 , emit: versions

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    mgefinder \\
        pair \\
        ${find} \\
        ${bam[0]} \\
        ${find} \\
        -o ${prefix}.mgefinder.pair \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mgefinder: \$(echo \$( mgefinder --version | tail -n 1 | sed "s/MGEfinder version://" ))
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.mgefinder.pair

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mgefinder: \$(echo \$( mgefinder --version | tail -n 1 | sed "s/MGEfinder version://" ))
    END_VERSIONS
    """
}
