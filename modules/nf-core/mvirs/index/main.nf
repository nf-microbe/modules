process MVIRS_INDEX {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mvirs:1.1.1--pyhdfd78af_0' :
        'biocontainers/mvirs:1.1.1--pyhdfd78af_0' }"

    input:
    tuple val(meta) , path(fasta)

    output:
    tuple val(meta), path("${fasta}.*") , emit: index
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    mvirs \\
        index \\
        -f ${fasta}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mvirs: \$(echo \$(mvirs 2>&1) | sed 's/.*Version: //; s/ .*//')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${fasta}.amb
    touch ${fasta}.ann
    touch ${fasta}.bwt
    touch ${fasta}.pac
    touch ${fasta}.sa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mvirs: \$(echo \$(mvirs 2>&1) | sed 's/.*Version: //; s/ .*//')
    END_VERSIONS
    """
}
