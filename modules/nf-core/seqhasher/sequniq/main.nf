process SEQHASHER_SEQUNIQ {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:5.3.0':
        'biocontainers/gawk:5.3.0' }"

    input:
    tuple val(meta), path(seq_hasher_tsv)

    output:
    tuple val(meta), path("${prefix}.dups.tsv")     , optional: true    , emit: dups
    tuple val(meta), path("${prefix}.seq-uniq.tsv")                     , emit: unique
    path "versions.yml"                                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    awk 'seen[\$2]++ {
        print \$1 > "${prefix}.dups.tsv"; next
    }
    {
        print \$1 > "${prefix}.seq-uniq.tsv"
    }' ${seq_hasher_tsv}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$( awk --version| head -n 1 | sed 's/GNU Awk //; s/, .*//' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.dups.tsv
    touch ${prefix}.seq-uniq.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$( awk --version| head -n 1 | sed 's/GNU Awk //; s/, .*//' )
    END_VERSIONS
    """
}
