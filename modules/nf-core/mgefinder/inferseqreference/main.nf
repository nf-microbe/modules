process MGEFINDER_INFERSEQREFERENCE {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'shub://bhattlab/MGEfinder-singularity:latest' :
        'shub://bhattlab/MGEfinder-singularity:latest' }"

    input:
    tuple val(meta) , path(pair)
    tuple val(meta2), path(fasta)

    output:
    tuple val(meta), path("${prefix}.mgefinder.inferseq_reference.tsv") , emit: tsv
    path "versions.yml"                                                 , emit: versions

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    mgefinder \\
        inferseq-reference \\
        ${pair} \\
        ${fasta} \\
        -o ${prefix}.mgefinder.inferseq_reference.tsv \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mgefinder: \$(echo \$( prophage_tracer_WGS.sh -h | sed -n "s/Prophage Tracer V//; s/\s.*//;  2p" ))
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.mgefinder.inferseq_reference.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mgefinder: \$(echo \$( prophage_tracer_WGS.sh -h | sed -n "s/Prophage Tracer V//; s/\s.*//;  2p" ))
    END_VERSIONS
    """
}
