process MGEFINDER_INFERSEQREFERENCE {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://carsonjm/mgefinder:0.0.1' :
        'docker://carsonjm/mgefinder:0.0.1' }"

    input:
    tuple val(meta) , path(fasta)
    tuple val(meta2), path(pair)

    output:
    tuple val(meta), path("${prefix}.mgefinder.inferseq_reference.tsv") , emit: tsv
    path "versions.yml"                                                 , emit: versions

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    def is_compressed = fasta.getExtension() == "gz" ? true : false
    def fasta_name = is_compressed ? fasta.getBaseName() : fasta
    """
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${fasta} > ${fasta_name}
    fi

    mgefinder \\
        inferseq-reference \\
        ${pair} \\
        ${fasta_name} \\
        -o ${prefix}.mgefinder.inferseq_reference.tsv \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mgefinder: \$(echo \$( mgefinder --version | tail -n 1 | sed "s/MGEfinder version://" ))
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.mgefinder.inferseq_reference.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mgefinder: \$(echo \$( mgefinder --version | tail -n 1 | sed "s/MGEfinder version://" ))
    END_VERSIONS
    """
}
