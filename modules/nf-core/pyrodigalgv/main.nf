process PYRODIGALGV {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pyrodigal-gv:0.3.2--pyh7e72e81_0':
        'biocontainers/pyrodigal-gv:0.3.2--pyh7e72e81_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${prefix}.pyrodigalgv.gff.gz")   , emit: gff
    tuple val(meta), path("${prefix}.pyrodigalgv.faa.gz")   , emit: faa
    tuple val(meta), path("${prefix}.pyrodigalgv.fna.gz")   , emit: fna
    path "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    pyrodigal-gv \\
        -i ${fasta} \\
        -a ${prefix}.pyrodigalgv.faa \\
        -d ${prefix}.pyrodigalgv.fna \\
        -f gff \\
        -o ${prefix}.pyrodigalgv.gff \\
        --jobs ${task.cpus} \\
        ${args}

    gzip ${prefix}.pyrodigalgv.faa ${prefix}.pyrodigalgv.fna ${prefix}.pyrodigalgv.gff

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pyrodigalgv: \$(pyrodigal-gv --version | sed 's/pyrodigal-gv v//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.pyrodigalgv.faa.gz
    echo "" | gzip > ${prefix}.pyrodigalgv.fna.gz
    echo "" | gzip > ${prefix}.pyrodigalgv.gff.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pyrodigalgv: \$(pyrodigal-gv --version | sed 's/pyrodigal-gv v//')
    END_VERSIONS
    """
}
