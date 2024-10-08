process VIROMEQC_VIROMEQC {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-b28a1a551d380ce8d57f9d83894ccb9559b44404:08a4cf815fcef0080ede5b3633202cbda8edf59b-0':
        'biocontainers/mulled-v2-b28a1a551d380ce8d57f9d83894ccb9559b44404:08a4cf815fcef0080ede5b3633202cbda8edf59b-0' }"

    input:
    tuple val(meta), path(fastq)
    path(viromeqc_db)

    output:
    tuple val(meta), path("${prefix}.viromeqc.tsv") , emit: tsv
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    viromeQC.py \\
        --input ${fastq} \\
        --index_dir ${viromeqc_db} \\
        --output ${prefix}.viromeqc.tsv \\
        --bowtie2_threads ${task.cpus} \\
        --diamond_threads ${task.cpus} \\
        --sample_name ${prefix} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        viromeqc: \$(echo \$(viromeQC.py --version 2>&1 | sed -n 's/^.*Version:\\s//; 1p' ))
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.viromeqc.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        viromeqc: \$(echo \$(viromeQC.py --version 2>&1 | sed -n 's/^.*Version:\\s//; 1p' ))
    END_VERSIONS
    """
}
