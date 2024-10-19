process HOSTILE_CLEAN {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hostile:1.1.0--pyhdfd78af_0':
        'biocontainers/hostile:1.1.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fastq)
    tuple val(meta2), path(index)

    output:
    tuple val(meta), path("*.clean.fastq.gz")   , emit: fastq
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ? "--aligner-args ${task.ext.args2}": ''
    prefix = task.ext.prefix ?: "${meta.id}"
    if (meta.single_end) {
        """
        export HOSTILE_CACHE_DIR=./${index}
        hostile clean \\
            --fastq1 ${fastq} \\
            --index ${index} \\
            --reorder \\
            --out-dir ${prefix}.hostile \\
            --threads ${task.cpus} \\
            ${args2} \\
            ${args}

        mv ${prefix}.hostile/*.fastq.gz .

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            hostile: \$(hostile --version)
        END_VERSIONS
        """
    } else {
        """
        export HOSTILE_CACHE_DIR=./${index}
        hostile clean \\
            --fastq1 ${fastq[0]} \\
            --fastq2 ${fastq[1]} \\
            --index ${index_name} \\
            --reorder \\
            --out-dir ${prefix}.hostile \\
            --threads ${task.cpus} \\
            ${args2} \\
            ${args}

        mv ${prefix}.hostile/*.fastq.gz .

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            hostile: \$(hostile --version)
        END_VERSIONS
        """
    }

    stub:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    if (meta.single_end) {
        """
        echo ${args}
        echo ${args2}

        echo "" | gzip > ${prefix}.clean.fastq.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            hostile: \$(hostile --version)
        END_VERSIONS
        """
    } else {
        """
        echo ${args}
        echo ${args2}

        echo "" | gzip > ${prefix}_1.clean.fastq.gz
        echo "" | gzip > ${prefix}_2.clean.fastq.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            hostile: \$(hostile --version)
        END_VERSIONS
        """
    }
}
