process SEQKIT_SAMPLE {
    tag "${meta.id}"
    label 'process_low'
    // File IO can be a bottleneck. See: https://bioinf.shenwei.me/seqkit/usage/#parallelization-of-cpu-intensive-jobs

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqkit:2.8.1--h9ee0642_0':
        'biocontainers/seqkit:2.8.1--h9ee0642_0' }"

    input:
    tuple val(meta), path(fastx)

    output:
    tuple val(meta), path("${prefix}*") , emit: fastx
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args ?: ''
    def args2       = task.ext.args2 ?: ''
    prefix          = task.ext.prefix ?: "${meta.id}"
    def extension   = "fastq"
    if ("$fastx" ==~ /.+\.fasta|.+\.fasta.gz|.+\.fa|.+\.fa.gz|.+\.fas|.+\.fas.gz|.+\.fna|.+\.fna.gz|.+\.fsa|.+\.fsa.gz/ ) {
        extension   = "fasta"
    }
    extension       = fastx.toString().endsWith('.gz') ? "${extension}.gz" : extension
    def call_gzip   = extension.endsWith('.gz') ? "| gzip -c $args2" : ''
    if("${prefix}.${extension}" == "$fastx") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    if (meta.single_end || meta.assembler != null) {
        """
        seqkit \\
            sample \\
            --threads ${task.cpus} \\
            ${args} \\
            ${fastx} \\
            ${call_gzip} \\
            > ${prefix}.${extension}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            seqkit: \$(seqkit version | cut -d' ' -f2)
        END_VERSIONS
        """
    } else {
        """
        seqkit \\
            sample \\
            --threads ${task.cpus} \\
            ${args} \\
            ${fastx[0]} \\
            --rand-seed 1 \\
            ${call_gzip} \\
            > ${prefix}_1.${extension}

        seqkit \\
            sample \\
            --threads ${task.cpus} \\
            ${args} \\
            ${fastx[1]} \\
            --rand-seed 1 \\
            ${call_gzip} \\
            > ${prefix}_2.${extension}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            seqkit: \$(seqkit version | cut -d' ' -f2)
        END_VERSIONS
        """
    }

    stub:
    prefix          = task.ext.prefix ?: "${meta.id}"
    def extension   = "fastq"
    if ("$fastx" ==~ /.+\.fasta|.+\.fasta.gz|.+\.fa|.+\.fa.gz|.+\.fas|.+\.fas.gz|.+\.fna|.+\.fna.gz|.+\.fsa|.+\.fsa.gz/ ) {
        extension   = "fasta"
    }
    extension = fastx.toString().endsWith('.gz') ? "${extension}.gz" : extension
    if("${prefix}.${extension}" == "$fastx") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    def is_compressed = extension.endsWith(".gz") ? true : false
    if (meta.single_end || meta.assembler != null) {
        """
        if [ "$is_compressed" == "true" ]; then
            echo "" | gzip > ${prefix}.${extension}
        else
            touch ${prefix}.${extension}
        fi

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            seqkit: \$(seqkit version | cut -d' ' -f2)
        END_VERSIONS
        """
    } else {
        """
        if [ "$is_compressed" == "true" ]; then
            echo "" | gzip > ${prefix}_1.${extension}
            echo "" | gzip > ${prefix}_2.${extension}
        else
            touch ${prefix}_1.${extension}
            touch ${prefix}_2.${extension}
        fi

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            seqkit: \$(seqkit version | cut -d' ' -f2)
        END_VERSIONS
        """
    }
}
