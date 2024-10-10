process SRA_SRATOOLS {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/fastp_sra-tools:ce7baca742fa975e' :
        'community.wave.seqera.io/library/fastp_sra-tools:2234a1a896a58751' }"

    input:
    tuple val(meta), val(accession)

    output:
    tuple val(meta), path("${accession}*.fastq.gz") , emit: fastq
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def args2   = task.ext.args2 ?: ''
    prefix      = task.ext.prefix ?: "${meta.id}"
    """
    ### Download FastQ files with sratools
    prefetch \\
        ${accession} \\
        --output-file ${accession}/${accession}.sra \\
        ${args}

    fasterq-dump \\
        ${accession} \\
        --threads ${task.cpus} \\
        ${args2} \\


    rm -rf ${accession}/${accession}.sra
    gzip ${accession}*.fastq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sratools: \$(prefetch --version 2>&1 | grep -Eo '[0-9.]+')
    END_VERSIONS
    """

    stub:
    def args    = task.ext.args ?: ''
    prefix      = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${accession}.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sratools: \$(prefetch --version 2>&1 | grep -Eo '[0-9.]+')
    END_VERSIONS
    """
}
