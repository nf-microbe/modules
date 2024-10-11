process SRA_ASPERACLI {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://carsonjm/asperacli_fastp_megahit_phables_seqkit:0.0.1' :
        'docker.io/carsonjm/asperacli_fastp_megahit_phables_seqkit:0.0.1' }"

    input:
    tuple val(meta), val(accession)

    output:
    tuple val(meta), path("${accession}*.fastq.gz") , emit: fastq
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args            = task.ext.args ?: ''
    def conda_prefix    = ['singularity', 'apptainer'].contains(workflow.containerEngine) ? "export CONDA_PREFIX=/opt/conda" : ""
    prefix              = task.ext.prefix ?: "${meta.id}"
    """
    ${conda_prefix}

    ### Download FastQ metadata
    curl -X 'GET' \\
        'https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${accession}&result=read_run&fields=library_layout,fastq_md5,fastq_aspera' \\
        -H 'accept: */*' > acc_metadata.tsv

    md5=\$(awk 'NR>1 { print \$2 }' acc_metadata.tsv)
    fastq=\$(awk 'NR>1 { print \$3 }' acc_metadata.tsv)
    layout=\$(awk 'NR>1 { print \$4 }' acc_metadata.tsv)

    ### Download FastQ files
    if [[ \$layout == "PAIRED" ]];
    then
        fastq_array=(\${fastq//;/ })
        md5_array=(\${md5//;/ })

        ascp \\
            -QT -l 300m -P 33001 \\
            -i \$CONDA_PREFIX/etc/aspera/aspera_bypass_dsa.pem \\
            era-fasp@\${fastq_array[0]} \\
            ${accession}_1.fastq.gz

        echo "\${md5_array[0]}  ${accession}_1.fastq.gz" > ${prefix}_1.fastq.gz.md5
        md5sum -c ${accession}_1.fastq.gz.md5

        ascp \\
            -QT -l 300m -P 33001 \\
            -i \$CONDA_PREFIX/etc/aspera/aspera_bypass_dsa.pem \\
            era-fasp@\${fastq_array[1]} \\
            ${accession}_2.fastq.gz

        echo "\${md5_array[1]}  ${accession}_2.fastq.gz" > ${prefix}_2.fastq.gz.md5
        md5sum -c ${accession}_2.fastq.gz.md5
    else
        ascp \\
            -QT -l 300m -P 33001 \\
            -i \$CONDA_PREFIX/etc/aspera/aspera_bypass_dsa.pem \\
            era-fasp@\$fastq \\
            ${accession}.fastq.gz

        echo "\$md5  ${accession}.fastq.gz" > ${prefix}.fastq.gz.md5
        md5sum -c ${accession}.fastq.gz.md5
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ascli: 4.14.0
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    ### Download FastQ metadata
    curl -X 'GET' \\
        'https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${accession}&result=read_run&fields=library_layout,fastq_md5,fastq_aspera' \\
        -H 'accept: */*' > acc_metadata.tsv

    md5=\$(awk 'NR>1 { print \$2 }' acc_metadata.tsv)
    fastq=\$(awk 'NR>1 { print \$3 }' acc_metadata.tsv)
    layout=\$(awk 'NR>1 { print \$4 }' acc_metadata.tsv)

    ### Touch empty FastQ files
    if [[ \$layout == "PAIRED" ]];
    then
        echo "" | gzip > ${accession}_1.fastq.gz
        echo "" | gzip > ${accession}_2.fastq.gz
    else
        echo "" | gzip > ${accession}.fastq.gz
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ascli: 4.14.0
    END_VERSIONS
    """
}
