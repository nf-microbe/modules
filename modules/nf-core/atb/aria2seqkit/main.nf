process ATB_ARIA2SEQKIT {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/seqkit_aria2_biopython:2bb0d76d26291e6d' :
        'community.wave.seqera.io/library/seqkit_aria2_biopython:2256c902d3a46191' }"

    input:
    tuple val(meta), val(sample_id)
    val atb_file_list

    output:
    tuple val(meta), path("${sample_id}.fasta.gz")  , emit: fasta
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args ?: ''
    def args2       = task.ext.args2 ?: ''
    prefix          = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p tmp/download tmp/seqkit
    ### Identify assembly to download based on sample_id
    atb_data=\$(zcat ${atb_file_list} | grep "${sample_id}" | cut -f4,5,6)
    read -a atb_array <<< "\$atb_data"

    ### Download AllTheBacteria assemblies
    aria2c \\
        \${atb_array[2]} \\
        --out=tmp/download/\${atb_array[1]} \\
        --max-connection-per-server=${task.cpus} \\
        --split=${task.cpus} \\
        --max-concurrent-downloads=${task.cpus} \\
        ${args}

    tar -xvf tmp/download/*.tar.xz -C tmp/download
    rm tmp/download/*.tar.xz

    ### Remove short contigs
    seqkit \\
        seq \\
        --threads ${task.cpus} \\
        ${args2} \\
        tmp/download/\${atb_array[0]} \\
        --out-file ${sample_id}.fasta.gz

    rm -rf tmp/download/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        aria2: \$(echo \$(aria2c --version 2>&1) | grep 'aria2 version' | cut -f3 -d ' ')
        seqkit: \$(seqkit version | cut -d' ' -f2)
    END_VERSIONS
    """

    stub:
    def args   = task.ext.args ?: ''
    def args2   = task.ext.args2 ?: ''
    prefix  = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${sample_id}.fasta.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        aria2: \$(echo \$(aria2c --version 2>&1) | grep 'aria2 version' | cut -f3 -d ' ')
        seqkit: \$(seqkit version | cut -d' ' -f2)
    END_VERSIONS
    """
}
