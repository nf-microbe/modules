process ALLTHEBACTERIA_ARIA2SEQKIT {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/seqkit_aria2_biopython:2bb0d76d26291e6d' :
        'community.wave.seqera.io/library/seqkit_aria2_biopython:2256c902d3a46191' }"

    input:
    tuple val(meta), val(url)

    output:
    tuple val(meta), path("tmp/seqkit/*.fasta.gz")  , emit: fastas
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def url_list    = url.collect { urls -> urls.toString() }.join(',')
    def args        = task.ext.args ?: ''
    def args2       = task.ext.args2 ?: ''
    prefix          = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p tmp/download tmp/seqkit
    IFS=',' read -r -a url_array <<< "${url_list}"
    printf '%s\\n' "\${url_array[@]}" > aria2_file.tsv

    ### Download AllTheBacteria assemblies
    aria2c \\
        --input-file=aria2_file.tsv \\
        --dir=tmp/download/ \\
        --max-connection-per-server=${task.cpus} \\
        --split=${task.cpus} \\
        --max-concurrent-downloads=${task.cpus} \\
        ${args}

    tar -xvf tmp/download/*.tar.xz -C tmp/download
    rm tmp/download/*.tar.xz

    ### Remove short contigs
    for dir in tmp/download/*; do
        for file in \$dir/*; do
            filename=\$(basename \$file)

            seqkit \\
                seq \\
                --threads ${task.cpus} \\
                ${args2} \\
                \$file \\
                --out-file tmp/seqkit/\${filename%.*}.fasta
        done
    done

    rm -rf tmp/download/

    #### compress fasta files
    gzip tmp/seqkit/*.fasta

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
    mkdir -p tmp/seqkit/
    echo "" | gzip > tmp/seqkit/test1.fasta.gz
    echo "" | gzip > tmp/seqkit/test2.fasta.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        aria2: \$(echo \$(aria2c --version 2>&1) | grep 'aria2 version' | cut -f3 -d ' ')
        seqkit: \$(seqkit version | cut -d' ' -f2)
    END_VERSIONS
    """
}
