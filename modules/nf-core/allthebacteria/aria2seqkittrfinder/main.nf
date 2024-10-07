process ALLTHEBACTERIA_ARIA2SEQKITTRFINDER {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/seqkit_aria2_biopython:2bb0d76d26291e6d' :
        'community.wave.seqera.io/library/seqkit_aria2_biopython:2256c902d3a46191' }"

    input:
    tuple val(meta), val(url)

    output:
    tuple val(meta), path("${prefix}.trfinder.fasta.gz")    , emit: fasta
    path "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def url_list    = url.collect { urls -> urls.toString() }
    def args        = task.ext.args ?: ''
    def args2       = task.ext.args2 ?: ''
    def args3       = task.ext.args3 ?: ''
    prefix          = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p tmp
    echo "${url_list.join('\n')}" > aria2_file.tsv
    ### Download AllTheBacteria assemblies
    aria2c \\
        --input-file=aria2_file.tsv \\
        --dir=tmp \\
        --max-connection-per-server=${task.cpus} \\
        --split=${task.cpus} \\
        --max-concurrent-downloads=${task.cpus} \\
        ${args}
    tar -xvf tmp/*.tar.xz -C tmp
    cat tmp/**/*.fa > ${prefix}.fasta
    rm -rf tmp
    ### Remove short contigs
    seqkit \\
        seq \\
        --threads ${task.cpus} \\
        ${args2} \\
        ${prefix}.fasta \\
        --out-file ${prefix}.seqkit.fasta
    rm ${prefix}.fasta
    ### Identify Terminal repeats
    trfinder.py \\
        --input ${prefix}.seqkit.fasta \\
        --prefix ${prefix} \\
        ${args3}
    gzip -f ${prefix}.trfinder.fasta
    rm -rf ${prefix}.seqkit.fasta
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        aria2: \$(echo \$(aria2c --version 2>&1) | grep 'aria2 version' | cut -f3 -d ' ')
        seqkit: \$(seqkit version | cut -d' ' -f2)
        python: \$( python --version | sed 's/Python //' )
        biopython: \$(python -c "import Bio; print(Bio.__version__)")
    END_VERSIONS
    """

    stub:
    def args   = task.ext.args ?: ''
    def args2   = task.ext.args2 ?: ''
    def args3   = task.ext.args3 ?: ''
    prefix  = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.trfinder.fasta.gz
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        aria2: \$(echo \$(aria2c --version 2>&1) | grep 'aria2 version' | cut -f3 -d ' ')
        seqkit: \$(seqkit version | cut -d' ' -f2)
        python: \$( python --version | sed 's/Python //' )
        biopython: \$(python -c "import Bio; print(Bio.__version__)")
    END_VERSIONS
    """
}
