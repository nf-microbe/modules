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
    tuple val(meta), path("${prefix}.trfinder.tsv")         , emit: tsv
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
    mkdir -p tmp/download tmp/seqkit tmp/trfinder
    echo "${url_list.join('\n')}" > aria2_file.tsv

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

    #rm -rf tmp/download/

    ### Identify Terminal repeats
    cd tmp/trfinder

    for file in ../seqkit/*.fasta; do
        filename=\$(basename \$file)

        trfinder.py \\
            --input \$file \\
            --prefix ATB_\${filename%.*} \\
            ${args3}
    done

    cd ../..
    #rm -rf tmp/seqkit/

    # combine and compress tr finder files
    cat tmp/trfinder/*.trfinder.fasta > ${prefix}.trfinder.fasta
    awk '(NR == 1) || (FNR > 1)' tmp/trfinder/*.trfinder.tsv > ${prefix}.trfinder.tsv

    #rm -rf tmp/trfinder/
    gzip -f ${prefix}.trfinder.fasta

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
    touch ${prefix}.trfinder.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        aria2: \$(echo \$(aria2c --version 2>&1) | grep 'aria2 version' | cut -f3 -d ' ')
        seqkit: \$(seqkit version | cut -d' ' -f2)
        python: \$( python --version | sed 's/Python //' )
        biopython: \$(python -c "import Bio; print(Bio.__version__)")
    END_VERSIONS
    """
}
