process ENA_ASPERACLISEQKITTRFINDER {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/aspera-cli_seqkit_biopython:b1c05ee4816a116b' :
        'community.wave.seqera.io/library/aspera-cli_seqkit_biopython:ce3fa7446f366f82' }"

    input:
    tuple val(meta), val(url)

    output:
    tuple val(meta), path("${prefix}.trfinder.fasta.gz")    , emit: fasta
    tuple val(meta), path("${prefix}.trfinder.tsv")         , emit: tsv
    path "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def url_list    = url.collect { urls -> urls.toString() }.join(',')
    def args        = task.ext.args ?: ''
    def args2       = task.ext.args2 ?: ''
    def conda_prefix= ['singularity', 'apptainer'].contains(workflow.containerEngine) ? "export CONDA_PREFIX=/opt/conda" : ""
    prefix          = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p tmp/download tmp/seqkit tmp/trfinder
    ${conda_prefix}
    IFS=',' read -r -a url_array <<< "${url_list}"

    ### Download ENA assemblies
    printf '%s\\n' "\${url_array[@]}" | xargs -I{} -n 1 -P ${task.cpus} bash -c \\
    'IFS=":" read -ra file <<< "{}"; \\
    mod_file=\$(echo \${file[1]:1} | tr / _); \\
    echo \${mod_file}; \\
    ascp \\
        -QT -l 300m -P 33001 \\
        -i \$CONDA_PREFIX/etc/aspera/aspera_bypass_dsa.pem \\
        era-fasp@{} \\
        tmp/download/\$mod_file.fna.gz'

    ### Remove short contigs
    for file in tmp/download/*; do
        filename=\$(basename \$file)

        seqkit \\
            seq \\
            --threads ${task.cpus} \\
            ${args} \\
            \$file \\
            --out-file tmp/seqkit/\${filename%%.*}.fasta
    done

    rm -rf tmp/download/

    ### Identify Terminal repeats
    cd tmp/trfinder

    for file in ../seqkit/*.fasta; do
        filename=\$(basename \$file)

        trfinder.py \\
            --input \$file \\
            --prefix ENA_\${filename%.*} \\
            ${args2}
    done

    cd ../..
    rm -rf tmp/seqkit/

    # combine and compress tr finder files
    cat tmp/trfinder/*.trfinder.fasta > ${prefix}.trfinder.fasta
    awk '(NR == 1) || (FNR > 1)' tmp/trfinder/*.trfinder.tsv > ${prefix}.trfinder.tsv

    rm -rf tmp/trfinder/
    gzip -f ${prefix}.trfinder.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ascli: 4.14.0
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
        ascli: 4.14.0
        seqkit: \$(seqkit version | cut -d' ' -f2)
        python: \$( python --version | sed 's/Python //' )
        biopython: \$(python -c "import Bio; print(Bio.__version__)")
    END_VERSIONS
    """
}
