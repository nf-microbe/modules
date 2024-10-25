process MVIRS_PARSER {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/hmmer_prodigal-gv_biopython_pandas:de876e9113694da7' :
        'community.wave.seqera.io/library/hmmer_prodigal-gv_biopython_pandas:aaaf7b9a9207df90' }"

    input:
    tuple val(meta) , path(contigs)
    tuple val(meta2), path(prophages)
    path integrases_hmm

    output:
    tuple val(meta), path("${prefix}.mvirs.summary.tsv")    , emit: summary
    tuple val(meta), path("${prefix}.mvirs.integrases.tsv") , emit: integrases
    path "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    mvirs_parser.py \\
        --mvirs_path ${prophages} \\
        --fna_path ${contigs} \\
        --hmm_path ${integrases_hmm} \\
        --prefix ${prefix}.mvirs \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed 's/Python //' )
        biopython: \$(python -c "import Bio; print(Bio.__version__)")
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
        prodigal-gv: \$(echo \$(prodigal-gv -v 2>&1) | sed -n 's/^.*Prodigal V//; s/-gv.*//; 1p')
        hmmer: \$(hmmsearch -h | grep HMMER | sed 's/.*HMMER //; s/ (.*//')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.mvirs.summary.tsv
    touch ${prefix}.mvirs.integrases.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed 's/Python //' )
        biopython: \$(python -c "import Bio; print(Bio.__version__)")
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
        prodigal-gv: \$(echo \$(prodigal-gv -v 2>&1) | sed -n 's/^.*Prodigal V//; s/-gv.*//; 1p')
        hmmer: \$(hmmsearch -h | grep HMMER | sed 's/.*HMMER //; s/ (.*//')
    END_VERSIONS
    """
}
