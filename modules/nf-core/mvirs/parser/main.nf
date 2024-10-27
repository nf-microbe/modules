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
    tuple val(meta3), path(hmmsearch)
    tuple val(meta4), path(contigs_faa)

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
        --mvirs ${prophages} \\
        --fna ${contigs} \\
        --hmmsearch ${hmmsearch} \\
        --faa ${contigs_faa} \\
        --prefix ${prefix}.mvirs \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed 's/Python //' )
        biopython: \$(python -c "import Bio; print(Bio.__version__)")
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
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
    END_VERSIONS
    """
}
