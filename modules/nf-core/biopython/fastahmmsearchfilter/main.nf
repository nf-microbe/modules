process BIOPYTHON_FASTAHMMSEARCHFILTER {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/hmmer_prodigal-gv_biopython_pandas:de876e9113694da7' :
        'community.wave.seqera.io/library/hmmer_prodigal-gv_biopython_pandas:aaaf7b9a9207df90' }"

    input:
    tuple val(meta) , path(fasta)
    tuple val(meta2), path(hmm_tbl)

    output:
    tuple val(meta), path("${prefix}.hmmfilter.fasta.gz")   , emit: fasta
    path "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    fastahmmsearchfilter.py \\
        --input_fasta ${fasta} \\
        --hmm_tbl ${hmm_tbl} \\
        --prefix ${prefix}.hmmfilter \\
        ${args}

    gzip ${prefix}.hmmfilter.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed 's/Python //' )
        biopython: \$(python -c "import Bio; print(Bio.__version__)")
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.hmmfilter.fasta.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed 's/Python //' )
        biopython: \$(python -c "import Bio; print(Bio.__version__)")
    END_VERSIONS
    """
}
