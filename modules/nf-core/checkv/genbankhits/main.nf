process CHECKV_GENBANKHITS {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/pyhmmer_pandas:9f16bf144bcd5b1a':
        'community.wave.seqera.io/library/pyhmmer_pandas:842847931f77e2ab' }"

    input:
    tuple val(meta) , path(aai)
    path ncbi
    path db

    output:
    tuple val(meta), path ("${prefix}_checkv_genbank_hits.tsv") , emit: genbank_hits
    path "versions.yml"                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    checkv_genbank_hits.py \\
        --checkv_aai ${aai} \\
        --checkv_db ${db} \\
        --ncbi_tsv ${ncbi} \\
        --output ${prefix}_checkv_genbank_hits.tsv \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed 's/Python //' )
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_checkv_genbank_hits.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed 's/Python //' )
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """
}
