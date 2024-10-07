process PYHMMER_PLASMIDCLASSIFY {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/pyhmmer_pandas:9f16bf144bcd5b1a':
        'community.wave.seqera.io/library/pyhmmer_pandas:842847931f77e2ab' }"

    input:
    tuple val(meta), path(faa)
    path plasmid_hmms

    output:
    tuple val(meta), path("${prefix}.plasmid_hmms.tsv")     , emit: hits
    tuple val(meta), path("${prefix}.plasmid_markers.tsv")  , emit: markers
    path "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    pyhmmer_plasmidclassify.py \\
        --input ${faa} \\
        --plasmid_hallmarks ${plasmid_hmms} \\
        --prefix ${prefix} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed 's/Python //' )
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
        pyhmmer: \$(python -c "import pyhmmer; print(pyhmmer.__version__)")
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.plasmid_hmms.tsv
    touch ${prefix}.plasmid_markers.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed 's/Python //' )
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
        pyhmmer: \$(python -c "import pyhmmer; print(pyhmmer.__version__)")
    END_VERSIONS
    """
}
