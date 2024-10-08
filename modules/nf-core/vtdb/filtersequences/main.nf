process VTDB_FILTERSEQUENCES {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-80c23cbcd32e2891421c54d1899665046feb07ef:77a31e289d22068839533bf21f8c4248ad274b60-0' :
        'biocontainers/mulled-v2-80c23cbcd32e2891421c54d1899665046feb07ef:77a31e289d22068839533bf21f8c4248ad274b60-0' }"

    input:
    tuple val(meta) , path(fasta)
    tuple val(meta2), path(faa)
    tuple val(meta5), path(completeness_data)
    path(sequences_to_keep)

    output:
    tuple val(meta), path("${prefix}_filtered.fasta.gz")        , emit: fasta
    tuple val(meta), path("${prefix}_filtered_proteins.faa.gz") , emit: proteins
    path "versions.yml"                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def seqs_to_keep = sequences_to_keep ? "--seqs_to_keep ${sequences_to_keep}" : '--seqs_to_keep NO_SEQS_TO_KEEP_FILE'
    def faa_args = faa ? "--input_proteins ${faa}" : '--input_proteins NO_PROTEINS_FILE'
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    filtersequences.py \\
        --input_fasta ${fasta} \\
        ${faa_args} \\
        --completeness_data ${completeness_data} \\
        ${seqs_to_keep} \\
        --output_fasta ${prefix}_filtered.fasta \\
        --output_proteins ${prefix}_filtered_proteins.faa \\
        ${args}

    gzip ${prefix}_filtered.fasta ${prefix}_filtered_proteins.faa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed 's/Python //' )
        pandas: \$(echo \$(pandas_version.py 2>&1))
        biopython: \$(echo \$(biopython_version.py 2>&1))
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}_filtered.fasta.gz
    echo "" | gzip > ${prefix}_filtered_proteins.faa.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed 's/Python //' )
        pandas: \$(echo \$(pandas_version.py 2>&1))
        biopython: \$(echo \$(biopython_version.py 2>&1))
    END_VERSIONS
    """
}
