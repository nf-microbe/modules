process PROPAGATE {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vrhyme:1.1.0--pyhdfd78af_1' :
        'biocontainers/vrhyme:1.1.0--pyhdfd78af_1' }"

    input:
    tuple val(meta) , path(fastq)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(coords)

    output:
    tuple val(meta), path("${prefix}.propagate.tsv")    , emit: results
    tuple val(meta), path("${prefix}.propagate.log")    , emit: log     , optional: true
    path  "versions.yml"                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def is_compressed = fasta.getExtension() == "gz" ? true : false
    def fasta_name = is_compressed ? fasta.getBaseName() : fasta
    """
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${fasta} > ${fasta_name}
    fi

    Propagate \\
        -f ${fasta_name} \\
        -v ${coords} \\
        -r ${fastq} \\
        -o ${prefix}.propagate \\
        -t ${task.cpus} \\
        ${args}

    mv ${prefix}.propagate/${prefix}.propagate.tsv ./${prefix}.propagate.tsv
    mv ${prefix}.propagate/${prefix}.propagate.log ./${prefix}.propagate.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        propagate: \$(Propagate --version 2>&1 | sed 's/^.*PropagAtE v//; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.propagate.tsv
    touch ${prefix}.propagate.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        propagate: \$(Propagate --version 2>&1 | sed 's/^.*PropagAtE v//; s/ .*\$//')
    END_VERSIONS
    """
}
