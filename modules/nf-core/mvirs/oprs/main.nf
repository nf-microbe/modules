process MVIRS_OPRS {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mvirs:1.1.1--pyhdfd78af_0' :
        'biocontainers/mvirs:1.1.1--pyhdfd78af_0' }"

    input:
    tuple val(meta) , path(fastq)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(index)

    output:
    tuple val(meta), path("${prefix}.mvirs.bam")        , emit: bam
    tuple val(meta), path("${prefix}.mvirs.fasta.gz")   , emit: fasta
    tuple val(meta), path("${prefix}.mvirs.oprs")       , emit: oprs
    tuple val(meta), path("${prefix}.mvirs.clipped")    , emit: clipped
    tuple val(meta), path("${prefix}.mvirs.prophages")  , emit: prophage
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    mvirs \\
        oprs \\
        -f ${fastq[0]} \\
        -r ${fastq[1]} \\
        -db ${index[0].getBaseName()} \\
        -o ${prefix}.mvirs \\
        -t ${task.cpus} \\
        ${args} || if grep -q "Finishing mVIRs" .command.log; then true; else false; fi

    ### Extract fasta header information into TSV
    grep "^>" ${prefix}.mvirs.fasta \\
        | sed 's/^>//; s/\t/ /; s/:/ /; s/-/ /g; s/OPRs=//; s/HSs=//; s/SF=//' \\
        | awk '{print \$1"\\t"\$2"\\t"\$3"\\t"\$4"\\t"\$5"\\t"\$6}' > ${prefix}.mvirs.prophages  || true

    gzip -f ${prefix}.mvirs.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mvirs: \$(echo \$(mvirs 2>&1) | sed 's/.*Version: //; s/ .*//')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.mvirs.bam
    echo "" | gzip > ${prefix}.mvirs.fasta.gz
    touch ${prefix}.mvirs.oprs
    touch ${prefix}.mvirs.clipped
    touch ${prefix}.mvirs.prophages

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mvirs: \$(echo \$(mvirs 2>&1) | sed 's/.*Version: //; s/ .*//')
    END_VERSIONS
    """
}
