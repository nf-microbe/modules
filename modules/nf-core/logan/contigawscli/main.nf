process LOGAN_CONTIGAWSCLI {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/seqkit_awscli_biopython:4bbd3cad980280ad' :
        'community.wave.seqera.io/library/seqkit_awscli_biopython:1ebd36d34b5de691' }"

    input:
    tuple val(meta), val(accession)

    output:
    tuple val(meta), path("${accession}.contigs.fa.gz") , emit: fasta
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    prefix  = task.ext.prefix ?: "${meta.id}"
    """
    ### Download Logan contigs
    aws s3 cp s3://logan-pub/c/${accession}/${accession}.contigs.fa.zst . --no-sign-request || true

    zstd -d ${accession}.contigs.fa.zst

    gzip ${accession}.contigs.fa
    rm -rf ${accession}.contigs.fa.zst

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awscli: \$( aws --version | sed 's/aws-cli\\///; s/ Python.*//' )
        seqkit: \$(seqkit version | cut -d' ' -f2)
        python: \$( python --version | sed 's/Python //' )
        biopython: \$(python -c "import Bio; print(Bio.__version__)")
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${accession}.contigs.fa.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awscli: \$( aws --version | sed 's/aws-cli\\///; s/ Python.*//' )
        seqkit: \$(seqkit version | cut -d' ' -f2)
        python: \$( python --version | sed 's/Python //' )
        biopython: \$(python -c "import Bio; print(Bio.__version__)")
    END_VERSIONS
    """
}
