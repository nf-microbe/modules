process LOGAN_UNITIGAWSCLIMULTIPLIER {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/seqkit_awscli_biopython:4bbd3cad980280ad' :
        'community.wave.seqera.io/library/seqkit_awscli_biopython:1ebd36d34b5de691' }"

    input:
    tuple val(meta), val(accession)

    output:
    tuple val(meta), path("${accession}.unitigs.fa.gz")         , emit: raw_fasta
    tuple val(meta), path("${accession}.filter.fasta.gz")       , emit: filtered_fasta      , optional: true
    tuple val(meta), path("${accession}.multiplier.fasta.gz")   , emit: multiplied_fasta    , optional: true
    tuple val(meta), path("${accession}.filt_mult.fasta.gz")    , emit: filt_mult_fasta     , optional: true
    path "versions.yml"                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args            = task.ext.args ?: ''
    def run_filt        = task.ext.args.contains("-a")
    def run_mult        = task.ext.args.contains("-m")
    def run_filt_mult   = task.ext.args.contains("-a") && task.ext.args.contains("-m")
    prefix              = task.ext.prefix ?: "${meta.id}"
    """
    ### Download Logan unitigs
    aws s3 cp s3://logan-pub/u/${accession}/${accession}.unitigs.fa.zst . --no-sign-request || true

    ### Multiply Logan unitigs based on average kmer abundance
    zstd -d ${accession}.unitigs.fa.zst
    rm ${accession}.unitigs.fa.zst

    multiplier.py \\
        -i ${accession}.unitigs.fa \\
        -o ${accession} \\
        ${args}

    gzip ${accession}.unitigs.fa

    if [ "${run_filt}" == "true" ]; then
        gzip ${accession}.filter.fasta
    fi

    if [ "${run_mult}" == "true" ]; then
        gzip ${accession}.multiplier.fasta
    fi

    if [ "${run_filt_mult}" == "true" ]; then
        gzip ${accession}.filt_mult.fasta
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awscli: \$( aws --version | sed 's/aws-cli\\///; s/ Python.*//' )
        seqkit: \$(seqkit version | cut -d' ' -f2)
        python: \$( python --version | sed 's/Python //' )
        biopython: \$(python -c "import Bio; print(Bio.__version__)")
    END_VERSIONS
    """

    stub:
    def args            = task.ext.args ?: ''
    def run_filt        = task.ext.args.contains("-a")
    def run_mult        = task.ext.args.contains("-m")
    def run_filt_mult   = task.ext.args.contains("-a") && task.ext.args.contains("-m")
    prefix              = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${accession}.unitigs.fa.gz

    if [ "${run_filt}" == "true" ]; then
        echo "" | gzip > ${accession}.filtered.fasta.gz
    fi

    if [ "${run_mult}" == "true" ]; then
        echo "" | gzip > ${accession}.multiplier.fasta.gz
    fi

    if [ "${run_filt_mult}" == "true" ]; then
        echo "" | gzip > ${accession}.filt_mult.fasta.gz
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awscli: \$( aws --version | sed 's/aws-cli\\///; s/ Python.*//' )
        seqkit: \$(seqkit version | cut -d' ' -f2)
        python: \$( python --version | sed 's/Python //' )
        biopython: \$(python -c "import Bio; print(Bio.__version__)")
    END_VERSIONS
    """
}
