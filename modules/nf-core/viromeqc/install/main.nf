process VIROMEQC_INSTALL {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-b28a1a551d380ce8d57f9d83894ccb9559b44404:08a4cf815fcef0080ede5b3633202cbda8edf59b-0':
        'biocontainers/mulled-v2-b28a1a551d380ce8d57f9d83894ccb9559b44404:08a4cf815fcef0080ede5b3633202cbda8edf59b-0' }"

    output:
    path("viromeqc_db/")    , emit: viromeqc_db
    path "versions.yml"     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    viromeQC.py \\
        --install \\
        --index_dir viromeqc_db \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        viromeqc: \$(echo \$(viromeQC.py --version 2>&1 | sed -n 's/^.*Version:\\s//; 1p' ))
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    """
    mkdir -p viromeqc_db
    touch viromeqc_db/SILVA_132_LSURef_tax_silva.clean.1.bt2
    touch viromeqc_db/SILVA_132_LSURef_tax_silva.clean.2.bt2
    touch viromeqc_db/SILVA_132_LSURef_tax_silva.clean.3.bt2
    touch viromeqc_db/SILVA_132_LSURef_tax_silva.clean.4.bt2
    touch viromeqc_db/SILVA_132_LSURef_tax_silva.clean.rev.1.bt2
    touch viromeqc_db/SILVA_132_LSURef_tax_silva.clean.rev.2.bt2
    touch viromeqc_db/SILVA_132_SSURef_Nr99_tax_silva.clean.1.bt2
    touch viromeqc_db/SILVA_132_SSURef_Nr99_tax_silva.clean.2.bt2
    touch viromeqc_db/SILVA_132_SSURef_Nr99_tax_silva.clean.3.bt2
    touch viromeqc_db/SILVA_132_SSURef_Nr99_tax_silva.clean.4.bt2
    touch viromeqc_db/SILVA_132_SSURef_Nr99_tax_silva.clean.rev.1.bt2
    touch viromeqc_db/SILVA_132_SSURef_Nr99_tax_silva.clean.rev.2.bt2
    touch viromeqc_db/amphora_bacteria.dmnd
    touch viromeqc_db/amphora_bacteria_294.dmnd

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        viromeqc: \$(echo \$(viromeQC.py --version 2>&1 | sed -n 's/^.*Version:\\s//; 1p' ))
    END_VERSIONS
    """
}
