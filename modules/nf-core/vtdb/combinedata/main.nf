process VTDB_COMBINEDATA {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-80c23cbcd32e2891421c54d1899665046feb07ef:77a31e289d22068839533bf21f8c4248ad274b60-0' :
        'biocontainers/mulled-v2-80c23cbcd32e2891421c54d1899665046feb07ef:77a31e289d22068839533bf21f8c4248ad274b60-0' }"

    input:
    tuple val(meta)     , path(fasta)
    tuple val(meta2)    , path(trfinder_stats)
    tuple val(meta3)    , path(genomad_scores)
    tuple val(meta4)    , path(genomad_genes)
    tuple val(meta5)    , path(genomad_taxa)
    tuple val(meta6)    , path(busco_hmms)
    tuple val(meta7)    , path(plasmid_hmms)
    tuple val(meta8)    , path(virus_hmms)
    tuple val(meta7)    , path(completeness)
    tuple val(meta8)    , path(contamination)
    tuple val(meta9)    , path(tantan)
    tuple val(meta10)   , path(nuc_stats)

    output:
    tuple val(meta), path("${prefix}_combined_data.tsv")    , emit: tsv
    path "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args                    = task.ext.args ?: ''
    def trfinder_input          = trfinder_stats ? "--trfinder ${trfinder_stats}" : ""
    def genomad_scores_input    = genomad_scores ? "--genomad_scores ${genomad_scores}" : ""
    def genomad_genes_input     = genomad_genes ? "--genomad_genes ${genomad_genes}" : ""
    def genomad_taxa_input      = genomad_taxa ? "--genomad_taxa ${genomad_taxa}" : ""
    def busco_hmms_input        = busco_hmms ? "--busco_hmms ${busco_hmms}" : ""
    def plasmid_hmms_input      = plasmid_hmms ? "--plasmid_hmms ${plasmid_hmms}" : ""
    def virus_hmms_input        = virus_hmms ? "--virus_hmms ${virus_hmms}" : ""
    def completeness_input      = completeness ? "--completeness ${completeness}" : ""
    def contamination_input     = contamination ? "--contamination ${contamination}" : ""
    def tantan_input            = tantan ? "--tantan ${tantan}" : ""
    def nucleotide_stats_input  = nuc_stats ? "--sequence_stats ${nuc_stats}" : ""
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    combinedata.py \\
        --input ${fasta} \\
        ${trfinder_input} \\
        ${genomad_scores_input} \\
        ${genomad_genes_input} \\
        ${genomad_taxa_input} \\
        ${busco_hmms_input} \\
        ${plasmid_hmms_input} \\
        ${virus_hmms_input} \\
        ${completeness_input} \\
        ${contamination_input} \\
        ${tantan_input} \\
        ${nucleotide_stats_input} \\
        --output ${prefix}_combined_data.tsv \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed 's/Python //' )
        biopython: \$(python -c "import Bio; print(Bio.__version__)")
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_combined_data.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed 's/Python //' )
        biopython: \$(python -c "import Bio; print(Bio.__version__)")
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """
}
