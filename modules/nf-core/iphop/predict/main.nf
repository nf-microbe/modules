process IPHOP_PREDICT {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/iphop:1.3.3--fa72498c1bff3fe4':
        'biocontainers/iphop:1.3.3--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)
    path iphop_db

    output:
    tuple val(meta), path("${prefix}_Host_prediction_to_genus.csv")     , emit: iphop_genus
    tuple val(meta), path("${prefix}_Host_prediction_to_genome.csv")    , emit: iphop_genome
    tuple val(meta), path("${prefix}_Detailed_output_by_tool.csv")      , emit: iphop_detailed_output
    path "versions.yml"                                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def is_compressed = fasta.getName().endsWith(".gz") ? true : false
    fasta_name = fasta.getName().replace(".gz", "")
    """
    if [ "$is_compressed" == "true" ]; then
        gzip -c -d ${fasta} > ${fasta_name}
    fi

    iphop \\
        predict \\
        --fa_file ${fasta_name} \\
        --out_dir iphop_results \\
        --db_dir ${iphop_db} \\
        --num_threads ${task.cpus} \\
        $args

    mv iphop_results/Host_prediction_to_genus_m*.csv ./${prefix}_Host_prediction_to_genus.csv
    mv iphop_results/Host_prediction_to_genome_m*.csv ./${prefix}_Host_prediction_to_genome.csv
    mv iphop_results/Detailed_output_by_tool.csv ./${prefix}_Detailed_output_by_tool.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        iphop: \$(echo \$(iphop --version 2>&1) | head -n 1 | sed 's/^.*iPHoP v//; s/: integrating.*\$//' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def min_score = args.contains('--min_score') ? args.split('--min_score ')[1] : '90'
    """
    touch ./${prefix}_Host_prediction_to_genus.csv
    touch ./${prefix}_Host_prediction_to_genome.csv
    touch ./${prefix}_Detailed_output_by_tool.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        iphop: \$(echo \$(iphop --version 2>&1) | head -n 1 | sed 's/^.*iPHoP v//; s/: integrating.*\$//' )
    END_VERSIONS
    """
}
