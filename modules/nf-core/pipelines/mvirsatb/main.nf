process PIPELINES_MVIRSATB {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/fastp_hmmer_mvirs_sra-tools_pruned:e51b03139b970e16' :
        'community.wave.seqera.io/library/fastp_hmmer_mvirs_sra-tools_pruned:0c067f9fe60f5eb1' }"

    input:
    tuple val(meta) , path(fasta)
    tuple val(meta2), path(faa)
    tuple val(meta3), val(run)

    output:
    tuple val(meta), path("${prefix}.tbl")                  , emit: hmmer_tbl
    tuple val(meta), path("${prefix}.mvirs.summary.tsv")    , emit: summary     , optional: true
    tuple val(meta), path("${prefix}.mvirs.integrases.tsv") , emit: integrases  , optional: true
    path "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix  = task.ext.prefix ?: "${meta.id}"
    """
    #--------------------------------
    # ASSEMBLY FILTERING
    #--------------------------------
    # Length filtering
    echo "Step 1: Length filtering FastA & FAA files"

    fastafaalengthfilter.py \\
        --input_fasta ${fasta} \\
        --input_faa ${faa} \\
        --prefix ${prefix}.lengthfilter \\
        --fasta_min_len ${params.assembly_min_len}

    # HMMER to identify integrases if FastA not empty
    if [ \$(cat ${prefix}.lengthfilter.fasta | wc -l) -gt 1 ]; then
        echo "Step 2: Running HMMER"

        hmmsearch \\
            ${params.hmmsearch_args} \\
            --cpu ${task.cpus} \\
            --tblout ${prefix}.tbl \\
            ${params.hmm_file} \\
            ${prefix}.lengthfilter.faa
    else
        echo "Pipeline stopped: Length filtered FastA empty"
        touch ${prefix}.tbl
    fi

    # HMM filtering if HMM hit detected
    if [ \$(cat ${prefix}.tbl | wc -l) -gt 13 ]; then
        echo "Step 3: Filtering FastA to only contigs with HMM hits"

        fastahmmsearchfilter.py \\
            --input_fasta ${prefix}.lengthfilter.fasta \\
            --hmm_tbl ${prefix}.tbl \\
            --prefix ${prefix}.hmmfilter

        # REMOVE LENGTH FILTERED FASTA
        rm ${prefix}.lengthfilter.fasta
    else
        echo "Pipeline stopped: No HMM hits detected"
    fi

    #--------------------------------
    # READ DOWNLOAD
    #--------------------------------
    # Download FastQs if HMM hit detected
    if [ -f ${prefix}.hmmfilter.fasta ]; then
        echo "Step 4: Downloading FastQ files associated with filtered FastA."

        prefetch \\
            ${params.prefetch_args} \\
            ${run}

        fasterq-dump \\
            ${params.fasterq_dump_args} \\
            --threads ${task.cpus} \\
            --outfile ${run}.fastq \\
            ${run}

        # REMOVE SRA FILE
        rm ${run}/${run}.sra
    fi

    #--------------------------------
    # READ QC
    #--------------------------------
    # QC reads with fastp
    if [ -f ${run}_1.fastq ] && [ -f ${run}_2.fastq ]; then
        echo "Step 5: QC'ing FastQ files."

        fastp \\
            --in1 ${run}_1.fastq \\
            --in2 ${run}_2.fastq \\
            --out1 ${prefix}_1.fastp.fastq \\
            --out2 ${prefix}_2.fastp.fastq \\
            --thread ${task.cpus} \\
            --detect_adapter_for_pe \\
            ${params.fastp_args}

        # REMOVE SRA-TOOLS FASTQs
        rm ${run}_1.fastq ${run}_2.fastq

        # REMOVE FAST HTML
        rm fastp.html
    else
        echo "Pipeline stopped: PE FastQ files not present."
    fi

    #--------------------------------
    # MVIRS
    #--------------------------------
    if [ -f ${prefix}_1.fastp.fastq ]; then
        echo "Step 6: Running mVIRs."

        # Create an index for mvirs
        mvirs \\
            index \\
            -f ${prefix}.hmmfilter.fasta

        # Run mvirs to identify split reads/OPRs
        mvirs \\
            oprs \\
            -f ${prefix}_1.fastp.fastq \\
            -r ${prefix}_2.fastp.fastq \\
            -db ${prefix}.hmmfilter.fasta \\
            -o ${prefix}.mvirs \\
            -t ${task.cpus} \\
            ${params.mvirs_oprs_args}

        # REMOVE FASTP FASTQs
        rm ${prefix}_1.fastp.fastq ${prefix}_2.fastp.fastq

        # REMOVE MVIRS BAMS
        rm ${prefix}.mvirs.bam

        # REMOVE MVIRS INDICES
        rm ${prefix}.hmmfilter.fasta.*
    else
        echo "Pipeline stopped: PE FastQ files not present."
    fi

    if [ \$(cat ${prefix}.mvirs.fasta | wc -l) -gt 1 ]; then
        echo "Step 7: Parsing mVIRs output."

        # Parse mvirs outputs
        mvirs_parser.py \\
            --mvirs ${prefix}.mvirs.fasta \\
            --fna ${prefix}.hmmfilter.fasta \\
            --hmmsearch ${prefix}.tbl \\
            --faa ${prefix}.lengthfilter.faa \\
            --prefix ${prefix}.mvirs \\
            ${params.mvirs_parser_args}
    else
        echo "Pipeline stopped: No mVIRs MGEs detected."
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed 's/Python //' )
        biopython: \$(python -c "import Bio; print(Bio.__version__)")
        hmmer: \$(hmmsearch -h | grep -o '^# HMMER [0-9.]*' | sed 's/^# HMMER *//')
        sratools: \$(prefetch --version 2>&1 | grep -Eo '[0-9.]+')
        fastp: \$(fastp --version 2>&1 | sed -e "s/fastp //g")
        mvirs: \$(echo \$(mvirs 2>&1) | sed 's/.*Version: //; s/ .*//')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """

    stub:
    def args   = task.ext.args ?: ''
    def args2   = task.ext.args2 ?: ''
    def args3   = task.ext.args3 ?: ''
    prefix  = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tbl
    touch ${prefix}.mvirs.summary.tsv
    touch ${prefix}.mvirs.integrases.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed 's/Python //' )
        biopython: \$(python -c "import Bio; print(Bio.__version__)")
        hmmer: \$(hmmsearch -h | grep -o '^# HMMER [0-9.]*' | sed 's/^# HMMER *//')
        sratools: \$(prefetch --version 2>&1 | grep -Eo '[0-9.]+')
        fastp: \$(fastp --version 2>&1 | sed -e "s/fastp //g")
        mvirs: \$(echo \$(mvirs 2>&1) | sed 's/.*Version: //; s/ .*//')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """
}
