//
// Subworkflow with functionality that is required for nf-microbe subworkflows
//

/*
========================================================================================
    SUBWORKFLOW DEFINITION
========================================================================================
*/

workflow UTILS_NFMICROBE_FUNCTIONS {

    main:

    ch_dummy_channel = Channel.empty()

    emit:
    dummy_emit = true

}

/*
========================================================================================
    FUNCTIONS
========================================================================================
*/

//
// Identify work directories to clean
//
def getWorkDirs(ch_to_clean, ch_downstream, print_output) {
    // combine channel to clean and dependent channel to clean only channels that have an output
    def ch_double_workdir1 = Channel.empty()
    def ch_double_workdir2 = Channel.empty()
    def ch_branch = ch_to_clean.combine(ch_downstream, by:0)
        .branch { it ->
            double_ch: it.size() > 3
            single_ch: true
        }
    ch_double_workdir1 = ch_branch.double_ch
        .map { meta, ch_to_clean1, ch_to_clean2, ch_dep ->
            return [ [ meta ], ch_to_clean1, ch_dep ]
        }
    ch_double_workdir2 = ch_branch.double_ch
        .map { meta, ch_to_clean1, ch_to_clean2, ch_dep ->
            return [ [ meta ], ch_to_clean2, ch_dep ]
        }

    def ch_workdirs = ch_double_workdir1.mix(ch_double_workdir2).mix(ch_branch.single_ch)
        // filter to retain work directory
        .map { meta, files_to_clean, dependent_files ->
            // do not clean directory if it is not a work directory
            if (( files_to_clean =~ /(^.*\/work\/[^\/]+\/[^\/]+\/).*/ )) {
                def dir_to_clean = ( files_to_clean =~ /(^.*\/work\/[^\/]+\/[^\/]+\/).*/ )[0][1]
                    return [ [ id: meta.id ], dir_to_clean ]
            }
        }
    // remove redundancy
    .unique()

    // print output for tests, otherwise should be false
    if (print_output) {
        ch_workdirs.view()
    }

    return ch_workdirs
}

//
// Filter channels to remove empty FastQ files
//
def rmEmptyFastQs(ch_fastqs, print_output, stub_run=workflow.stubRun) {
    def ch_nonempty_fastqs = ch_fastqs
        .filter { meta, fastq ->
                if ( meta.single_end ) {
                    try {
                        fastq.countFastq(limit: 10) > 0
                    } catch (java.util.zip.ZipException e) {
                        log.debug "[rmEmptyFastQs]: ${fastq} is not in GZIP format, this is likely because it was cleaned with --remove_intermediate_files"
                        true
                    } catch (EOFException) {
                        log.debug "[rmEmptyFastQs]: ${fastq} has an EOFException, this is likely an empty gzipped file."
                    }
                } else {
                    try {
                        fastq[0].countFastq( limit: 10 ) > 0
                    } catch (java.util.zip.ZipException e) {
                        log.debug "[rmEmptyFastQs]: ${fastq} is not in GZIP format, this is likely because it was cleaned with --remove_intermediate_files"
                        true
                    } catch (EOFException) {
                        log.debug "[rmEmptyFastQs]: ${fastq[0]} has an EOFException, this is likely an empty gzipped file."
                    }
                }
            }

    if (stub_run) {
        ch_nonempty_fastqs = ch_fastqs
    }

    if (print_output) {
        ch_nonempty_fastqs.view()
    }

    return ch_nonempty_fastqs
}

//
// Filter channels to remove empty FastA files
//
def rmEmptyFastAs(ch_fastas, print_output, stub_run=workflow.stubRun) {
    def ch_nonempty_fastas = ch_fastas
        .filter { meta, fasta ->
            try {
                fasta.countFasta( limit: 10 ) > 0
            } catch (java.util.zip.ZipException e) {
                log.debug "[rmEmptyFastAs]: ${fasta} is not in GZIP format, this is likely because it was cleaned with --remove_intermediate_files"
                true
            } catch (EOFException) {
                log.debug "[rmEmptyFastAs]: ${fasta} has an EOFException, this is likely an empty gzipped file."
            }
        }

    if (stub_run) {
        ch_nonempty_fastas = ch_fastas
    }

    if (print_output) {
        ch_nonempty_fastas.view()
    }

    return ch_nonempty_fastas
}


