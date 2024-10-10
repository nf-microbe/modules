// Import modules
include { PHABLES_INSTALL       } from '../../../modules/nf-core/phables/install/main'
include { PHABLES_RUN           } from '../../../modules/nf-core/phables/run/main'

workflow FASTQGFA_PHABLES_FASTA {
    take:
    fastq_gz            // [ [ meta.id, meta.single_end ], [ reads_1.fastq.gz, reads_2.fastq.gz ] ]
    gfa_gz              // [ [ meta.id, meta.assembler ], gfa.gz ]
    phables_config      // [ phables.config ]
    phables_db          // [ phables_db ]

    main:
    ch_versions = Channel.empty()

    // if genomad_db exists, skip PHABLES_INSTALL
    if (phables_db) {
        ch_phables_db = phables_db
    } else {
        //
        // MODULE: download phables database
        //
        PHABLES_INSTALL(
        )
        ch_phables_db   = PHABLES_INSTALL.out.phables_db
        ch_versions     = ch_versions.mix(PHABLES_INSTALL.out.versions)
    }

    // combine fastq and gfa for phables
    ch_phables_input = fastq_gz.map { meta, fastq -> [ meta.id, fastq ] }
        .combine(gfa_gz.map { meta, gfa -> [ meta.id, meta, gfa ] }, by:0)
        .multiMap { meta_id, fastq, meta_gfa, gfa ->
            fastq : [ meta_gfa, fastq ]
            gfa   : [ meta_gfa, gfa ]
        }

    //
    // MODULE: Run phables to extend assemblies
    //
    PHABLES_RUN(
        ch_phables_input.fastq,
        ch_phables_input.gfa,
        phables_config,
        ch_phables_db
    )
    ch_phables_fasta_gz = PHABLES_RUN.out.fasta
    ch_versions = ch_versions.mix(PHABLES_RUN.out.versions)

    emit:
    fasta_gz    = ch_phables_fasta_gz   // [ [ meta,  ], extended_contigs.fna.gz ]    , FASTA file containing extended contigs
    versions    = ch_versions.unique()  // [ versions.yml ]
}
