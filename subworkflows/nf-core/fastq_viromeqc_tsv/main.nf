include { VIROMEQC_INSTALL  } from '../../../modules/nf-core/viromeqc/install/main'
include { VIROMEQC_VIROMEQC } from '../../../modules/nf-core/viromeqc/viromeqc/main'

workflow FASTQ_VIROMEQC_TSV {

    take:
    fastq_gz    // [ [ meta ]   , [ reads_1.fastq.gz, reads_2.fastq.gz ]
    viromeqc_db // [ viromeqc_db ]

    main:
    ch_versions = Channel.empty()

    // if viromeqc_db exists, skip VIROMEQC_INSTALL
    if (viromeqc_db){
        ch_genomad_db = viromeqc_db
    } else {
        //
        // MODULE: download viromeQC database
        //
        VIROMEQC_INSTALL(
        )
        ch_viromeqc_db  = VIROMEQC_INSTALL.out.viromeqc_db
        ch_versions     = ch_versions.mix(VIROMEQC_INSTALL.out.versions)
    }

    //
    // MODULE: Estimate viral enrichment
    //
    VIROMEQC_VIROMEQC(
        fastq_gz,
        ch_viromeqc_db
    )
    ch_versions         = ch_versions.mix(VIROMEQC_VIROMEQC.out.versions)

    emit:
    enrichment_tsv  = VIROMEQC_VIROMEQC.out.tsv // channel: [ [ meta ], TSV ]
    versions        = ch_versions                      // channel: [ versions.yml ]
}

