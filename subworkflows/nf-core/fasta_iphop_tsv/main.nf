// Import modules
include { IPHOP_DOWNLOAD    } from '../../../modules/nf-core/iphop/download'
include { IPHOP_PREDICT     } from '../../../modules/nf-core/iphop/predict'

workflow FASTA_IPHOP_TSV {
    take:
    fasta_gz    // [ [ meta ], fasta.gz ]
    iphop_db    // [ iphop_db ]

    main:
    ch_versions             = Channel.empty()

    // if iphop_db exists, skip IPHOP_DOWNLOAD
    if (iphop_db){
        ch_iphop_db = iphop_db
    } else {
        //
        // MODULE: download iPHOP database
        //
        IPHOP_DOWNLOAD(
        )
        ch_iphop_db = IPHOP_DOWNLOAD.out.iphop_db
        ch_versions = ch_versions.mix(IPHOP_DOWNLOAD.out.versions)
    }

    //
    // MODULE: Run geNomad end-to-end
    //
    IPHOP_PREDICT(
        fasta_gz,
        ch_iphop_db
    )
    ch_versions = ch_versions.mix(IPHOP_PREDICT.out.versions)

    emit:
    ch_iphop_genus_tsv  = IPHOP_PREDICT.out.iphop_genus     // [ [ meta ], Host_prediction_to_genus_m*.csv ]
    ch_iphop_genome_tsv = IPHOP_PREDICT.out.iphop_genome    // [ [ meta ], Host_prediction_to_genome_m*.csv ]
    versions            = ch_versions                       // channel: [ versions.yml ]
}
