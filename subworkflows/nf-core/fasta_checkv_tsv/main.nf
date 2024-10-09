// Import modules
include { CHECKV_DOWNLOADDATABASE   } from '../../../modules/nf-core/checkv/downloaddatabase'
include { CHECKV_ENDTOEND           } from '../../../modules/nf-core/checkv/endtoend'
include { CHECKV_GENBANKHITS        } from '../../../modules/nf-core/checkv/genbankhits'

workflow FASTA_CHECKV_TSV {
    take:
    virus_fasta_gz  // channel: [ [ meta ], fasta.gz ]
    checkv_db       // path: checkv_db
    ncbi_info       // path: ncbi_info

    main:
    ch_versions = Channel.empty()

    // if checkv DB exists, skip CHECKV_DOWNLOADDATABASE
    if (checkv_db){
        ch_checkv_db = checkv_db
    } else {
        //
        // MODULE: download standard CheckV database
        //
        CHECKV_DOWNLOADDATABASE(
        )
        ch_checkv_db    = CHECKV_DOWNLOADDATABASE.out.checkv_db
        ch_versions     = ch_versions.mix(CHECKV_DOWNLOADDATABASE.out.versions)
    }

    //
    // MODULE: Run CheckV end to end
    //
    CHECKV_ENDTOEND(
        virus_fasta_gz,
        ch_checkv_db.first()
    )
    ch_versions = ch_versions.mix(CHECKV_ENDTOEND.out.versions)

    //
    // MODULE: Identify CheckV Genbank hits
    //
    CHECKV_GENBANKHITS(
        CHECKV_ENDTOEND.out.aai,
        ncbi_info,
        ch_checkv_db.first()
    )
    ch_versions     = ch_versions.mix(CHECKV_GENBANKHITS.out.versions)

    emit:
    quality_summary_tsv = CHECKV_ENDTOEND.out.quality_summary   // channel: [ [ meta ], quality_summary.tsv ]
    contamination_tsv   = CHECKV_ENDTOEND.out.contamination     // channel: [ [ meta ], contamination.tsv ]
    completeness_tsv    = CHECKV_ENDTOEND.out.completeness      // channel: [ [ meta ], completeness.tsv ]
    genbank_hits_tsv    = CHECKV_GENBANKHITS.out.genbank_hits   // path: [ [ meta ], genbank_hits.tsv ]
    versions            = ch_versions                           // path: [ versions.yml ]
}
