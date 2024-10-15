// Import modules
include { MVIRS_INDEX   } from '../../../modules/nf-core/mvirs/index'
include { MVIRS_OPRS    } from '../../../modules/nf-core/mvirs/oprs'

workflow FASTQFASTA_MVIRS_TSV {

    take:
    fastq_gz            // channel: [ [ meta.id, meta.single_end ], [ path(fastq_1), path(fastq_2) ] ]
    fasta_gz            // channel: [ [ meta.id ], path(fasta)]

    main:

    ch_versions = Channel.empty()

    //
    // MODULE: Index reference genome/assembly
    //
    MVIRS_INDEX(
        fasta_gz
    )
    ch_versions = ch_versions.mix(MVIRS_INDEX.out.versions)

    // join fastQ and Fasta
    ch_mvirs_oprs_input = fastq_gz
        .combine(fasta_gz, by:0)
        .combine(MVIRS_INDEX.out.index, by:0)
        .multiMap { meta, fastq, fasta, index ->
            fastq:  [ meta, fastq ]
            fasta:  [ meta, fasta ]
            index:  [ meta, index ]
        }

    //
    // MODULE: Align reads to reference genome/assembly
    //
    MVIRS_OPRS(
        ch_mvirs_oprs_input.fastq,
        ch_mvirs_oprs_input.fasta,
        ch_mvirs_oprs_input.index
    )
    ch_versions = ch_versions.mix(MVIRS_OPRS.out.versions)

    emit:
    prophage    = MVIRS_OPRS.out.prophage   // channel: [ [ meta ], prophage.tsv ]
    versions    = ch_versions               // channel: [ versions.yml ]
}
