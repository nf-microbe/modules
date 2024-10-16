// Import modules
include { MVIRS_INDEX   } from '../../../modules/nf-core/mvirs/index'
include { MVIRS_OPRS    } from '../../../modules/nf-core/mvirs/oprs'

workflow FASTQFASTA_MVIRS_TSV {

    take:
    fastq_gz            // channel: [ [ meta.id, meat.group, meta.single_end ], [ path(fastq_1), path(fastq_2) ] ]
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
    ch_mvirs_oprs_input = fastq_gz.map { meta, fastq -> [ [ id:meta.group ], meta, fastq ] }
        .combine(fasta_gz, by:0)
        .combine(MVIRS_INDEX.out.index, by:0)
        .multiMap { meta_group, meta_fastq, fastq, fasta, index ->
            fastq:  [ meta_fastq, fastq ]
            fasta:  [ meta_fastq, fasta ]
            index:  [ meta_fastq, index ]
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
