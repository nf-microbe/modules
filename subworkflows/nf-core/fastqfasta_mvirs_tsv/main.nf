// Import modules
include { MVIRS_INDEX   } from '../../../modules/nf-core/mvirs/index'
include { MVIRS_OPRS    } from '../../../modules/nf-core/mvirs/oprs'

// Import subworkflows
include { rmEmptyFastAs; rmEmptyFastQs  } from '../utils_nfmicrobe_functions'

workflow FASTQFASTA_MVIRS_TSV {

    take:
    fastq_gz            // channel: [ [ meta.id, meat.group, meta.single_end ], [ path(fastq_1), path(fastq_2) ] ]
    fasta_gz            // channel: [ [ meta.id ], path(fasta) ]

    main:

    ch_versions = Channel.empty()

    // identify FastA files with associated FastQ files
    ch_mvirs_inputs  = fastq_gz.map { meta, fastq -> [ meta.group, meta, fastq ] }
        .join(fasta_gz.map { meta, fasta -> [ meta.id, meta, fasta ] } )
        .map { meta_id, meta_fastq, fastq, meta_fasta, fasta ->
            [ meta_fasta, fasta ]
        }

    //
    // MODULE: Index reference genome/assembly
    //
    MVIRS_INDEX(
        ch_mvirs_inputs
    )
    ch_versions = ch_versions.mix(MVIRS_INDEX.out.versions)

    // join fastQ, Fasta, and index
    ch_mvirs_oprs_input = MVIRS_INDEX.out.index.map { meta, index -> [ meta.id, index ] }
        .combine(fastq_gz.map { meta, fastq -> [ meta.group, meta, fastq ] }, by:0)
        .combine(fasta_gz.map { meta, fasta -> [ meta.id, fasta ] }, by:0)
        .multiMap { meta_group, index, meta_fastq, fastq, fasta ->
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

    // remove channels with empty proviruses
    ch_mvirs_provirus_fasta_gz  = rmEmptyFastAs(MVIRS_OPRS.out.fasta, false)

    emit:
    mvirs_prophages = MVIRS_OPRS.out.prophage       // channel: [ [ meta ], prophage.tsv ]
    mvirs_fasta     = ch_mvirs_provirus_fasta_gz    // channel: [ [ meta ], provirus.fasta.gz ]
    mvirs_clipped   = MVIRS_OPRS.out.clipped        // channel: [ [ meta ], clipped.tsv ]
    mvirs_oprs      = MVIRS_OPRS.out.oprs           // channel: [ [ meta ], oprs.tsv ]
    versions        = ch_versions                   // channel: [ versions.yml ]
}
