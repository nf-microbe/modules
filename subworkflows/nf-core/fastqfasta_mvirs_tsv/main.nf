// Import modules
include { MVIRS_INDEX   } from '../../../modules/nf-core/mvirs/index'
include { MVIRS_OPRS    } from '../../../modules/nf-core/mvirs/oprs'
include { MVIRS_PARSER  } from '../../../modules/nf-core/mvirs/parser'

// Import subworkflows
include { rmEmptyFastAs; rmEmptyFastQs  } from '../utils_nfmicrobe_functions'

workflow FASTQFASTA_MVIRS_TSV {

    take:
    fastq_gz            // channel: [ [ meta.id, meat.group, meta.single_end ], [ path(fastq_1), path(fastq_2) ] ]
    fasta_gz            // channel: [ [ meta.id ], path(fasta)]
    integrase_hmm       // path: "assets/hmms/integrases/integrases.hmm"

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

    // remove channels with empty proviruses
    ch_mvirs_provirus_fasta_gz  = rmEmptyFastAs(MVIRS_OPRS.out.fasta, false)

    // join fasta and mvirs prophages
    ch_mvirs_drs_input = ch_mvirs_provirus_fasta_gz.map { meta, fastq -> [ [ id:meta.group ], meta, fastq ] }
        .combine(fasta_gz, by:0)
        .multiMap { meta_group, meta, provirus, contigs ->
            contigs:  [ meta, contigs ]
            provirus: [ meta, provirus ]
        }

    //
    // MODULE: Detect direct repeats in proviruses
    //
    MVIRS_PARSER(
        ch_mvirs_drs_input.contigs,
        ch_mvirs_drs_input.provirus,
        integrase_hmm
    )
    ch_versions = ch_versions.mix(MVIRS_PARSER.out.versions)

    emit:
    mvirs_prophages = MVIRS_OPRS.out.prophage       // channel: [ [ meta ], prophage.tsv ]
    mvirs_fasta     = MVIRS_OPRS.out.fasta          // channel: [ [ meta ], provirus.fasta.gz ]
    mvirs_clipped   = MVIRS_OPRS.out.clipped        // channel: [ [ meta ], clipped.tsv ]
    mvirs_oprs      = MVIRS_OPRS.out.oprs           // channel: [ [ meta ], oprs.tsv ]
    mvirs_summary   = MVIRS_PARSER.out.summary      // channel: [ [ meta ], summary.tsv ]
    mvirs_integrases= MVIRS_PARSER.out.integrases   // channel: [ [ meta ], integrases.tsv ]
    versions        = ch_versions                   // channel: [ versions.yml ]
}
