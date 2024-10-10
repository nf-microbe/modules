// Import modules
include { COVERM_CONTIG as COVERM_CONTIG_COBRA  } from '../../../modules/nf-core/coverm/contig'
include { COBRAMETA                             } from '../../../modules/nf-core/cobrameta'
include { SEQKIT_FX2TAB as SEQKIT_FX2TAB_COBRA  } from '../../../modules/nf-core/seqkit/fx2tab'
include { SEQKIT_SEQ as SEQKIT_SEQ_COBRA        } from '../../../modules/nf-core/seqkit/seq'

workflow FASTQFASTA_COBRA_FASTA {

    take:
    fastq_gz            // channel: [ [ meta.id, meta.single_end ], [ path(fastq_1), path(fastq_2) ] ]
    fasta_gz            // channel: [ [ meta.id, meta.single_end, meta.assembler ], path(fasta)]
    assembly_logs       // channel: [ [ meta.id, meta.single_end, meta.assembler ], path(logs) ]

    main:

    ch_versions = Channel.empty()

    //
    // MODULE: Identify sequences longer than 10kb
    //
    SEQKIT_SEQ_COBRA(
        fasta_gz
    )
    ch_versions = ch_versions.mix(SEQKIT_SEQ_COBRA.out.versions)

    //
    // MODULES: Convert fasta to tabular format
    //
    SEQKIT_FX2TAB_COBRA(
        SEQKIT_SEQ_COBRA.out.fastx
    )
    ch_versions             = ch_versions.mix(SEQKIT_FX2TAB_COBRA.out.versions)

    // join fastq and fasta by meta.id
    ch_coverm_input = fastq_gz.map { meta, fastq -> [ meta.id, fastq ] }
        .combine(fasta_gz.map { meta, fasta -> [ meta.id, meta, fasta ] }, by:0)
        .multiMap { meta_id, fastq, meta_fasta, fasta ->
            fastq: [ meta_fasta, fastq ]
            fasta: [ meta_fasta, fasta ]
        }

    //
    // MODULE: Align reads to their corresponding assembly
    //
    COVERM_CONTIG_COBRA(
        ch_coverm_input.fastq,
        ch_coverm_input.fasta
    )
    ch_versions = ch_versions.mix(COVERM_CONTIG_COBRA.out.versions.first())

    // prepare input for cobra
    ch_cobra_input = fasta_gz
        .combine(SEQKIT_FX2TAB_COBRA.out.text, by:0)
        .combine(COVERM_CONTIG_COBRA.out.tsv, by:0)
        .combine(COVERM_CONTIG_COBRA.out.bam, by:0)
        .combine(assembly_logs, by:0)
        .multiMap { meta, fasta, query, coverage, bam, log ->
            fasta:      [ meta, fasta ]
            coverage:   [ meta, coverage ]
            query:      [ meta, query ]
            bam:        [ meta, bam ]
            log:        [ meta, log ]
        }

    //
    // MODULE: Extend query contigs
    //
    COBRAMETA(
        ch_cobra_input.fasta,
        ch_cobra_input.coverage,
        ch_cobra_input.query,
        ch_cobra_input.bam,
        ch_cobra_input.log,
        ch_cobra_input.fasta.map { meta, fasta -> meta.assembler },
    )

    emit:
    fasta_gz    = COBRAMETA.out.fasta   // channel: [ [ meta.id, meta.assembler, meta.extension ], [ fasta.gz ] ]
    versions    = ch_versions           // channel: [ versions.yml ]
}
