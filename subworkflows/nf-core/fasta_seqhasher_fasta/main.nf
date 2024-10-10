// Import modules
include { CSVTK_CONCAT          } from '../../../modules/nf-core/csvtk/concat'
include { SEQHASHER_SEQHASHER   } from '../../../modules/nf-core/seqhasher/seqhasher'
include { SEQHASHER_SEQUNIQ     } from '../../../modules/nf-core/seqhasher/sequniq'
include { SEQKIT_GREP           } from '../../../modules/nf-core/seqkit/grep'

workflow FASTA_SEQHASHER_FASTA {

    take:
    fasta       // channel: [ val(meta), [ acc ] ]
    seqhasher   // channel: file(seqhasher binary)

    main:

    ch_versions = Channel.empty()

    //
    // MODULE: Calcuate seqhashes for each fasta file
    //
    SEQHASHER_SEQHASHER(
        fasta,
        seqhasher
    )
    ch_seqhasher_tsv    = SEQHASHER_SEQHASHER.out.tsv
    ch_versions         = ch_versions.mix(SEQHASHER_SEQHASHER.out.versions.first())

    //
    // MODULE: Combine seqhashes into a single file and remove duplicates
    //
    CSVTK_CONCAT(
        ch_seqhasher_tsv.map { meta, tsv -> [ [id:'combined_seqhasher'], tsv ] }.groupTuple(sort:'deep'),
        "tsv",
        "tsv"
    )
    ch_seqhasher_combined_tsv   = CSVTK_CONCAT.out.csv
    ch_versions                 = ch_versions.mix(CSVTK_CONCAT.out.versions.first())

    //
    // MODULE: Remove duplicates from the combined seqhasher file
    //
    SEQHASHER_SEQUNIQ(
        ch_seqhasher_combined_tsv
    )
    ch_seqhasher_unique_tsv = SEQHASHER_SEQUNIQ.out.unique
    ch_versions             = ch_versions.mix(SEQHASHER_SEQUNIQ.out.versions.first())

    //
    // MODULE: Extract non-duplicated sequences from the original fasta files
    //
    SEQKIT_GREP(
        fasta,
        ch_seqhasher_unique_tsv.map { meta, tsv -> tsv }.first()
    )
    ch_unique_seqs_fasta_gz = SEQKIT_GREP.out.filter
    ch_versions             = ch_versions.mix(SEQKIT_GREP.out.versions.first())

    emit:
    unique_seqs_fasta_gz    = SEQKIT_GREP.out.filter    // channel: [ [ meta.id ], fasta.gz ]
    versions                = ch_versions               // channel: [ versions.yml ]
}

