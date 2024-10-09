// Import modules
include { SEQHASHER_SEQUNIQ } from '../../../modules/nf-core/seqhasher/sequniq'
include { SEQKIT_CONCAT     } from '../../../modules/nf-core/seqkit/concat'
include { SEQKIT_GREP       } from '../../../modules/nf-core/seqkit/grep'
include { VCLUST_ALIGN      } from '../../../modules/nf-core/vclust/align'
include { VCLUST_CLUSTER    } from '../../../modules/nf-core/vclust/cluster'
include { VCLUST_PREFILTER  } from '../../../modules/nf-core/vclust/prefilter'

workflow FASTA_VCLUST_FASTATSV {
    take:
    fasta_gz    // [ [ meta ], fasta.gz ]
    vclust      // path/to/vclust/binary

    main:
    ch_versions = Channel.empty()

    // combine all input fasta files into a single file
    ch_seqkit_concat_input  = fasta_gz
        .map { meta, fasta -> [ [ id: "all_samples" ], fasta ] }
        .groupTuple(sort: 'deep')

    //
    // MODULE: Concatenate all input FastA files
    //
    SEQKIT_CONCAT(
        ch_seqkit_concat_input
    )
    ch_versions = ch_versions.mix(SEQKIT_CONCAT.out.versions)

    //
    // MODULE: Prefilter alignments
    //
    VCLUST_PREFILTER(
        SEQKIT_CONCAT.out.fastx,
        vclust
    )
    ch_versions = ch_versions.mix(VCLUST_PREFILTER.out.versions)

    //
    // MODULE: Align sequences
    //
    VCLUST_ALIGN(
        SEQKIT_CONCAT.out.fastx,
        VCLUST_PREFILTER.out.prefilter,
        vclust
    )
    ch_versions = ch_versions.mix(VCLUST_ALIGN.out.versions)

    //
    // MODULE: Cluster sequences
    //
    VCLUST_CLUSTER(
        VCLUST_ALIGN.out.ani,
        VCLUST_ALIGN.out.ids,
        vclust
    )
    ch_versions = ch_versions.mix(VCLUST_CLUSTER.out.versions)

    //
    // MODULE: Dereplicate cluster output to retain only cluster rep
    //
    SEQHASHER_SEQUNIQ(
        VCLUST_CLUSTER.out.clusters
    )
    ch_versions = ch_versions.mix(SEQHASHER_SEQUNIQ.out.versions)

    //
    // MODULE: Extract rep sequences from FastA file
    //
    SEQKIT_GREP(
        SEQKIT_CONCAT.out.fastx,
        SEQHASHER_SEQUNIQ.out.unique.map { meta, tsv -> tsv }
    )
    ch_versions = ch_versions.mix(SEQKIT_GREP.out.versions)

    emit:
    cluster_reps_fasta_gz   = SEQKIT_GREP.out.filter    // [ [ meta ], fasta.gz ]
    versions                = ch_versions               // channel: [ versions.yml ]
}
