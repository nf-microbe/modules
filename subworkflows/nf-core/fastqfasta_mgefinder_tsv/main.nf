// Import modules
include { BWA_INDEX                     } from '../../../modules/nf-core/bwa/index'
include { BWA_MEM                       } from '../../../modules/nf-core/bwa/mem'
include { MGEFINDER_FORMATBAM           } from '../../../modules/nf-core/mgefinder/formatbam'
include { MGEFINDER_FIND                } from '../../../modules/nf-core/mgefinder/find'
include { MGEFINDER_PAIR                } from '../../../modules/nf-core/mgefinder/pair'
include { MGEFINDER_INFERSEQREFERENCE   } from '../../../modules/nf-core/mgefinder/inferseqreference'

workflow FASTQFASTA_MGEFINDER_TSV {

    take:
    fastq_gz            // channel: [ [ meta.id, meta.single_end ], [ path(fastq_1), path(fastq_2) ] ]
    fasta_gz            // channel: [ [ meta.id ], path(fasta)]

    main:

    ch_versions = Channel.empty()

    //
    // MODULE: Index reference genome/assembly
    //
    BWA_INDEX(
        fasta_gz
    )
    ch_versions = ch_versions.mix(BWA_INDEX.out.versions)

    // join fastQ and Fasta
    ch_bwa_mem_input = fastq_gz.map { meta, fastq -> [ [ id:meta.group ], meta, fastq ] }
        .combine(fasta_gz, by:0)
        .combine(BWA_INDEX.out.index, by:0)
        .multiMap { meta_group, meta_fastq, fastq, fasta, index ->
            fastq:  [ meta_fastq, fastq ]
            fasta:  [ meta_fastq, fasta ]
            index:  [ meta_fastq, index ]
        }

    //
    // MODULE: Align reads to reference genome/assembly
    //
    BWA_MEM(
        ch_bwa_mem_input.fastq,
        ch_bwa_mem_input.index,
        ch_bwa_mem_input.fasta,
        false
    )
    ch_versions = ch_versions.mix(BWA_MEM.out.versions)

    //
    // MODULE: Format BAM
    //
    MGEFINDER_FORMATBAM(
        BWA_MEM.out.bam
    )
    ch_versions = ch_versions.mix(MGEFINDER_FORMATBAM.out.versions)

    //
    // MODULE: Find insertion sites
    //
    MGEFINDER_FIND(
        MGEFINDER_FORMATBAM.out.bam
    )
    ch_versions = ch_versions.mix(MGEFINDER_FIND.out.versions)

    ch_mgefinder_find = MGEFINDER_FORMATBAM.out.bam
        .filter { meta, find ->
            find.countLines( limit: 10 ) > 1
        }

    // join fasta, find bam, and find file
    ch_pair_input = fasta_gz
        .combine(MGEFINDER_FORMATBAM.out.bam.map { meta, bam -> [ [ id: meta.group ], meta, bam ] }, by:0)
        .combine(ch_mgefinder_find.map { meta, find -> [ [ id: meta.group ], find ] }, by:0)
        .multiMap { meta_group, fasta, meta_bam, bam, find ->
            fasta:  [ meta_bam, fasta ]
            bam:    [ meta_bam, bam ]
            find:   [ meta_bam, find ]
        }

    //
    // MODULE: Pair termini with each other
    //
    MGEFINDER_PAIR(
        ch_pair_input.fasta,
        ch_pair_input.bam,
        ch_pair_input.find
    )
    ch_versions = ch_versions.mix(MGEFINDER_PAIR.out.versions.first())

    ch_mgefinder_pair = MGEFINDER_PAIR.out.pair
        .filter { meta, pair ->
            pair.countLines( limit: 10 ) > 1
        }

    // join FastA and pair file
    ch_inferseqreference_input = fasta_gz
        .combine(ch_mgefinder_pair.map { meta, pair -> [ [ id: meta.group ], meta, pair ] }, by:0)
        .multiMap { meta_group, fasta, meta_pair, pair ->
            fasta:  [ meta_pair, fasta ]
            pair:   [ meta_pair, pair ]
        }

    //
    // MODULE: infer the identity of insertions
    //
    MGEFINDER_INFERSEQREFERENCE(
        ch_inferseqreference_input.fasta,
        ch_inferseqreference_input.pair
    )

    emit:
    prophage    = MGEFINDER_INFERSEQREFERENCE.out.tsv   // channel: [ [ meta ], inferse_reference.tsv ]
    versions    = ch_versions                           // channel: [ versions.yml ]
}
