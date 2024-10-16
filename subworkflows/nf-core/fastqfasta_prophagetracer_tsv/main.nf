// Import modules
include { BWA_INDEX             } from '../../../modules/nf-core/bwa/index'
include { BWA_MEM               } from '../../../modules/nf-core/bwa/mem'
include { SAMBAMBA_MARKDUP      } from '../../../modules/nf-core/sambamba/markdup'
include { PROPHAGETRACER_WGS    } from '../../../modules/nf-core/prophagetracer/wgs'

workflow FASTQFASTA_PROPHAGETRACER_TSV {

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
            meta_fastq.id   = meta_fastq.id + '-' + meta_group.id
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
    // MODULE: Mark duplicates
    //
    SAMBAMBA_MARKDUP(
        BWA_MEM.out.bam
    )
    ch_versions = ch_versions.mix(SAMBAMBA_MARKDUP.out.versions)

    // join bam and Fasta
    ch_prophage_tracer_input = fasta_gz
        .combine(SAMBAMBA_MARKDUP.out.bam.map { meta, bam -> [ [ id: meta.group ], meta, bam ] }, by:0)
        .multiMap { meta_group, fasta, meta_bam, bam ->
            fasta:  [ meta_bam, fasta ]
            bam:  [ meta_bam, bam ]
        }

    //
    // MODULE: Identify prophages
    //
    PROPHAGETRACER_WGS(
        ch_prophage_tracer_input.fasta,
        ch_prophage_tracer_input.bam
    )
    ch_versions = ch_versions.mix(PROPHAGETRACER_WGS.out.versions.first())

    emit:
    prophage    = PROPHAGETRACER_WGS.out.prophage   // channel: [ [ meta ], prophage.tsv ]
    versions    = ch_versions                       // channel: [ versions.yml ]
}
