// Import modules
include { HOSTILE_FETCH } from '../../../modules/nf-core/hostile/fetch/main'
include { HOSTILE_CLEAN } from '../../../modules/nf-core/hostile/clean/main'

// Run subworkflow
workflow FASTQ_HOSTILE_FASTQ {

    take:
    fastq_gz        // channel: [ [ meta ], [ reads_1.fastq.gz, reads_2.fastq.gz ] ] (MANDATORY)
    hostile_index   // val: /path/to/host_index.bt2 (MANDATORY)

    main:
    def ch_versions = Channel.empty()


    //
    // MODULE: Download hostile reference
    //
    HOSTILE_FETCH(
        Channel.of([ [ id:'hostile_index' ], hostile_index ])
    )
    def ch_hostile_index    = HOSTILE_FETCH.out.index.first()
    ch_versions             = ch_versions.mix(HOSTILE_FETCH.out.versions)

    //
    // MODULE: Remove human reads
    //
    HOSTILE_CLEAN(
        fastq_gz,
        ch_hostile_index
    )
    ch_versions = ch_versions.mix(HOSTILE_CLEAN.out.versions)

    emit:
    fastq_gz    = HOSTILE_CLEAN.out.fastq   // channel: [ val(meta), [ reads_1.fastq.gz, reads_2.fastq.gz ] ]
    versions    = ch_versions               // channel: [ software_versions.yml ]
}

