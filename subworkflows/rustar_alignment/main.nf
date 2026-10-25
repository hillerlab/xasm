/*
Copyright (c) 2026 The Hiller Lab at the Senckenberg Gessellschaft für Naturforschung
Distributed under the terms of the Apache License, Version 2.0.
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { RUSTAR_ALIGN as RUSTAR_ALIGN_1PASS } from '../../modules/custom/rustar/align/main'
include { RUSTAR_ALIGN as RUSTAR_ALIGN_2PASS } from '../../modules/custom/rustar/align/main'
include { SAMTOOLS_INDEX } from '../../modules/custom/samtools/index/main'
include { JOIN_JUNCTIONS } from '../../modules/custom/junctions/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow RUSTAR_ALIGNMENT {
    take:
        reads
        gtf    // channel: [ val(meta), path(gtf) ]

    main:
        ch_versions = Channel.empty()

        ch_trimmed_reads = reads
        ch_trimmed_reads
            .multiMap { meta, reads, index ->
                first_pass: [meta, reads, index]
                second_pass: [meta, reads, index]
            }
            .set { ch_trimmed_split }

        gtf
            .multiMap { meta, file ->
                first_pass: file
                second_pass: file
            }
            .set { gtf_split }

        if (params.star_twopass_junctions_file) {
            // INFO: pre-curated junctions provided -> skip first pass and junction merging, align once
            ch_junctions_file = Channel.value(file(params.star_twopass_junctions_file, checkIfExists: true))
        } else {
            if (params.star_ignore_gtf_for_mapping) {
                ch_rustar_first_pass_out = RUSTAR_ALIGN_1PASS(
                    ch_trimmed_split.first_pass,
                    Channel.value([]),
                    Channel.value([]),
                    params.star_ignore_gtf_for_mapping,
                    params.star_seq_platform ?: '',
                    params.star_seq_center ?: '',
                    params.star_seq_library ?: '',
                    params.star_machine_type ?: '',
                    params.star_keep_first_pass_bam ?: false,
                    false // delete fastq
                )
            } else {
                ch_rustar_first_pass_out = RUSTAR_ALIGN_1PASS(
                    ch_trimmed_split.first_pass,
                    gtf_split.first_pass,
                    Channel.value([]),
                    params.star_ignore_gtf_for_mapping,
                    params.star_seq_platform ?: '',
                    params.star_seq_center ?: '',
                    params.star_seq_library ?: '',
                    params.star_machine_type ?: '',
                    params.star_keep_first_pass_bam ?: false,
                    false, // delete fastq
                )
            }

            ch_splice_junctions = ch_rustar_first_pass_out.spl_junc_tab

            // Collect all junction files into a single item
            ch_all_junctions = ch_splice_junctions
                .map { meta, junctions -> junctions } // Extract just the files
                .collect() // Collect all files into a single list
                .map { files -> [ [id: 'merged_junctions'], files ] } // Add a generic meta map

            ch_filtered_junctions = JOIN_JUNCTIONS(
                ch_all_junctions,
                params.junction_min_junction_length,
                params.junction_min_read_coverage
            )

            // Now ch_filtered_junctions contains a single merged junction file
            ch_junctions_file = ch_filtered_junctions.filtered_junctions
                .map { meta, junctions -> junctions } // Extract just the file
        }

        ch_second_pass_input = ch_trimmed_split.second_pass
            .combine(ch_junctions_file)

        if (params.star_ignore_gtf_for_mapping) {
            ch_rustar_second_pass_out = RUSTAR_ALIGN_2PASS(
                ch_second_pass_input.map { meta, reads, index, _junctions -> [meta, reads, index] },
                Channel.value([]),
                ch_second_pass_input.map { _meta, _reads, _index, junctions -> junctions },
                params.star_ignore_gtf_for_mapping,
                params.star_seq_platform ?: '',
                params.star_seq_center ?: '',
                params.star_seq_library ?: '',
                params.star_machine_type ?: '',
                true, // keep bam -> we need it for Aletsch
                params.star_delete_fastq_after_alignment
            )
        } else {
            ch_rustar_second_pass_out = RUSTAR_ALIGN_2PASS(
                ch_second_pass_input.map { meta, reads, index, _junctions -> [meta, reads, index] },
                gtf_split.second_pass,
                ch_second_pass_input.map { _meta, _reads, _index, junctions -> junctions },
                params.star_ignore_gtf_for_mapping,
                params.star_seq_platform ?: '',
                params.star_seq_center ?: '',
                params.star_seq_library ?: '',
                params.star_machine_type ?: '',
                true, // keep bam -> we need it for Aletsch
                params.star_delete_fastq_after_alignment
            )
        }

        ch_versions = ch_versions.mix(RUSTAR_ALIGN_2PASS.out.versions.first())
        ch_log_final = ch_rustar_second_pass_out.log_final

        // ruSTAR does not index its own BAM (no samtools in the image), so the
        // (meta, bam, bai) tuple STAR_ALIGNMENT emits is rebuilt here.
        SAMTOOLS_INDEX(ch_rustar_second_pass_out.bam_sorted_aligned)
        ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

        ch_rustar_second_pass_out.bam_sorted_aligned
            .join(SAMTOOLS_INDEX.out.bai)
            .set { ch_bam_sorted }

    emit:
        bams = ch_bam_sorted
        bedgraph = ch_rustar_second_pass_out.bedgraph
        junctions = ch_junctions_file
        percent_mapped = ch_log_final.map { meta, log -> [ meta, getStarPercentMapped(params, log) ] }
        log_final = ch_log_final
        bam_size = ch_rustar_second_pass_out.bam_size
        versions = ch_versions // channel: [ versions.yml ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Function that parses and returns the alignment rate from the ruSTAR log output
// (Log.final.out is STAR-compatible, so this mirrors star_alignment/main.nf)
//
def getStarPercentMapped(params, align_log) {
    def percent_aligned = 0
    def pattern = /Uniquely mapped reads %\s*\|\s*([\d\.]+)%/
    align_log.eachLine { line ->
        def matcher = line =~ pattern
        if (matcher) {
            percent_aligned = matcher[0][1].toFloat()
        }
    }

    return percent_aligned
}
