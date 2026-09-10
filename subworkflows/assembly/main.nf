/*
Copyright (c) 2026 The Hiller Lab at the Senckenberg Gessellschaft für Naturforschung
Distributed under the terms of the Apache License, Version 2.0.
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { BAMSPLIT_CHROM } from '../../modules/custom/bamsplit/chrom/main'
include { ALETSCH } from '../../modules/custom/aletsch/run/main'
include { STRINGTIE3 } from '../../modules/custom/stringtie3/main'
include { BEAVER } from '../../modules/custom/beaver/run/main'
include { TRANSMETA } from '../../modules/custom/transmeta/main'
include { RENAME_GTF } from '../../modules/custom/rename/gtf/main'
include { RENAME_GTF as RENAME_METASSEMBLY_GTF } from '../../modules/custom/rename/gtf/main'
include { RENAME_BAM } from '../../modules/custom/rename/bam/main'
include { REMOVE_BAMS } from '../../modules/custom/remove/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// INFO: derive the chromosome name from the per-chr file [bamsplit names files by chromosome]
def getChrName(Path f) {
    return f.name
        .replaceFirst(/\.bam\.(bai|csi)$/, '')
        .replaceFirst(/\.(bam|bai|csi)$/, '')
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    LOCAL ASSEMBLY
    - run the selected local assembler [aletsch, stringtie3] on the given bam channel
    - share either per-chr chunking or direct input through the joined beaver inputs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow LOCAL_ASSEMBLY {
    take:
        ch_inputs      // channel: [ val(meta), path(bam), path(bai) ]

    main:
        ch_versions = Channel.empty()

        switch (params.assembler) {
            case 'aletsch':
                ch_assembler = ALETSCH(
                    ch_inputs
                )
                ch_gtf = ch_assembler.gtf
                ch_counts = ch_assembler.assembled_transcripts
                break

            case 'stringtie3':
                ch_assembler = STRINGTIE3(
                    ch_inputs.map { meta, bam, bai -> [meta, bam] }
                )
                ch_gtf = ch_assembler.transcript_gtf
                ch_counts = ch_assembler.assembled_transcripts
                break

            default:
                error("Unknown assembler '${params.assembler}'; options are: 'aletsch', 'stringtie3'")
        }

        ch_renamed_gtf = RENAME_GTF(
            ch_gtf
        )

        ch_versions = ch_versions.mix(ch_assembler.versions)

    emit:
        gtf = ch_renamed_gtf.gtf
        counts = ch_counts
        versions = ch_versions // channel: [ versions.yml ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow ASSEMBLY {
    take:
        ch_bams // channel: [ val(meta), path(bam), path(bai) ]
        annotation_gtf // channel: [ val(meta), path(gtf) ]

    main:
        ch_versions = Channel.empty()
        ch_final_gtf = Channel.empty()
        ch_features = Channel.empty()
        ch_counts = Channel.empty()
        ch_grouped_gtfs = Channel.empty()
        ch_grouped_bams = Channel.empty()
        ch_name_map = Channel.empty()
        // INFO: barrier accumulator for the full-BAM cleanup below; carries one item per
        // LOCAL_ASSEMBLY task across both the per-chr and per-sample beaver branches
        ch_local_versions = Channel.empty()

        /*
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            ASSEMBLY BY CHROMOSOME
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */

        if (params.assembly_by_chr) {

            /*
            ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                SPLIT BAM PER-CHR
                - assemble locally [per sample+chr]
                - metassemble per-chr
                - merge into a single assembly
            ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            */

            // INFO: split bam per-chr -> assemble locally [per sample+chr] -> metassemble per-chr -> merge into a single assembly
            ch_split = BAMSPLIT_CHROM(
                ch_bams.map { meta, bam, bai -> [meta, bam] }
            )

            // INFO: explode the per-chr bam title into one item per chromosome, enriching meta
            // with the chromosome name and keeping the original sample id, so that bam/bai can
            // be paired 1:1 and chunks get unique assay prefixes [unmapped chunk is dropped]
            ch_chunk_bam = ch_split.bam
                .flatMap { meta, bams ->
                    bams.collectMany { bam ->
                        def chr = getChrName(bam)
                        def excluded = params.assembly_exclude_chromosomes ?: []
                        if (chr == 'unmapped' || chr in excluded) { [] } else {
                            [ [ meta + [ id: "${meta.id}_${chr}", sample: meta.id, chr: chr ], bam ] ]
                        }
                    }
                }
            ch_chunk_bai = ch_split.bai
                .flatMap { meta, bais ->
                    bais.collectMany { bai ->
                        def chr = getChrName(bai)
                        def excluded = params.assembly_exclude_chromosomes ?: []
                        if (chr == 'unmapped' || chr in excluded) { [] } else {
                            [ [ meta + [ id: "${meta.id}_${chr}", sample: meta.id, chr: chr ], bai ] ]
                        }
                    }
                }
                .mix(
                    ch_split.csi
                        .flatMap { meta, csis ->
                            csis.collectMany { csi ->
                                def chr = getChrName(csi)
                                def excluded = params.assembly_exclude_chromosomes ?: []
                                if (chr == 'unmapped' || chr in excluded) { [] } else {
                                    [ [ meta + [ id: "${meta.id}_${chr}", sample: meta.id, chr: chr ], csi ] ]
                                }
                            }
                        }
                )

            ch_chunk_input = ch_chunk_bam
                .join(ch_chunk_bai)
                .map { meta, bam, bai -> [meta, bam, bai] }

            /*
            ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                LOCAL ASSEMBLY + METASSEMBLY PER-CHR
                - beaver: assemble locally [per sample+chr] -> metassemble the per-chr gtf
                - transmeta: metassemble the aligned bams per-chr, using the annotation as guide
            ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            */

            if (params.metassembler == 'beaver') {
                // Aletsch 1.1.x does not flush the final reference when given a
                // BAM with only one populated chromosome. Keep chromosome
                // fan-out from bamsplit, but select that chromosome from the
                // original full BAM. Other local assemblers consume split BAMs.
                ch_local_inputs = ch_chunk_input
                if (params.assembler == 'aletsch') {
                    ch_chunk_metas = ch_chunk_bam.map { meta, _bam -> [meta.sample, meta] }
                    ch_full_bams = ch_bams.map { meta, bam, bai -> [meta.id, bam, bai] }
                    ch_local_inputs = ch_chunk_metas
                        // Multiple chromosome records share a sample key. `combine(by:)`
                        // preserves every chromosome; `join` retains only one duplicate key.
                        .combine(ch_full_bams, by: 0)
                        .map { sample, chunk_meta, bam, bai -> [chunk_meta, bam, bai] }
                }

                ch_local = LOCAL_ASSEMBLY(
                    ch_local_inputs
                )

                ch_renamed_gtf = ch_local.gtf

                // INFO: collect the per sample + chr gtf and group them by chr for the metassembler,
                // keeping the sample flags of the first sample in the group
                ch_grouped_gtfs = ch_renamed_gtf
                    .map { meta, gtf -> [meta.chr, meta.single_end, meta.strandedness, gtf] }
                    .groupTuple()
                    .map { chr, ses, strs, gtfs -> [ chr, ses.first(), strs.first(), gtfs ] }

                // INFO: local assemblers run per-chr -> aggregate transcript counts back to the sample level
                ch_counts = ch_local.counts
                    .map { meta, n -> [meta.sample, meta.single_end, meta.strandedness, n as int] }
                    .groupTuple()
                    .map { sample, single_ends, strandednesses, counts ->
                        [ [ id: sample, single_end: single_ends.first(), strandedness: strandednesses.first() ], counts.sum() as int ]
                    }

                ch_local_versions = ch_local_versions.mix(ch_local.versions)

                ch_versions = ch_versions.mix(ch_local.versions)
            } else {
                // INFO: rename the per-chr bams with the sample prefix so that they do not
                // collide when grouped and their outputs can be mapped back to samples
                ch_renamed_bams = RENAME_BAM(
                    ch_chunk_bam
                )

                ch_grouped_bams = ch_renamed_bams.bam
                    .map { meta, bam -> [meta.chr, meta.single_end, meta.strandedness, bam] }
                    .groupTuple()
                    .map { chr, ses, strs, bams ->
                        // transmeta takes one -s flag per chunk; a mixed chunk can only
                        // honour the first sample's strandedness, so surface it.
                        if (strs.unique().size() > 1) {
                            log.warn "[ASSEMBLY] ${chr}: mixed strandedness ${strs.unique()} in one transmeta chunk; using '${strs.first()}'"
                        }
                        [ chr, ses.first(), strs.first(), bams ]
                    }

                // INFO: map each bam name back to its sample id to aggregate transmeta counts
                ch_name_map = ch_renamed_bams.bam
                    .map { meta, bam -> [
                        bam.baseName,
                        [ id: meta.sample, single_end: meta.single_end, strandedness: meta.strandedness ]
                    ] }

                ch_versions = ch_versions.mix(RENAME_BAM.out.versions)
            }

            ch_versions = ch_versions.mix(BAMSPLIT_CHROM.out.versions)
        } else {

          /*
          ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              ASSEMBLY BY SAMPLE
          ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          */

            if (params.metassembler == 'beaver') {
                ch_local = LOCAL_ASSEMBLY(
                    ch_bams
                )

                ch_renamed_gtf = ch_local.gtf

                ch_grouped_gtfs = ch_renamed_gtf
                    .map { meta, gtf -> [ 'metassembly', meta, gtf ] }
                    .groupTuple()
                    .map { key, metas, gtfs -> [ null, metas.first().single_end, metas.first().strandedness, gtfs ] }

                ch_counts = ch_local.counts

                ch_local_versions = ch_local_versions.mix(ch_local.versions)

                ch_versions = ch_versions.mix(ch_local.versions)
            } else {
                ch_renamed_bams = RENAME_BAM(
                    ch_bams.map { meta, bam, bai -> [meta, bam] }
                )

                ch_grouped_bams = ch_renamed_bams.bam
                    .map { meta, bam -> [ 'metassembly', meta, bam ] }
                    .groupTuple()
                    .map { key, metas, bams ->
                        def strands = metas.collect { it.strandedness }.unique()
                        if (strands.size() > 1) {
                            log.warn "[ASSEMBLY] mixed strandedness ${strands} in one transmeta chunk; using '${strands.first()}'"
                        }
                        [ null, metas.first().single_end, strands.first(), bams ]
                    }

                ch_name_map = ch_renamed_bams.bam
                    .map { meta, bam -> [
                        bam.baseName,
                        [ id: meta.sample ?: meta.id, single_end: meta.single_end, strandedness: meta.strandedness ]
                    ] }

                ch_versions = ch_versions.mix(RENAME_BAM.out.versions)
            }
        }

        /*
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            FULL-BAM CLEANUP
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */

        if (!params.aletsch_keep_bam && !params.star_make_coverage) {
            // The full BAM is read by every LOCAL_ASSEMBLY task and by BAMSPLIT_CHROM, so
            // wait for the last local task before deleting it; a sibling rm can otherwise
            // race still-running consumers. Covers both the per-chr and per-sample
            // branches; with transmeta ch_local_versions stays empty, collect() emits
            // nothing and REMOVE_BAMS never runs [fail-safe: the bam is kept].
            ch_cleanup_done = ch_local_versions.collect().map { _versions -> 'done' }
            REMOVE_BAMS(
                ch_bams.map { meta, bam, bai -> [meta, [bam, bai]] }
                    .combine(ch_cleanup_done)
                    .map { meta, bams, _done -> [meta, bams] }
            )
            ch_versions = ch_versions.mix(REMOVE_BAMS.out.versions)
        }

        /*
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            METASSEMBLER
            - beaver [default] | transmeta
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */

        switch (params.metassembler) {
            case 'beaver':
                ch_metassembly = BEAVER(
                    ch_grouped_gtfs.map { chr, se, st, gtfs ->
                        chr ? [ [ id: "beaver_${chr}", chr: chr ], gtfs ] : [ [ id: 'metassembly' ], gtfs ]
                    }
                )
                ch_features = ch_metassembly.csv
                break

            case 'transmeta':
                ch_annotation_file = annotation_gtf.map { meta, gtf -> gtf }.first()

                // INFO: pair the annotation with each metassembly chunk, keeping the same
                // emission order as the bams channel so that they are zipped 1:1
                ch_annotation_inputs = ch_grouped_bams
                    .combine(ch_annotation_file)
                    .map { chr, se, st, bams, gtf -> [ [ id: chr ? "transmeta_${chr}" : 'transmeta' ], gtf ] }

                ch_transmeta_inputs = ch_grouped_bams
                    .map { chr, se, st, bams ->
                        [ [ id: chr ? "transmeta_${chr}" : 'transmeta', chr: chr, single_end: se, strandedness: st ], bams ]
                    }

                ch_metassembly = TRANSMETA(
                    ch_transmeta_inputs,
                    ch_annotation_inputs
                )

                ch_features = Channel.empty()

                // INFO: transmeta reports the assembled transcripts per sample -> aggregate back to the sample level
                ch_counts = ch_metassembly.counts
                    .map { meta, tsv -> tsv }
                    .splitCsv(sep: '\t')
                    .map { row -> [ row[0].toString(), row[1].toString() as int ] }
                    .join(ch_name_map)
                    .map { name, n, sample_meta -> [ sample_meta, n ] }
                    .groupTuple()
                    .map { sample_meta, ns -> [ sample_meta, ns.sum() as int ] }
                break

            default:
                error("Unknown metassembler '${params.metassembler}'; options are: 'beaver', 'transmeta'")
        }

        /*
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            CONCATENATE INTO A SINGLE METASSEMBLY GTF
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */

        if (params.assembly_by_chr) {
            // INFO: chromosome-scoped assemblers commonly restart transcript IDs at 0.
            // Prefix each result before concatenation so IDs remain globally unique.
            ch_renamed_metassembly = RENAME_METASSEMBLY_GTF(ch_metassembly.gtf)
            ch_final_gtf = ch_renamed_metassembly.gtf
                .map { meta, gtf -> gtf }
                .collectFile(
                    name: "${params.prefix ?: 'metassembly'}.gtf",
                    storeDir: "${params.output_dir}/06_${params.metassembler}",
                    newLine: true,
                )
                .map { gtf -> [ [ id: gtf.baseName ], gtf ] }
            ch_versions = ch_versions.mix(ch_renamed_metassembly.versions)
        } else {
            ch_final_gtf = ch_metassembly.gtf
        }

        ch_versions = ch_versions.mix(ch_metassembly.versions)

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        OUTPUT
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    emit:
        gtf = ch_final_gtf
        features = ch_features
        counts = ch_counts
        versions = ch_versions // channel: [ versions.yml ]
}
