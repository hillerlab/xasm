/*
Copyright (c) 2026 The Hiller Lab at the Senckenberg Gessellschaft für Naturforschung
Distributed under the terms of the Apache License, Version 2.0.
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { GTF_REMOVE_DIRT } from '../../modules/custom/gtf/clean/main'
include { GXF2BED } from '../../modules/custom/gxf2bed/main'
include { ISOTOOLS_ORPHAN } from '../../modules/custom/isotools/orphan/main'
include { ISOTOOLS_ORPHAN as ISOTOOLS_ORPHAN_DENOVO } from '../../modules/custom/isotools/orphan/main'
include { ISOTOOLS_FUSION } from '../../modules/custom/isotools/fusion/main'

include { XLOCI_INTRON as XLOCI_EXTRACT_INTRONS } from '../../modules/custom/xloci/intron/main.nf'
include { INTRONIC as IIC_PREDICT_SPLICEOSOME } from '../../modules/custom/intronic/main.nf'
include { ISOTOOLS_CLASSIFY_INTRON } from '../../modules/custom/isotools/classify/intron/main.nf'
include { ISOTOOLS_CLASSIFY_INTRON as ISOTOOLS_CLASSIFY_INTRON_TWOPASS } from '../../modules/custom/isotools/classify/intron/main.nf'
include { ISOTOOLS_INTRON_RETENTION } from '../../modules/custom/isotools/intron/main.nf'
include { ISOTOOLS_INTRON_RETENTION as ISOTOOLS_INTRON_RETENTION_TWOPASS } from '../../modules/custom/isotools/intron/main.nf'
include { ISOTOOLS_NMD } from '../../modules/custom/isotools/nmd/main.nf'

include { STRIP_OCCURRENCES as STRIP_RETENTIONS } from '../../modules/custom/strip/main'
include { STRIP_OCCURRENCES as STRIP_RETENTIONS_TWOPASS } from '../../modules/custom/strip/main'
include { STRIP_OCCURRENCES as STRIP_STRONG_RTS } from '../../modules/custom/strip/main'
include { STRIP_OCCURRENCES as STRIP_WEAK_RTS } from '../../modules/custom/strip/main'
include { STRIP_OCCURRENCES as STRIP_ARTIFACTS } from '../../modules/custom/strip/main'
include { PUBLISH as PUBLISH_FINAL_TRANSCRIPTS } from '../../modules/custom/publish/main'
include { RENAME_FINAL_TRANSCRIPTS } from '../../modules/custom/rename/final/main'

include { SORT_BED as SORT_BED_FL_TRANSCRIPTS } from '../../modules/custom/sort/main'
include { SORT_BED as SORT_BED_SCRAPS } from '../../modules/custom/sort/main'
include { SORT_BED as SORT_BED_FUSIONS } from '../../modules/custom/sort/main'
include { SORT_BED as SORT_BED_TWOPASS } from '../../modules/custom/sort/main'
include { SORT_BED as SORT_BED_FINAL } from '../../modules/custom/sort/main'
include { SORT_BED as SORT_BED_NMD } from '../../modules/custom/sort/main'
include { SORT_BED as SORT_BED_TRUNCATIONS } from '../../modules/custom/sort/main'
include { SORT_BED as SORT_BED_TRUNCATIONS_XORF } from '../../modules/custom/sort/main'
include { SORT_BED as SORT_BED_XORF } from '../../modules/custom/sort/main'
include { SORT_BED as SORT_BED_TOGA } from '../../modules/custom/sort/main'
include { SORT_BED as SORT_BED_REFERENCE } from '../../modules/custom/sort/main'
include { SORT_BED as SORT_BED_STRONG_RTS } from '../../modules/custom/sort/main'
include { SORT_BED as SORT_BED_WEAK_RTS } from '../../modules/custom/sort/main'
include { SORT_BED as SORT_BED_ARTIFACTS } from '../../modules/custom/sort/main'
include { SORT_BED as SORT_BED_RETENTIONS } from '../../modules/custom/sort/main'
include { XORF_RUN } from '../xorf/main'
include { XORF_RUN as XORF_RUN_ON_TWOPASS_RETENTIONS } from '../xorf/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    LOCAL SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow POLISH_TWOPASS {
    take:
        hq                       // channel: [ meta, bed ] first-pass XORF HQ (post-strip)
        truncations              // channel: [ meta, bed ] first-pass XORF truncation discards
        retentions               // channel: [ meta, bed ] first-pass retention discards
        intron_iic               // channel: [ meta, iic ] first-pass intronIC evidence
        genome                   // channel: [ meta, fasta ]
        annotation               // channel: [ meta, bed ]
        repeats                  // channel: [ meta, bed ]
        splice_scores            // channel: [ meta, bigwig directory ]

    main:
        ch_xorf_twopass_versions = Channel.empty()

        // XORF output ids embed the source transcript: meta.id =
        // "<prefix>_flnc@<chr>". POLISH names the ORF input `<prefix>_flnc`
        // and the subworkflow's groupKey appends `@<chr>`; the first-pass IIC
        // evidence is keyed by the plain sample prefix, so drop the suffix to
        // join.
        // ponytail: assumes the XORF emit id is exactly <prefix>_flnc@<chr>;
        // if the POLISH naming or the groupKey scheme changes, this regex is
        // the single point to fix.
        ch_classify_inputs = hq
            .map { meta, bed ->
                def source_id = meta.id.replaceFirst(/@.*$/, '').replaceFirst(/_flnc$/, '')
                tuple(source_id, meta, bed)
            }
            .join(
                intron_iic.map { meta, iic -> tuple(meta.id, meta, iic) }
            )
            .map { id, hq_meta, bed, iic_meta, iic ->
                tuple(hq_meta, bed, iic)
            }

        // NOTE: --reference is the reference intron source. First pass uses the
        // annotation alone. Second pass unions the annotation with the ORF
        // models so an intron is supported if either source has it (novel
        // ORF-supported splice sites stay categorized; annotation-supported
        // structure is not dropped).
        //
        // The two beds are collapsed to ONE file: passing hq per-transcript
        // staged reads and annotation under the same basename triggers
        // Nextflow's "input file name collision". The merged element is then
        // replicated to one per reads element, because plain multi-channel
        // process inputs zip pairwise and a single annotation element would
        // truncate classify to one task.
        ch_toga_unsorted = hq
            .map { _meta, bed -> bed }
            .concat(annotation.map { _meta, bed -> bed })
            .collectFile(
                name: "${params.prefix ?: 'polish'}.twopass_toga.bed",
                newLine: false,
                sort: true
            )
            .map { bed -> [ [ id: bed.baseName ], bed ] }

        SORT_BED_TOGA(ch_toga_unsorted)
        ch_toga = SORT_BED_TOGA.out.sorted

        ch_toga_n = hq
            .map { meta, _bed -> meta }
            .combine(ch_toga)
            .map { _meta, ann_meta, ann_bed -> [ ann_meta, ann_bed ] }

        ISOTOOLS_CLASSIFY_INTRON_TWOPASS(
            ch_classify_inputs,
            genome,
            ch_toga_n,
            repeats,
            splice_scores
        )

        // Candidates stay the first-pass discards, joined with first-pass IIC
        // exactly as before the ORF rewrite (STRIP_OCCURRENCES appends `.discard`
        // to the metadata id; join on the original id).
        ch_retention_candidates = retentions
            .map { meta, bed ->
                def source_id = meta.id.replaceFirst(/\.discard$/, '')
                tuple(source_id, meta + [ id: "${source_id}_twopass" ], bed)
            }
            .join(
                intron_iic.map { meta, iic -> tuple(meta.id, meta, iic) }
            )
            .map { id, retention_meta, bed, iic_meta, iic ->
                [retention_meta, bed]
            }

        ISOTOOLS_INTRON_RETENTION_TWOPASS(
            ch_retention_candidates,
            ISOTOOLS_CLASSIFY_INTRON_TWOPASS.out.tsv
        )

        STRIP_RETENTIONS_TWOPASS(
            ch_retention_candidates,
            ISOTOOLS_INTRON_RETENTION_TWOPASS.out.descriptor,
            'RETENTION'
        )

        // After deciding which retentions to keep we need to predict
        // ORFs for a second-round hard-filter.
        ch_xorf_hq_twopass_retentions = Channel.empty()
        ch_xorf_trunc_twopass = Channel.empty()
        if (params.xorf_call_orfs) {
            XORF_RUN_ON_TWOPASS_RETENTIONS(
                STRIP_RETENTIONS_TWOPASS.out.hq,
                genome,
                'xorf_twopass_retentions'
            )
            ch_xorf_twopass_versions = XORF_RUN_ON_TWOPASS_RETENTIONS.out.versions
            // With selenocysteine masking the merged id is `<t>_flnc@<masked|UNMSK>`
            // and which group wins the submodule's collect is order-dependent, so
            // canonicalize to the deterministic `<t>_flnc` before it keys filenames.
            // Upstream RENAME now emits only *.renamed.bed, so bed is a single file.
            // XORF.out.files is post-strip HQ (do_polishing default), so this
            // merge waits on STRIP_TRUNCATIONS rather than CONCAT_RENAMED.
            ch_xorf_hq_twopass_retentions = XORF_RUN_ON_TWOPASS_RETENTIONS.out.files.map { meta, bed, tsv ->
                if (bed instanceof List) {
                    error "XORF emitted List bed for ${meta.id}: ${bed} — expected single file " +
                          "(check params.xorf_skip_joined_concat / rename module)"
                }
                [ meta + [ id: meta.id.tokenize('@')[0] ], bed ]
            }
            ch_xorf_trunc_twopass = XORF_RUN_ON_TWOPASS_RETENTIONS.out.truncations
        }

        // ORF HQ and its retention discards are disjoint. Mixing and
        // sorting also falls back naturally to HQ when there are no retentions.
        ch_merged_hq = hq
            .concat(ch_xorf_hq_twopass_retentions)
            .map { meta, bed -> bed }
            .collectFile(
                name: "${params.prefix ?: 'polish'}.twopass.bed",
                newLine: false
            )
            .map { bed -> [ [ id: "${bed.baseName}_clean" ], bed ] }

        SORT_BED_TWOPASS(ch_merged_hq)

        ch_merged_truncations = truncations
            .concat(ch_xorf_trunc_twopass)
            .map { _meta, bed -> bed }
            .collectFile(
                name: "${params.prefix ?: 'polish'}.truncations.bed",
                newLine: false
            )
            .map { bed -> [ [ id: bed.baseName ], bed ] }

        SORT_BED_TRUNCATIONS(ch_merged_truncations)

        // NMD must not start until both the merged HQ and the combined
        // truncation discards are ready. ifEmpty covers the no-truncation case
        // (collectFile emits nothing; SORT_BED_TRUNCATIONS never runs).
        ch_hq_ready = SORT_BED_TWOPASS.out.sorted
            .combine(SORT_BED_TRUNCATIONS.out.sorted.ifEmpty([[:], []]))
            .map { meta, bed, _tmeta, _tbed -> [meta, bed] }

        ch_versions = ISOTOOLS_CLASSIFY_INTRON_TWOPASS.out.versions
            .mix(ISOTOOLS_INTRON_RETENTION_TWOPASS.out.versions)
            .mix(STRIP_RETENTIONS_TWOPASS.out.versions)
            .mix(SORT_BED_TOGA.out.versions)
            .mix(SORT_BED_TWOPASS.out.versions)
            .mix(SORT_BED_TRUNCATIONS.out.versions)
            .mix(ch_xorf_twopass_versions)

    emit:
        hq           = ch_hq_ready
        truncations  = SORT_BED_TRUNCATIONS.out.sorted
        retentions   = STRIP_RETENTIONS_TWOPASS.out.discard
        introns      = ISOTOOLS_CLASSIFY_INTRON_TWOPASS.out.tsv
        versions     = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow POLISH {
    take:
        metassembly               // channel: [ meta, gtf/bed ] 
        genome                    // channel: [ fasta ]
        annotation                // channel: [ meta, bed ]
        annevo                    // channel: [ bed ] ANNEVO annotation, 0 or 1 element
        tiberius                  // channel: [ bed ] TIBERIUS annotation, 0, 1 or 2 elements
        repeats                   // path: /path/to/repeats.{bed/gff/gtf}
        splice_scores_dir         // path: /path/to/splice/scores/dir
        do_twopass_polish         // val: boolean
        prefix                    // val: string

    main:
        ch_versions = Channel.empty()
        ch_genome = genome.map { genome -> [ [id:genome.baseName], genome ] }

        // INFO: every iso-* consumer below reads ONE reference annotation. When ANNEVO
        // or TIBERIUS provide extra ones, all are collapsed into a single
        // coordinate-sorted bed here so the modules keep staging a single
        // --reference/--ref path.
        // ponytail: guarded so runs without ANNEVO/TIBERIUS stay bit-identical to before.
        if (params.annevo_annotation || params.annevo_predict || params.tiberius_annotation || params.tiberius_predict) {
            ch_reference_unsorted = annotation
                .map { _meta, bed -> bed }
                .concat(annevo)
                .concat(tiberius)
                .collectFile(
                    name: "${prefix ?: 'polish'}.reference.bed",
                    newLine: false
                )
                .map { bed -> [ [ id: bed.baseName ], bed ] }

            SORT_BED_REFERENCE(ch_reference_unsorted)
            ch_reference = SORT_BED_REFERENCE.out.sorted
            ch_versions = ch_versions.mix(SORT_BED_REFERENCE.out.versions)
        } else {
            ch_reference = annotation
        }

        ch_metassembly = metassembly.branch { meta, f ->
            gtf: f.name.endsWith('.gtf')
            bed: f.name.endsWith('.bed')
        }

        GTF_REMOVE_DIRT ( ch_metassembly.gtf )
        GXF2BED         ( GTF_REMOVE_DIRT.out.gtf )

        // only one branch carries data per run; the empty one spawns no tasks
        ch_cleaned_transcripts = GXF2BED.out.bed.mix( ch_metassembly.bed )

        ISOTOOLS_FUSION(
            ch_cleaned_transcripts,
            ch_reference
        )

        ch_splice_scores = splice_scores_dir ?
            Channel.fromPath(splice_scores_dir, checkIfExists: true)
                .map { it -> [ [ id: it.baseName ], it ] } :
            Channel.of([[], []])

        ch_full_length_transcripts = Channel.empty()
        ch_scraps = Channel.empty()
        if (annotation) {
            ISOTOOLS_ORPHAN(
                ISOTOOLS_FUSION.out.pass,
                ch_reference,
                ch_splice_scores
            )

            ch_full_length_transcripts = ISOTOOLS_ORPHAN.out.hq
            ch_scraps = ISOTOOLS_ORPHAN.out.scraps
        } else {
            ISOTOOLS_ORPHAN_DENOVO(
                GXF2BED.out.bed, // INFO: because fusion depends on annotation
                Channel.of([[], []]),
                ch_splice_scores
            )

            ch_full_length_transcripts = ISOTOOLS_ORPHAN_DENOVO.out.hq
            ch_scraps = ISOTOOLS_ORPHAN_DENOVO.out.scraps
        }

        // INFO: checkpoint to sort all beds
        SORT_BED_FL_TRANSCRIPTS(ch_full_length_transcripts)
        ch_full_length_transcripts = SORT_BED_FL_TRANSCRIPTS.out.sorted

        SORT_BED_SCRAPS(ch_scraps)
        SORT_BED_FUSIONS(ISOTOOLS_FUSION.out.fusion)

        XLOCI_EXTRACT_INTRONS(ch_genome, ch_full_length_transcripts)
        IIC_PREDICT_SPLICEOSOME(XLOCI_EXTRACT_INTRONS.out.tsv)

        if (repeats) {
            Channel.value([
                [ id: "repeats" ],
                file(repeats, checkIfExists: true)
            ]).set { ch_repeats }
        } else {
            ch_repeats = Channel.value([[:], []])
        }

        ch_full_length_transcripts
          .map { meta, bed -> tuple(meta.id, meta, bed) }
          .join(
            IIC_PREDICT_SPLICEOSOME.out.iic
              .map { meta, iic -> tuple(meta.id, meta, iic) }
          )
          .map { id, read_meta, bed, iic_meta, iic ->
            tuple(read_meta, bed, iic)
          }
          .set { ch_classify_inputs }

        ISOTOOLS_CLASSIFY_INTRON(
          ch_classify_inputs,
          ch_genome,
          ch_reference,
          ch_repeats,
          ch_splice_scores
        )

        ISOTOOLS_INTRON_RETENTION(
          ch_full_length_transcripts,
          ISOTOOLS_CLASSIFY_INTRON.out.tsv
        )

        STRIP_RETENTIONS(
            ch_full_length_transcripts,
            ISOTOOLS_INTRON_RETENTION.out.descriptor,
            "RETENTION"
        )
        STRIP_STRONG_RTS(
            STRIP_RETENTIONS.out.hq
                .map { meta, bed -> [ [ id: meta.id + '_no_retentions' ], bed ] },
            ISOTOOLS_INTRON_RETENTION.out.descriptor,
            "HAS_STRONG_RT"
        )
        STRIP_WEAK_RTS(
            STRIP_STRONG_RTS.out.hq
                .map { meta, bed -> [ [ id: meta.id + '_no_strong_rts' ], bed ] },
            ISOTOOLS_INTRON_RETENTION.out.descriptor,
            "HAS_WEAK_RT"
        )
        STRIP_ARTIFACTS(
            STRIP_WEAK_RTS.out.hq
                .map { meta, bed -> [ [ id: meta.id + '_no_weak_rts' ], bed ] },
            ISOTOOLS_INTRON_RETENTION.out.descriptor,
            "HAS_ARTIFACT"
        )

        STRIP_ARTIFACTS.out.hq
            .map { meta, bed -> [ [ id: prefix + '_flnc' ], bed ] }
            .set { ch_fl_hq_transcripts }

        // Strip keeps whatever order it was given. bigtools bedtobigbed
        // requires one contiguous run per chromosome, in raw byte order.
        SORT_BED_STRONG_RTS(STRIP_STRONG_RTS.out.hq)
        SORT_BED_WEAK_RTS(STRIP_WEAK_RTS.out.hq)
        SORT_BED_ARTIFACTS(STRIP_ARTIFACTS.out.hq)

        ch_xorf_hq = Channel.empty()
        if (params.xorf_call_orfs) {
            XORF_RUN(
                ch_fl_hq_transcripts,
                ch_genome,
                'xorf'
            )
            // With selenocysteine masking the merged id is `<t>_flnc@<masked|UNMSK>`
            // and which group wins the submodule's collect is order-dependent, so
            // canonicalize to the deterministic `<t>_flnc` before it keys filenames.
            // Upstream RENAME now emits only *.renamed.bed, so bed is a single file.
            ch_xorf_hq = XORF_RUN.out.files.map { meta, bed, tsv ->
                if (bed instanceof List) {
                    error "XORF emitted List bed for ${meta.id}: ${bed} — expected single file " +
                          "(check params.xorf_skip_joined_concat / rename module)"
                }
                [ meta + [ id: meta.id.tokenize('@')[0] ], bed ]
            }
            ch_versions = ch_versions.mix(XORF_RUN.out.versions)
        }

        ch_final_hq = Channel.empty()
        ch_final_truncations = Channel.empty()
        ch_final_retentions = STRIP_RETENTIONS.out.discard

        if (do_twopass_polish) {
            // Fail fast if twopass requested but ORF channel is empty.
            // main.nf validates FULL_RUN, but FROM_POLISHING bypasses it and XORF can
            // still emit empty if all ORFs filtered — without this, POLISH_TWOPASS
            // silently yields empty or retention-only output.
            // Skip guard in stubMode where channels are intentionally empty.
            ch_xorf_hq_for_twopass = workflow.stubRun ?
                ch_xorf_hq :
                ch_xorf_hq.ifEmpty {
                    error "POLISH_TWOPASS requires XORF ORF predictions but ch_xorf_hq is empty. " +
                          "Ensure params.xorf_call_orfs=true and XORF produced output " +
                          "(check xorf filters / truncation / selenocysteine masking)."
                }
            // POLISH_TWOPASS joins its inputs on the sample prefix (XORF ids are
            // `<prefix>_flnc`, so `source_id` strips back to `<prefix>`), but the
            // first-pass IIC and retention channels carry the metassembly basename.
            // In `--from polish` that differs from `prefix` (e.g. `HLmyoMyot9A.clean`
            // vs `HLmyoMyot9A`), so the joins collapse to zero records and the
            // classify/retention/strip tasks are silently skipped. Normalize both to
            // the prefix key here; no-op on the full path where basename == prefix.
            POLISH_TWOPASS(
                ch_xorf_hq_for_twopass,
                XORF_RUN.out.truncations,
                STRIP_RETENTIONS.out.discard
                    .map { meta, bed -> [ meta + [ id: "${prefix}.discard" ], bed ] },
                IIC_PREDICT_SPLICEOSOME.out.iic
                    .map { meta, iic -> [ meta + [ id: prefix ], iic ] },
                ch_genome,
                ch_reference,
                ch_repeats,
                ch_splice_scores
            )

            ch_final_hq = POLISH_TWOPASS.out.hq
            ch_final_truncations = POLISH_TWOPASS.out.truncations
            ch_final_retentions = POLISH_TWOPASS.out.retentions
            ch_versions = ch_versions.mix(POLISH_TWOPASS.out.versions)
        } else if (params.xorf_call_orfs) {
            // ORF-only path (default: xorf=true, twopass=false) previously ignored
            // XORF and fed first-pass HQ to NMD. Now correctly route ORF predictions
            // to NMD/publish. Sort pre-NMD for determinism (mirrors twopass/flnc sorted).
            ch_xorf_hq_for_nmd = workflow.stubRun ?
                ch_xorf_hq :
                ch_xorf_hq.ifEmpty {
                    error "xorf_call_orfs=true but XORF produced no HQ BED (ch_xorf_hq empty). " +
                          "Check XORF filtering / truncation."
                }
            SORT_BED_XORF(ch_xorf_hq_for_nmd)
            SORT_BED_TRUNCATIONS_XORF(XORF_RUN.out.truncations)
            ch_final_hq = SORT_BED_XORF.out.sorted
            ch_final_truncations = SORT_BED_TRUNCATIONS_XORF.out.sorted
            ch_versions = ch_versions
                .mix(SORT_BED_XORF.out.versions)
                .mix(SORT_BED_TRUNCATIONS_XORF.out.versions)
        } else {
            ch_final_hq = ch_fl_hq_transcripts
        }

        SORT_BED_RETENTIONS(ch_final_retentions)

        // Guard NMD input: fail fast instead of silently publishing nothing
        // Skip in stubRun where dummy channels may be empty.
        ch_final_hq_guarded = workflow.stubRun ?
            ch_final_hq :
            ch_final_hq.ifEmpty {
                error "No HQ transcripts available for NMD (ch_final_hq empty). " +
                      "Check upstream filtering / XORF / twopass."
            }

        ISOTOOLS_NMD(
            ch_final_hq_guarded
        )

        // iso-nmd splits the HQ set: `reads` is the NMD-passing continuation
        // (published as 10_final HQ), `nmd` is the NMD-positive discard
        // (published under 10_final/nmd). Neither glob preserves coordinate
        // order; bedToBigBed requires sorted input. `nmd` is optional — a run
        // with zero NMD hits never fires SORT_BED_NMD.
        SORT_BED_NMD(ISOTOOLS_NMD.out.nmd)

        // An additional renaming step is required to override UNKNOWN tags
        // due to blasting low-confidence ORFs. We perform CDS overlap and
        // derive names from the component. All-UNKNOWN components remain
        // the same.
        RENAME_FINAL_TRANSCRIPTS(ISOTOOLS_NMD.out.reads)
        SORT_BED_FINAL(RENAME_FINAL_TRANSCRIPTS.out.reads)
        ch_final_hq = SORT_BED_FINAL.out.sorted
        PUBLISH_FINAL_TRANSCRIPTS(ch_final_hq)

        ch_versions = ch_versions
            .mix(ISOTOOLS_NMD.out.versions)
            .mix(SORT_BED_FINAL.out.versions)
            .mix(SORT_BED_NMD.out.versions)
            .mix(SORT_BED_STRONG_RTS.out.versions)
            .mix(SORT_BED_WEAK_RTS.out.versions)
            .mix(SORT_BED_ARTIFACTS.out.versions)
            .mix(SORT_BED_RETENTIONS.out.versions)
            .mix(RENAME_FINAL_TRANSCRIPTS.out.versions)

    emit:
        hq             = ch_final_hq
        truncations    = ch_final_truncations
        retentions     = SORT_BED_RETENTIONS.out.sorted
        strong_rts     = SORT_BED_STRONG_RTS.out.sorted
        weak_rts       = SORT_BED_WEAK_RTS.out.sorted
        artifacts      = SORT_BED_ARTIFACTS.out.sorted
        introns        = ISOTOOLS_CLASSIFY_INTRON.out.tsv
        scraps         = SORT_BED_SCRAPS.out.sorted
        fusions        = SORT_BED_FUSIONS.out.sorted
        nmd            = ISOTOOLS_NMD.out.nmd
        versions       = ch_versions
}
