<p align="center">
  <p align="center">
    <img width=100 align="center" src="../figures/logo.png" >
  </p>

<p align="center">
  <picture>
    <source
      media="(prefers-color-scheme: dark)"
      srcset="../figures/hillerlab-dark.png"
    >
    <source
      media="(prefers-color-scheme: light)"
      srcset="../figures/hillerlab-light.png"
    >
    <img
      width="200"
      alt="Hiller Lab"
      src="../figures/hillerlab-light.png"
    >
  </picture>
</p>

  <span>
    <h1 align="center">
        xasm
    </h1>
  </span>

  <span>
    <h2 align="center">
        CHANGELOG
    </h2>
  </span>

  <p align="center">
    <a href="https://github.com/hillerlab/xasm" reference="_blank">
      <img alt="GitHub License" src="https://img.shields.io/github/license/hillerlab/xasm?color=blue">
    </a>
  </p>

  <p align="center">
    <samp>
        <span> The Hiller Lab at the Senckenberg Research Institute </span>
        <br>
        <br>
        <a href="https://nbisweden.github.io/workshop-RNAseq/2011/lab_assembly.html">metassembly</a> .
        <a href="https://github.com/hillerlab/xasm/blob/main/assets/pipeline/xasm.mermaid">pipeline</a> .
        <a href="https://hillerlab.com/">us</a>
    </samp>
  </p>

</p>

---

All notable changes to xasm will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [0.1.12] - 2026-09-10

### Added

- Strandedness inference on the CBQ path. When CBQ reads exist before trimming (`bqtools_encode_fastqs` or native `.cbq`), `METASSEMBLE` forks them into `bqc sniff strand` in parallel with `BQC` trimming and writes the inferred `forward`/`reverse`/`unstranded` into each sample's `strandedness` metadata. The reference CDS transcriptome is extracted once with `xloci` (`XLOCI_EXTRACT_CDS`), turned into one reusable Salmon 2.x index by `bqc sniff index` (`BQC_SNIFF_INDEX`), and consumed per sample by `BQC_SNIFF_STRAND`; all three publish under `00_prepare/strand/`. `--strand_transcriptome` skips the extraction, `--strand_salmon_index` skips extraction + build, and `--infer_strandedness false` turns the branch off. Requires a `bqc` image built with `--features sniff-strand` and bqc >= v0.0.4.
- The run samplesheet gained a tenth column, the per-sample strandedness (appended after `assembled_count`; previous columns unchanged).
- `PREPARE_INDEXES` now emits a `strand_index` channel alongside `star_index`/`deacon_index`.

### Changed

- `STRINGTIE3` (`--fr`/`--rf`) and `TRANSMETA` (`-s`) now see inferred strandedness instead of the hardcoded `unstranded` for CBQ samples. `ASSEMBLY` warns when a transmeta chunk mixes orientations, since transmeta takes a single `-s` per chunk. ALETSCH is unchanged: its protocol field is not strandedness and it already consumes the `XS` tag STAR emits. `COVERAGE` stays unstranded by design.
- FASTQ/fastp samples and `bqtools_encode_before_alignment` runs cannot be sniffed (no CBQ before trimming) and keep `unstranded`; `validateRun` warns once.
- Version bumped to `0.1.12` in the pipeline manifest.

---

## [0.1.11] - 2026-09-04

### Added

- Optional CPU ESMFold2-Fast pLDDT scores in XORF BLAST via `--xorf_esm` (requires `--xorf_call_orfs`; off by default). Weights are downloaded once per XORF run unless `--xorf_esmfold_local_weights` points at a Hugging Face hub cache. Bumped the `modules/xorf` submodule to upstream v0.1.0.
- Optional TIBERIUS annotation as a further gene model for polishing, wired like ANNEVO. `--tiberius_annotation <file>` takes a ready-made TIBERIUS GTF/GFF3 (or GTF/GFF/BED); `--tiberius_predict` runs TIBERIUS on the genome instead (requires `--tiberius_model_cfg`, a Hiller alias from `/opt/tiberius/models.tsv`). Unlike ANNEVO the two flags may be combined, and TIBERIUS may run alongside ANNEVO — every active source contributes a BED. The TIBERIUS subworkflow and its predict/merge modules are vendored from `hillerlab/core`.
- `POLISH` now merges up to three BED12s (reference + ANNEVO + TIBERIUS, the latter carrying 0–2 elements) into the single coordinate-sorted `<prefix>.reference.sorted.bed` fed to `iso-orphan` (`--ref`), `iso-fusion` (`--ref`), `iso-classify` (`--toga`) and the two-pass `--toga` union. The merge stays guarded: runs without ANNEVO/TIBERIUS hand the annotation channel through untouched and stay bit-identical to before.
- `test-tiberius` e2e profile: reuses the `test` profile with `--tiberius_annotation <fixture>` and a separate output dir, and asserts the conversion + merge tasks ran while prediction did not.

### Changed

- Version bumped to `0.1.11` in the pipeline manifest.

### Notes

- The TIBERIUS branch is container-only (`ghcr.io/hillerlab/tiberius:latest`): no conda environment is shipped, so `-profile conda` cannot run `--tiberius_predict`. The image tag floats — pin it if a run must be bit-reproducible.

---

## [0.1.10] - 2026-09-04

### Fixed

- Documented that `10_final/truncations/` and `10_final/nmd/` are only published when such discards exist: empty categories emit nothing through STRIP → collectFile → SORT_BED by design (verified the non-empty publish path end-to-end with the real selectors). Completed the `10_final/` tree in `usage.md` (`scraps/`, `fusions/`, conditional `truncations/`/`nmd/`) and marked the conditional dirs in the `README.md` tree.

### Changed

- Version bumped to `0.1.10` in the pipeline manifest.

---

## [0.1.9] - 2026-09-04

### Added

- Optional ANNEVO annotation as a second gene model for polishing. `--annevo_annotation <file>` takes a ready-made ANNEVO GFF3 (or GTF/GFF/BED); `--annevo_predict` runs ANNEVO on the genome instead (requires `--annevo_lineage`, ignored when a file is given). The ANNEVO subworkflow, its prediction/decoding modules, and `FXSPLIT` are vendored from `hillerlab/core@9a407a0`.
- `POLISH` merges the ANNEVO BED12 with the reference annotation into a single coordinate-sorted `<prefix>.reference.sorted.bed` and feeds that one file to `iso-orphan` (`--ref`), `iso-fusion` (`--ref`), `iso-classify` (`--toga`) and the two-pass `--toga` union, so an intron supported by either source stays categorized. The merge is guarded by `--annevo_annotation`/`--annevo_predict`: runs without ANNEVO hand the annotation channel through untouched and stay bit-identical to before.
- `test-annevo` e2e profile: reuses the `test` profile with `--annevo_annotation <fixture>` and a separate output dir, and asserts the conversion + merge tasks ran while prediction did not.

### Changed

- Version bumped to `0.1.9` in the pipeline manifest.

### Notes

- The ANNEVO branch is container-only: neither ANNEVO nor `FXSPLIT` ships a conda environment, so `-profile conda` cannot run `--annevo_predict`.

---

## [0.1.8] - 2026-09-02

### Added

- `RENAME_FINAL_TRANSCRIPTS` resolves `@UNKNOWN` protein tags left by low-confidence XORF BLAST hits before publishing. Transcripts are packed into CDS-overlap components with `py-packbed` (`ghcr.io/alejandrogzi/py-packbed`, container-only — the wheel is not on bioconda); within a component carrying `UNKNOWN` tags, every `UNKNOWN` adopts the protein name supported by the most members (ties: alphabetical), all-`UNKNOWN` and `UNKNOWN`-free components are flushed unchanged, and untagged names never vote. The record count is preserved — only names change — and a `10_final/*.renamed.tsv` map (one `old<TAB>new` pair per changed transcript, so `wc -l` is the change count) is published alongside the final BED. The embedded script self-checks every branch via `--selftest` before each run; e2e goldens assert the new task and published map for `test`, `test-sb`, and `test-tm`.
- The `py-packbed` image entrypoints `python3`, which would swallow Nextflow's `/bin/bash -ue` wrapper, so the module passes `--entrypoint /usr/bin/env` for docker/podman engines. Upstream note: dropping `ENTRYPOINT ["python3"]` from `py-packbed.Dockerfile` would make the image Nextflow-friendly out of the box.

### Changed

- Version bumped to `0.1.8` in the pipeline manifest.

---

## [0.1.7] - 2026-08-26

### Fixed

- Twopass `SORT_BED_TWOPASS` no longer starts before XORF truncation stripping. The vendored XORF `files` emit was `CONCAT_RENAMED` even when `xorf_do_polishing` ran detach / iso-utr / strip, so the HQ merge (and then NMD) raced those steps and published unstripped ORFs. XORF `0.0.46` (`hillerlab/xorf#22`) emits stripped HQ as `files` plus `truncations` / `duplicates`; polish now waits on that emit.
- Combined truncation discards from the first-pass XORF run and the twopass-retention XORF run are sorted and published under `10_final/truncations`. NMD is gated on that combine so it cannot start until both HQ merge and truncation union are ready.
- `SORT_BED_NMD` now sorts `ISOTOOLS_NMD.out.nmd` (NMD-positive transcripts, `10_final/nmd`). The NMD-passing HQ continuation is sorted by `SORT_BED_FINAL` and remains the published `10_final` BED. A run with zero NMD hits simply skips `SORT_BED_NMD` (`nmd` is optional).

### Changed

- XORF gitlink bumped to `b0a6364` (`hillerlab/xorf#22`, merged to master).
- Version bumped to `0.1.7` in the pipeline manifest.

---

## [0.1.6] - 2026-08-25

### Fixed

- XORF intermediate outputs are published again. The per-instance publish selectors for `09_polish/xorf/` and `09_polish/xorf_twopass_retentions/` anchored on the wrapper name alone, but xorf module processes carry an extra `:XORF:` workflow segment in their qualified names — the anchored patterns never matched and nothing was published. Selectors now include the segment, and the e2e golden asserts the published `09_polish/xorf/` files (line counts) so a regression fails `verify.py` instead of silently dropping outputs.

### Changed

- Version bumped to `0.1.6` in the pipeline manifest.

---

## [0.1.5] - 2026-08-24

### Added

- Second-pass XORF ORF calling on rescued retention transcripts (`POLISH_TWOPASS`). After twopass re-review rescues retention discards, they now get their own XORF ORF round before merging into the twopass HQ set; the merged ORF predictions feed NMD and the final output.
- XORF outputs from the twopass retention run publish under `09_polish/xorf_twopass_retentions/` (mirrored `00_concat`…`04_results` + `XORF_PIPELINE_INFO`), keeping `09_polish/xorf/` exclusively first-pass. `XORF_RUN` takes an optional output-subdir argument; first-pass callers are unchanged.
- `--allow-missing` for twopass retention re-review (`ISOTOOLS_INTRON_RETENTION_TWOPASS`), warning instead of panicking when read/retention introns are missing from the reference.

### Fixed

- Twopass XORF step wiring: missing `XORF_RUN` alias import, `ch_genome` scoping error (used the `genome` take argument), and a `ch_versions` self-reference before initialization in `POLISH_TWOPASS`. The step previously aborted at runtime whenever `--do_twopass_polish --xorf_call_orfs` was set.
- Twopass join keys now normalize to `prefix` in the `--from polish` checkpoint (`POLISH_TWOPASS` joins XORF hq ids `<prefix>_flnc` against first-pass IIC/retention channels keyed by the metassembly basename, which differs from `prefix` in checkpointed runs).

### Changed

- E2E golden task counts updated for the second XORF instance (`WGET_SAMBA_WEIGHTS`, `GENOMEMASK_SELENO` run once per instance).
- Version bumped to `0.1.5` in the pipeline manifest.

---

## [0.1.4] - 2026-08-21

### Fixed

- XORF ORF predictions now reach the final output on the default path (`--xorf_call_orfs true`, `--do_twopass_polish false`). `POLISH` previously assigned `ch_final_hq = ch_fl_hq_transcripts` and only overwrote it inside the twopass branch, so `ISOTOOLS_NMD` and `PUBLISH_FINAL_TRANSCRIPTS` received first-pass HQ transcripts and the entire XORF run was silently discarded (only its versions entry survived). The ORF-only branch now sorts the ORF BED (`SORT_BED_XORF`) and feeds it to NMD/publish; the pre-ORF transcripts are used only when ORF calling is disabled.
- Empty-channel failures no longer pass silently. If XORF produces no predictions while twopass polish is requested — or when there is nothing to feed NMD — the run aborts with an explicit message instead of `POLISH_TWOPASS` degrading to empty or retention-only output and publishing nothing. Guards are skipped under `-stub` runs.
- `RENAME_PREDICTIONS` emitted both the staged input BED and `${prefix}.renamed.bed` (its `path("*.bed")` glob matched the staged input), so `out.files` carried a two-element list per group and `CONCAT_RENAMED` concatenated both copies, duplicating every record. The module now emits only `*.renamed.bed` / `*.renamed.tsv`, and `POLISH` errors if a list ever reaches it instead of silently keeping `bed[0]` (which could pick the pre-rename file).
- The `polish` checkpoint validates again: `FROM_POLISHING` calls `validateFromPolishing()`, which now rejects `--do_twopass_polish` without `--xorf_call_orfs` exactly like a full run.
- `ISOTOOLS_NMD` stub creates the `nmd/` directory and files matching its output globs (`touch nmd/*` failed on stub runs).

### Changed

- `modules/xorf` bumped with the rename-output fix above (submodule commit `96af8dd`).
- `subworkflows/polish/main.nf`: long lines wrapped, mixed indentation normalized; no functional changes.
- Version bumped to `0.1.4` in the pipeline manifest.

---

## [0.1.3] - 2026-08-17

### Added

- STAR can consume the CBQ/bqc path. `bqtools_encode_fastqs` or native `.cbq` inputs now run `bqc → deacon → bqtools decode → STAR` instead of requiring `--aligner ruSTAR`. New `BQTOOLS_DECODE` module expands a paired `.cbq` back to `${id}_1.fastq.gz` / `${id}_2.fastq.gz`. `bqtools_encode_before_alignment` remains ruSTAR-only.
- `test-bqc-star` end-to-end profile exercises encode → bqc → deacon → decode → STAR on the existing fixture and is checked against the same samplesheet numbers as `test` (pipeline-integrity CI and the nightly matrix).

### Changed

- `bqc_adapter_auto_detect` now defaults to `true` (`--auto-detect`). `bqc` still aborts when adapter evidence is ambiguous; turn it off or pass explicit `bqc_adapter_r1` / `bqc_adapter_r2`. The e2e `test-rustar` and `test-bqc-star` profiles set it back to `false` because the synthetic fixture looks like a pooled library.
- `BQTOOLS_ENCODE` (and the new decode module) pull `ghcr.io/hillerlab/bqtools:latest` instead of a digest pin.
- Two-pass `--toga` is now the concatenation of the reference annotation BED and the XORF ORF BED (was ORF-only since 0.1.2). An intron is supported in the second pass if either source has it. The `test` golden now expects `chrTestA: 2 / chrTestB: 2` final transcripts and that the fixture's first-pass retention discard is kept (`10_final/retentions/test_twopass.discard.bed`).
- Version bumped to `0.1.3` in the pipeline manifest.

---

## [0.1.2] - 2026-08-14

### Fixed

- Two-pass polishing failed with a process "input file name collision" when `do_twopass_polish` is enabled: `POLISH_TWOPASS` passed the XORF ORF predictions (`hq`) as the `--toga` reference annotation, and the bed was staged as both the reads and the annotation input of `ISOTOOLS_CLASSIFY_INTRON_TWOPASS`. The ORF set is now collapsed to a single merged annotation bed (`<prefix>.hq_annotation.bed`) that is replicated to one element per transcript — nextflow pairs multi-channel process inputs pairwise, so a single merged element would truncate classification to one transcript.
- `-stub-run` failed in the polish path on several process stubs that did not produce their declared outputs: `ISOTOOLS_FUSION` (`touch ${prefix}/*` without creating the directory, plus a stray `_${meta.chr}` suffix), `ISOTOOLS_ORPHAN` (undefined `prefix` and an undeclared `meta_scraps` output value), `SORT_BED` (`touch *.bed` re-touched the input, never creating `<prefix>.sorted.bed`), `STRIP_OCCURRENCES` (same pattern for `*.striped.bed`/`*.discard.bed`, plus an undeclared `meta2` output value), `ISOTOOLS_INTRON_RETENTION` (`touch *.tsv` only re-touched the staged input tsv, which Nextflow never binds as an output, so the descriptor channel stayed empty), and `BEDGRAPHTOBIGWIG` (placeholder `.BEDGRAPHTOBIGWIG` name while the output declares `*.bw`). Each stub now uses its script prefix and creates the declared output files, so a stub run executes the full metassembly→polish first-pass graph. (The vendored XORF gitlink's stubs still lack several outputs in stub mode.)

### Changed

- The two-pass `--toga` reference is now the merged XORF ORF set (`<prefix>.hq_annotation.bed`) instead of the reference annotation, so introns supported by ORF evidence in the query set are not re-categorized as unsupported during second-pass reclassification. This is consistent with the first pass operating on ORF-validated transcripts (see 0.1.1). The e2e golden was updated accordingly: the `test` profile now expects `chrTestA: 3 / chrTestB: 2` final transcripts and that `10_final/retentions/test_twopass.discard.bed` is not emitted (the fixture's retention candidate is rescued into the final set).
- Version bumped to `0.1.2` in the pipeline manifest.

---

## [0.1.1] - 2026-08-10

### Added

- ORF calling with XORF, gated by `xorf_call_orfs` (default `true`). The pinned `xorf` submodule is wired in as a gitlink and its ORF subworkflow is called through a new `XORF_RUN` wrapper (`subworkflows/xorf/main.nf`) that adapts channel shapes and prepares the protein database: a custom `xorf_custom_database` (`.dmnd`/`.dmnd.gz` used as-is, `.fa`/`.fasta` appended to SwissProt and reindexed with diamond) or the default zenodo download. The full ORF chain runs with per-step parameters (`xorf_chunk_size`, `xorf_skip_netstart`, `xorf_predict_min_score_max_predictions`, `xorf_predict_max_predictions`, `xorf_predict_threshold`, `xorf_predict_keep_raw`, `xorf_do_polishing`, `xorf_rename_*`, `xorf_samba_weights`/`xorf_samba_local_weights`), publishing under `09_polish/xorf/{00_concat,01_renamed,02_merged,03_duplicates,04_results}`; the `rename_predictions.py` helper is vendored into `bin/` and `assets/scripts/xorf.sh` manages the submodule.
- Two-pass polishing now operates on ORF-validated transcripts: when `do_twopass_polish` is set, `xorf_call_orfs` is required (validated at startup) and the twopass HQ input is XORF's ORF predictions rather than the raw first-pass HQ, so intron reclassification and retention rescuing run against transcripts with ORF evidence.
- Nonsense-mediated decay filtering: new `ISOTOOLS_NMD` module runs `iso-nmd` (premature-termination-codon detection) on the final HQ set; NMD candidates publish to `10_final/nmd` and the published final transcript set is the NMD-passing reads.

### Changed

- Output layout: single-exon scraps are published under `10_final/scraps` (previously `09_polish`) and fusion calls under `10_final/fusions` (previously `09_polish/fusions`), grouping every final transcript category under `10_final`.
- XORF input transcripts are named `<prefix>_flnc`; the deterministic id flows through XORF outputs and two-pass artifact names (e.g. `09_polish/twopass/classify/<prefix>_flnc@<prefix>_flnc.reference_introns.tsv`).
- Repository layout: the end-to-end fixture moved from `test_data/` to `assets/test/test_data/` and the test harness from `tests/` to `assets/ci/`; all `nextflow.config` profiles, CI workflows, `run.sh`/`verify.py`, the fixture builder, and docs were updated accordingly. The e2e golden suite now also asserts the XORF chain (fixture protein database, ORF task counts, `09_polish/xorf/` artifacts).
- The `INTRONIC` process label was raised from `process_medium` to `process_high`.
- Version bumped to `0.1.1` in the pipeline manifest.

### Removed

- The `publish.yml` GitHub Actions workflow (automatic release publishing) — no longer needed.

---

## [0.1.0] - 2026-08-07

### Added

- ruSTAR as an alternative aligner (`aligner = "ruSTAR"`): a Rust reimplementation of STAR with the same CLI flags and STAR-compatible log output, reading both FASTQ and CBQ natively. Because the aligner image ships no `samtools`, BAM indexing is delegated to the new `SAMTOOLS_INDEX` module. The two-pass scheme and the pre-curated-junction shortcut apply as in the STAR path.
- CBQ (columnar BINSEQ) read support. Native `.cbq` files in `input_dir` are detected automatically; `bqtools_encode_fastqs` encodes FASTQ inputs before QC (QC then runs via `bqc` instead of `fastp`), and `bqtools_encode_before_alignment` keeps the `fastp`+deacon path and encodes to `.cbq` only for the aligner (the two are mutually exclusive, configurable in `params.json` and `main.nf` help).
- `BQTOOLS_ENCODE` module to collapse a paired FASTQ set into a single `.cbq`.
- `BQC` module: CBQ-native all-in-one QC (adapter, trimming, filtering) with a structured JSON report and `bqc_*` parameters. Paired counts are doubled so `fastp_min_trimmed_reads` means the same thing on both paths.
- deacon upgraded to 0.16.0, which emits `.cbq` on the CBQ decontamination path.
- `test-rustar` end-to-end profile exercising the whole CBQ path (`bqtools_encode_fastqs` → `bqc` → deacon → ruSTAR → `samtools` index → cleanup) on the existing fixture, wired into the nightly matrix and CI profile validation, with a smoke checker that asserts the CBQ processes run and the STAR/FASTP path does not.

### Changed

- New `aligner` parameter (`STAR` | `ruSTAR`, default `STAR`) separating alignment configuration; version bumped to `0.1.0` in the pipeline manifest.
- Coverage tracks are refused up front under ruSTAR (`rustar-aligner` rejects `--outWigStrand Unstranded`), so `aligner = ruSTAR` requires `star_make_coverage = false`.
- Alignment under ruSTAR rebuilds the `(meta, bam, bai)` tuple via `SAMTOOLS_INDEX` instead of treating STAR-style inline indexing as the norm.

### Fixed

- Aletsch per-chromosome cleanup race: with `assembly_by_chr = true` the shared full BAM is read by every per-chromosome local-assembly task, and the `!aletsch_keep_bam && !star_make_coverage` cleanup deleted it after whichever chromosome finished first. The in-task deletion was removed and `REMOVE_BAMS` now waits for the last local assembly before deleting, so `--star_make_coverage false` is safe without `--aletsch_keep_bam`.
- `REMOVE_BAMS` referenced a `biocontainers/bash` image that does not exist on the configured registry; it now uses a valid ubuntu base.
- CHROMSIZE stub failed in `-stub-run`: the stub fell back to the always-empty `meta.id` and created a file where the real `chromsize -o` creates a directory. It now mirrors the script (`genome.baseName`, `mkdir -p`).

---

## [0.0.20] - 2026-08-04

### Added

- Per-chromosome assembly (`assembly_by_chr`, enabled by default). Aligned BAMs are now split by reference sequence with the new `BAMSPLIT_CHROM` module, local assembly and metassembly run per chromosome, and the results are merged into a single metassembly GTF with globally unique transcript IDs. The `assembly_exclude_chromosomes` parameter lets you keep chromosomes in the alignment while omitting them from the assembly fan-out.
- StringTie3 as an alternative local assembler (`assembler = "stringtie3"`), including nascent RNA options (`stringtie3_enable_nascent_assembly`, `stringtie3_include_nascent_rna`).
- TransMeta as an alternative meta-assembler (`metassembler = "transmeta"`), with annotation-guided merging controlled by `transmeta_use_annotation`.
- Two-pass polishing (`do_twopass_polish`). Retention discards from the first pass are reclassified against the first-pass intronIC evidence, UTR retentions are ignored in the second pass, and rescued transcripts are merged back into the clean HQ set. Output lands under `09_polish/twopass/`.
- `star_twopass_junctions_file` parameter to provide a pre-curated junction file (STAR sjdb 4-column format), which skips the STAR first pass and the junction-merging step entirely.
- GFF support: annotations in `.gff` format are now accepted alongside GTF and BED.
- Aletsch insert-size fallbacks (`aletsch_fallback_insert_size`, `aletsch_fallback_insert_std`) that replace a zero profiled insert-size mean when compact per-chromosome inputs fall below Aletsch's internal sample minimum.
- `--intron-track` output for the iso-classify process, emitted during the 0.0.19b milestone.
- End-to-end test suite: a committed 256 KB synthetic fixture (`test_data/e2e/`) with reads, genome, annotation, splice scores, and repeats; `tests/e2e/run.sh` and `tests/e2e/verify.py` with golden counts; new `test`, `test-sb`, and `test-tm` profiles (the old small fixture is now `test-polish`); and a nightly end-to-end workflow (`e2e-nightly.yml`) alongside the expanded CI smoke test.
- Pinned container image digests for the beaver, bed2gtf, chromsize, genepred, isox-py, xloci, and bamsplit modules to make runs reproducible.
- New `REMOVE_BAMS` and `RENAME_BAM` housekeeping modules, and a README logo via `.gitattributes`.

### Changed

- Version bumped to `0.0.20` in the pipeline manifest, which now also lists the authors.
- Beaver and TransMeta outputs are prefixed per chromosome (`<prefix>_<chr>` when `assembly_by_chr` is on) so per-chromosome results can be traced back to their source sample.
- Aletsch now derives its library type from `single_end`/`paired_end` metadata rather than the strandedness field, and supports a per-chromosome `-l` flag. Transcript counts are computed with `awk` instead of `grep -w`.
- The polished GTF/BED conversion path now also handles GFF input in `prepare_indexes` and the `--from polish` checkpoint.
- STAR alignment flows straight to the second pass with a pre-curated junction file when `star_twopass_junctions_file` is set.
- `params.json` was reorganized with a dedicated assembly section and a new polishing section (`from`, `assembly_by_chr`, `assembler`, `metassembler`, `do_twopass_polish`, etc.).
- CI pins Nextflow `25.10.0`, validates every test profile, and runs the default end-to-end smoke test on each push.

### Fixed

- Aletsch runs on single-chromosome chunked input: when a BAM contains only one populated chromosome, the process now selects that chromosome from the original full BAM instead of consuming the split chunk, which Aletsch 1.1.x would otherwise fail to flush.
- RENAME_GTF no longer deletes the input GTF (or its realpath target) after renaming.
- Stub behavior of the Aletsch module now produces the expected GTF output shape.
- Default strandedness metadata is reported as `unstranded` rather than `paired_end`.

---

## [0.0.19] - 2026-07-23

### Added

- New `--from bigbed` checkpoint that allows resuming the pipeline directly from the BigBed conversion step, skipping both meta-assembly and polishing entirely. This is useful when BED files are already available from a prior run and only the BigBed conversion (or re-conversion with different settings) is needed.
- `SORT_BED` module (`modules/custom/sort/main.nf`) that performs coordinate-based sorting (`sort -k1,1 -k2,2n -k3,3n`) on BED files prior to BigBed conversion. Sorting is now applied to full-length transcripts, scraps, and fusions within the polishing subworkflow, ensuring BED input compatibility with `bedToBigBed`.
- BigBed conversion step added to the `--from polish` checkpoint. Previously, BigBed conversion only ran during the full pipeline; now users resuming from polishing also get BigBed output for all categories (HQ, retentions, strong RTs, weak RTs, artifacts, fusions, scraps).
- `prefix` parameter (`params.prefix`) to override the output filename prefix used by the Beaver process. When set, Beaver output files will use this prefix instead of the default.
- `test` profile in `nextflow.config` for running xasm against a minimal synthetic dataset (`test_data/`). The profile configures a single-chromosome test genome, a test annotation GTF, and a test metassembly GTF, with BigBed conversion disabled by default.
- Test data files: `test_data/test_genome.fa`, `test_data/test_annotation.gtf`, and `test_data/test_metassembly.gtf`.
- `.gitignore` entries to preserve `test_data/` and its contents while ignoring other test-related files.

### Changed

- Disabled automatic work directory cleanup (`cleanup = false`) to preserve intermediate files for debugging and pipeline resumption.
- Updated `workflow.onComplete` handler to correctly distinguish between BED and BigBed output directories (`10_final` vs `11_bbs`) when reporting pipeline results, and to handle the case where BigBed conversion is active.
- Updated `params.json` and the task handler script (`assets/scripts/xasm.sh`) to expose the new `all_bed_path` parameter for the `--from bigbed` checkpoint and the `prefix` parameter.
- PublishDir pattern for BigTools processes corrected from `*.bed` to `*.bb` in `nextflow.config`.
- Added `FROM_BIGBED` publishDir configuration to write BigBed output to `11_bbs`.

### Fixed

- Groovy syntax errors in the `onComplete` handler where `and`/`or` keywords were used instead of the correct `&&`/`||` operators, which would cause runtime failures on pipeline completion.
- Fixed the `onComplete` handler referencing a non-existent `10_results` directory instead of the correct `10_final` directory.

---

## [0.0.18a] - 2026-06-16

### Added

- BigBed conversion (`bedToBigBed`) now runs by default for all transcript categories (HQ, retentions, strong RTs, weak RTs, artifacts, fusions, scraps) during the full pipeline run.
- New `11_bbs` output directory for storing BigBed (`.bb`) files alongside the existing `10_final` BED output.
- New `SORT_BED` instances for full-length transcripts, scraps, and fusions in the polishing subworkflow to ensure sorted BED input for downstream tools.

### Changed

- Updated the shell task handler (`assets/scripts/xasm.sh`) to assert correct input format.
- Changed `process_medium` resource label for the intron process.

---

## [0.0.18] - 2026-06-10

### Added

- Retained intronic transcript (RT) and artifact detection and classification implementation in the polishing subworkflow.
- Explicit publish step for classified transcript categories.

### Changed

- Adjusted resource allocations for multi-threaded processes to improve stability under concurrent execution.

### Fixed

- Retention transcripts now emit sorted (striped) output correctly.
- xloci intron extraction now applies `--unmask` as intended.
- Array handling fixes for GTF/BED input channels in the polishing subworkflow.
- Versioning channel propagation to `--from` checkpoints.
- Error strategy updated to `retry` for transient failures in external tool processes.

---

## [0.0.17] - 2026-05-20

### Added

- `--from polish` checkpoint for resuming the pipeline from the polishing step, skipping meta-assembly.
- JSON-based parameter input support.

### Changed

- Major refactor of the pipeline workflow structure and subworkflow organization.

---

## [0.0.16] - 2026-05-10

### Added

- Detach and strip operations for transcript classification.
- isotoools UTR processing integration.
- Array-based input support for select processes.

### Changed

- Reverted isotoools UTR behavior to previous state while retaining the detach/strip infrastructure.

### Fixed

- Compliance with `nf-core` linting standards.

---

## [0.0.15] - 2026-04-28

### Added

- GTF and BED file acceptance as input formats.
- Fusion transcript output as a default pipeline output.
- Pipeline DAG visualization output.
- Gene prediction linting step.
- Copyright headers.
