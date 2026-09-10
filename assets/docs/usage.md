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
        USER GUIDE
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

This guide explains how to run xasm, what every parameter does, and what you get
back in the results folder. If something is confusing,
let us know in a [GitHub issue](https://github.com/alejandrogzi/xasm/issues).

> [!IMPORTANT]
> - **What is xasm?** An end-to-end pipeline that assembles RNA-seq reads into
>  a polished transcriptome. 
> - **How does it work?** xasm takes the RNA-seq reads of **many samples of the
>  same species**, locally aligns them, assembles transcripts per sample, and
>  then **meta-assembles** those transcript sets into a single consensus
>  transcriptome.

> **Requirements:** Nextflow ≥ 24.10.5, Java, and a container engine
> (Docker or Apptainer/Singularity are the supported options). Conda/mamba
> profiles also exist but are not the recommended path.

---

## 1. What xasm does

xasm takes the RNA-seq reads of **many samples of the same species** and
builds one shared, polished transcriptome assembly from them.

The idea: no single sample tells the whole story. Each sample only covers a
part of the expressed genome, so xasm aligns every sample, assembles
transcripts per sample, and then **meta-assembles** those transcript sets into
a single consensus transcriptome. Finally it polishes that consensus by
classifying every intron and removing transcripts that look like artifacts
(read-throughs, retained introns, assembly junk), and gives you UCSC
genome-browser tracks of the result.

```
[ QC + trim ] → [ decontaminate ] → [ STAR 2-pass align ] → [ assemble per sample ]
→ [ meta-assemble ] → [ polish ] → [ final transcripts + browser tracks ]
```

---

## 2. Quick start

```bash
git clone https://github.com/hillerlab/xasm.git
cd xasm
```

Copy `params.json`, fill in your paths, then run:

```bash
# Docker
nextflow run main.nf -params-file params.json -profile docker

# Apptainer / Singularity
nextflow run main.nf -params-file params.json -profile apptainer
```

Everything else is optional. `params.json` carries the **scientific** settings
(genome paths, assembly choices). `nextflow.config` carries the
**infrastructure** settings (memory, CPUs, cluster behaviour) — you normally
never edit that file.

### Smoke-testing your setup first

The repo ships tiny test fixtures so you can check that your container
engine and profiles work before touching real data:

```bash
assets/ci/e2e/run.sh test          # default Aletsch + Beaver
assets/ci/e2e/run.sh test-sb       # StringTie 3 + Beaver
assets/ci/e2e/run.sh test-tm       # TransMeta
assets/ci/e2e/run.sh test-bqc-star # bqc + decode + STAR
assets/ci/e2e/run.sh all           # the whole matrix
```

Set `TEST_ENGINE=apptainer` to switch container engines. Results land in
separate directories under `test_results/`.

---

## 3. Your inputs

| Input | Format | Needed? | What it is used for |
|-------|--------|---------|---------------------|
| `--input_dir` | directory of `.f*q.gz` files | **yes** | The raw reads. Paired files must share a prefix and end in `_1`/`_2` (e.g. `sampleA_1.fastq.gz` + `sampleA_2.fastq.gz`). Single-end reads work too. The sample id is the prefix. |
| `--genome` | `.fa`, `.fasta`, `.2bit`, or `.fa.gz` | **yes** | The reference genome. `.2bit` and gzipped files are converted automatically. Used for STAR index, Deacon index, and all downstream steps. |
| `--annotation` | `.gtf`, `.gff`, `.bed`, or `.gz` | **yes** | The reference annotation. It guides the STAR index, iso-orphan detection, and TransMeta. |
| `--splice_scores_dir` | directory of SpliceAI `.bw` tracks | recommended | Pre-computed SpliceAI scores (donor/acceptor, plus/minus). The polish step uses them to tell real splice sites from artifacts. You can skip it, but the polish step is much weaker without it. |
| `--repeats` | `.bed` / `.gff` / `.gtf` | recommended | Repeats annotation, used to mask repeat-derived artifacts during intron classification. |
| `--annevo_annotation` | ANNEVO `.gff3` (or GTF/GFF/BED) | no | A second, *ab initio* gene model, merged with `--annotation` for polishing. See [4.10a](#410a-annevo-annotation-optional). |
| `--annevo_predict` | flag | no | Run ANNEVO on the genome instead of supplying a file. Requires `--annevo_lineage`. |
| `--tiberius_annotation` | TIBERIUS `.gff3` (or GTF/GFF/BED) | no | A third, *ab initio* gene model, merged with `--annotation` for polishing. See [4.10b](#410b-tiberius-annotation-optional). |
| `--tiberius_predict` | flag | no | Run TIBERIUS on the genome. Requires `--tiberius_model_cfg`. May be combined with `--tiberius_annotation`. |
| `--output_dir` | path | **yes** (defaults to `./results`) | Where results go. |
| `--prefix` | string | no | Prefix for output file names (e.g. `hg38` → `hg38.gtf`). If unset, `metassembly` is used. |

> A "full" run therefore needs `input_dir`, `genome`, `annotation`, and
> `output_dir`. The pipeline refuses to start without the first three.

---

## 4. Parameter reference

All parameters can be passed on the command line (`--param value`) or, better,
inside a JSON file (`-params-file params.json`). Defaults below come from
`nextflow.config`; the shipped `params.json` overrides a few of them — those
are marked **\***.

### 4.1 Inputs and outputs

| Parameter | Default | What it does |
|-----------|---------|--------------|
| `input_dir` | – | Directory containing the fastq files (`*{1,2}.f*q.gz`). **Required.** |
| `output_dir` | `./results` | Where every result directory is written. |
| `genome` | – | Reference genome: fasta, 2bit, or gzipped fasta. **Required.** |
| `annotation` | – | Reference annotation: GTF/GFF/BED (optionally gzipped). **Required.** |
| `splice_scores_dir` | – | Directory of SpliceAI bigwig tracks used during polishing. |
| `repeats` | – | Repeats file (BED/GFF/GTF) used to mask repeat-derived introns. |
| `prefix` | – | Prefix for the final assembly file names. |
| `env` | – | Reserved, currently unused. Ignore it. |

### 4.2 Checkpoints (resume mid-pipeline)

| Parameter | Default | What it does |
|-----------|---------|--------------|
| `from` | – | Resume from a checkpoint instead of the full pipeline. Options: `polish` or `bigbed`. See [section 7](#7-checkpoints-resuming-from-the-middle). |
| `polish_path` | – | The metassembly GTF/BED to polish when `from = "polish"`. |
| `all_bed_path` | – | Directory of final BEDs to convert to BigBed when `from = "bigbed"`. |

### 4.3 Linting

fastq lint runs `fq lint` sanity checks on the reads and fails fast if a
file looks corrupt. It costs a full pass over the data, so it is skipped after the
steps where fastp/deacon already guarantee correctness.

| Parameter | Default | What it does |
|-----------|---------|--------------|
| `fq_skip_linting_at_start` | `false` | Skip linting of the raw input fastqs. Turn on to save time when your reads are known-good. |
| `fq_skip_linting_after_trimming` | `true` | Skip linting of the trimmed reads. |
| `fq_skip_linting_after_deacon` | `true` | Skip linting of the decontaminated reads. |
| `fq_extra_args` | `--disable-validator P001` | Reserved; declared but not consumed by any step. |

### 4.4 Trimming (fastp)

Every sample is trimmed with fastp. This is also where a **sample filter** is
applied: if a sample ends up with fewer than `fastp_min_trimmed_reads` reads
after trimming, it is dropped from the rest of the run (with a warning) —
broken or tiny samples are not worth assembling.

| Parameter | Default | What it does |
|-----------|---------|--------------|
| `fastp_allow_gap_overlap_trimming` | `true` | Enable fastp's overlap-based trimming for paired reads (reads that overlap each other get trimmed consistently). |
| `fastp_allow_correction` | `true` | Enable fastp's read-overlap error correction (fixes mismatches where mates overlap). |
| `fastp_discard_trimmed_pass` | `false` | Drop reads that survive trimming instead of keeping them (useful only for special experiments; keep `false`). |
| `fastp_save_trimmed_fail` | `false` | Keep the reads that *failed* trimming. |
| `fastp_save_merged` | `false` | Keep the merged read pair output when read pairs fully overlap. |
| `fastp_min_trimmed_reads` | `10000` | Minimum number of reads a sample must have **after** trimming. Samples below this are discarded. Lower this (e.g. `100`) for test data. |
| `fastp_deduplicate_reads` | `true` | Deduplicate PCR/optical duplicates (`--dedup`). Huge memory cost on high-depth data — some users disable it. |
| `fastp_deduplication_accuracy` | `3` | Accuracy (1–6) of the duplicate detection; higher = slower but more exact. |

### 4.5 Decontamination (Deacon)

Deacon removes reads that look like they come from a **different organism**
than your genome (bacteria, viruses, other contamination). There are three
ways to get the Deacon index it filters against:

1. **Provide your own** (`deacon_index_path`) — fastest, fully reproducible.
2. **Download one** (`deacon_index_path` + `deacon_download_index = true`) —
   the path is treated as a URL.
3. **Build one from your genome** (`deacon_make_single_index = true`) — the
   pipeline indexes your genome and, optionally, subtracts a background of
   common contaminants (FDA-ARGOS + RefSeq viral) so the filter does not
   mistake your own genome for contamination.

| Parameter | Default | What it does |
|-----------|---------|--------------|
| `deacon_index_path` | – | Path (or URL) of an existing Deacon index. |
| `deacon_download_index` | `false` | Treat `deacon_index_path` as a URL and download it. |
| `deacon_make_single_index` | `true`\* | Build a fresh Deacon index from your genome (`false` in `nextflow.config`; the shipped `params.json` enables it). |
| `deacon_single_index_use_background` | `true`\* | When building an index, also subtract a background set of common contaminants (FDA-ARGOS + RefSeq viral). Recommended. |
| `deacon_multi_index_additional_genome_paths` | – | Extra genomes to include when building a multi-species index (currently unreachable with containers). |
| `deacon_keep_fastp_fastq` | `false` | Keep the trimmed fastqs in addition to the decontaminated ones (uses extra disk). |
| `deacon_entropy_threshold` | `0.5` | Sequence-complexity threshold for the index (lower = filter more aggressively). The background index is always built at 0.5. |
| `deacon_background_download_url` | Zenodo URL | Where to download the contaminant background index from. Only used when building an index *with* background. |

### 4.6 Alignment (STAR, two-pass)

xasm aligns every sample twice: the **first pass** finds splice junctions in
the data itself, the junctions from all samples are merged and filtered, and
the **second pass** re-aligns using those junctions. This "STAR-2pass" scheme
greatly improves detection of novel splice sites.

| Parameter | Default | What it does |
|-----------|---------|--------------|
| `aligner` | `STAR` | Which aligner runs the two-pass scheme. `STAR` reads gzipped FASTQ only (CBQ is decoded after deacon); `ruSTAR` (a Rust reimplementation) reads FASTQ **and** CBQ natively. See [4.6a](#46a-cbq-reads-bqtools--bqc). |
| `star_index_path` | – | Path to an existing STAR index (or `.tar.gz` of one). If unset, the index is built from your genome. |
| `star_ignore_gtf_for_index` | `false` | Build the STAR index **without** the annotation (`--sjdbGTFfile`). |
| `star_ignore_gtf_for_mapping` | `true` | Ignore the annotation during alignment, relying on the junctions found in pass one. This is the recommended 2-pass behaviour. |
| `star_seq_center` | `custom` | Sequencing center recorded in the BAM header. |
| `star_seq_platform` | `illumina` | Sequencing platform recorded in the BAM header. |
| `star_seq_library` | `default` | Library type recorded in the BAM header. |
| `star_machine_type` | `default` | Machine type recorded in the BAM header. |
| `star_min_mapped_reads` | `5` | Reserved; declared but not consumed by any step. |
| `star_extra_align_args` | – | Any extra STAR arguments you want appended to the alignment command. |
| `star_save_unaligned` | `false` | Reserved; declared but not consumed by any step. |
| `star_keep_first_pass_bam` | `false` | Keep the BAM of the first alignment pass (extra disk; normally deleted). |
| `star_delete_fastq_after_alignment` | `true` | Delete the aligned fastqs after the second pass to save disk. |
| `star_twopass_junctions_file` | – | Path to a pre-curated junction file (STAR sjdb format, 4 columns). Providing one **skips the first pass** entirely and aligns once with your junctions. |
| `star_make_coverage` | `true` | Produce coverage tracks (per-sample BigWig → merged BigWig in `05_coverage`). Disable to save time/disk if you do not need them. |

### 4.6a CBQ reads (bqtools + bqc)

CBQ is the columnar BINSEQ read format. A paired sample collapses into a
**single `.cbq` file** carrying both mates, which on realistic libraries is
roughly **20% smaller than the equivalent gzipped FASTQ pair** — the reason to
use it. `fastp` and `fq lint` cannot read CBQ, so on the CBQ path quality
control is done by `bqc` instead; the choice is derived from the input format
and is not configurable.

All of this is **opt-in and off by default**. STAR cannot read CBQ natively, so
on the STAR path xasm decodes after deacon (`bqc → deacon → decode → STAR`).
`bqtools_encode_before_alignment` still requires `--aligner ruSTAR` — encoding
just before STAR would only be decoded again.

There are three ways to get CBQ, and they differ in *when* the encode happens:

| How | Path | Use when |
|-----|------|----------|
| Put `.cbq` files in `input_dir` | `bqc → deacon → ruSTAR` **or** `bqc → deacon → decode → STAR` | Your reads are already CBQ. Set `--aligner` to match. |
| `bqtools_encode_fastqs = true` | `fastq → cbq → bqc → deacon → ruSTAR` **or** `fastq → cbq → bqc → deacon → decode → STAR` | You want the storage saving across QC and deacon, and (with STAR) `bqc` as the trimmer. |
| `bqtools_encode_before_alignment = true` | `fastq → fastp → deacon → cbq → ruSTAR` | You want to keep the validated fastp+deacon path and only cut the peak disk of the alignment queue. ruSTAR only. |

The two encode options are mutually exclusive. Mixing `.cbq` and FASTQ samples
in one `input_dir` is supported — each sample is routed by its own extension.

| Parameter | Default | What it does |
|-----------|---------|--------------|
| `bqtools_encode_fastqs` | `false` | Encode FASTQ inputs to `.cbq` before QC. Your original files in `input_dir` are never deleted. |
| `bqtools_encode_before_alignment` | `false` | Encode to `.cbq` after decontamination, just for the aligner (**ruSTAR only**). The intermediate deacon FASTQs **are** deleted (they are work-directory files, and nothing downstream removes them otherwise). |
| `bqc_min_length` | `15` | `--min-length`: reject reads shorter than this. Always sent — `bqc` refuses to run with no operation configured. |
| `bqc_quality_tail` | `15` | `--quality-tail`: trim 3' bases until the trailing window reaches this Phred score. |
| `bqc_poly_g` | `true` | `--poly-g`: trim G-rich 3' tails (two-colour chemistry). |
| `bqc_max_n` | – | `--max-n`: reject reads with more than this many ambiguous bases. |
| `bqc_adapter_auto_detect` | `true` | `--auto-detect`: infer adapters from the data. Unlike fastp's detection, `bqc` *aborts the task* when the evidence is ambiguous (pooled libraries, concatenated runs). Turn off with `--bqc_adapter_auto_detect false` or pass explicit sequences below. |
| `bqc_adapter_r1` / `bqc_adapter_r2` | – | Explicit adapter sequences per mate. |
| `bqc_extra_args` | – | Extra arguments passed verbatim to `bqc workflow`. |
| `infer_strandedness` | `true` | Infer each sample's RNA-seq strandedness from its raw `.cbq` before trimming (`bqc sniff strand`). The value (`forward`/`reverse`/`unstranded`) is written into the sample metadata and published as the last column of the run samplesheet. Requires CBQ at the trim point (`bqtools_encode_fastqs` or native `.cbq`); FASTQ/fastp and `bqtools_encode_before_alignment` samples stay `unstranded`. `undetermined` — bqc's answer when mapping evidence is below its gates — also falls back to `unstranded`. |
| `strand_transcriptome` | – | Ready-made transcriptome FASTA for the strand index; skips the XLOCI CDS extraction step. |
| `strand_salmon_index` | – | Ready-made Salmon 2.x index directory; skips both the CDS extraction and the index build. The `00_prepare/strand/bqc_salmon_index/` published by a previous run can be reused here. |

> **Strandedness inference (CBQ path).** When CBQ reads exist before trimming,
> xasm extracts the reference CDS transcriptome with `xloci`, builds one
> reusable Salmon 2.x index with `bqc sniff index`, and runs `bqc sniff strand`
> per sample in parallel with `BQC` trimming. The inferred strandedness flows
> into the assemblers that consume it: `stringtie3` (`--fr`/`--rf`) and
> `transmeta` (`-s`). **Aletsch is unaffected** — its third sample-info field is
> the sequencing protocol, not strandedness, and it reads the `XS` tag STAR
> already emits. Coverage tracks stay unstranded. The CDS FASTA, index and
> per-sample reports are published under `00_prepare/strand/`.
>
> The `bqc` image must be built with the `sniff-strand` feature and be at least
> v0.0.4 (which added `bqc sniff index`), i.e.
> `cargo install bqc --features sniff-strand`.

> Reads-after-trimming counts: `bqc` reports *records* (one paired record = both
> mates) while `fastp` reports *reads*. xasm doubles the paired CBQ counts so
> `fastp_min_trimmed_reads` means the same thing on both paths.

**ruSTAR is not a complete STAR drop-in.** `rustar-aligner 0.2.0` differs from
STAR in three ways that xasm works around, all verified against the binary:

| Divergence | Consequence |
|---|---|
| `--limitSjdbInsertNsj` unimplemented | Omitted from the ruSTAR configuration. |
| `--outSAMattributes` has no `ch` tag | ruSTAR BAMs carry `NH HI AS nM NM MD`; the STAR path additionally tags chimeric alignments with `ch`. |
| `--outWigStrand Unstranded` rejected | **Coverage tracks are not available with ruSTAR.** Setting `aligner = ruSTAR` with `star_make_coverage = true` is refused up front; run with `--star_make_coverage false`. |

An unknown flag *aborts* the task rather than being ignored, so if you set
`star_extra_align_args`, confirm the flags exist with `rustar-aligner --help`.

> **No need to pass `--aletsch_keep_bam true` anymore.** With `assembly_by_chr
> = true`, Aletsch consumes the *shared* full BAM for every chromosome; the
> `!aletsch_keep_bam && !star_make_coverage` cleanup used to delete it after
> whichever chromosome finished first, so the siblings failed with "the file
> does not exist". The cleanup now waits for the last chromosome's local
> assembly before removing the BAM, so `--star_make_coverage false` (which
> ruSTAR requires) is safe without `--aletsch_keep_bam`.

### 4.7 Junction filtering

The junctions found across all samples in pass one are merged into one set;
only junctions that are long enough and supported by enough reads survive.

| Parameter | Default | What it does |
|-----------|---------|--------------|
| `junction_min_read_coverage` | `2` | Minimum number of reads supporting a junction to keep it. |
| `junction_min_junction_length` | `50` | Minimum intron length for a junction to be kept (shorter introns are unreliable). |

### 4.8 Keeping intermediate files

The pipeline deletes most intermediate files to save disk. Flip these to
`true` if you want to keep them for debugging or downstream use.

| Parameter | Default | What it does |
|-----------|---------|--------------|
| `bigtools_keep_bedgraph` | `false` | Keep the per-sample bedGraph coverage files. |
| `wigtools_keep_bigwig` | `false` | Keep the per-sample BigWig files (only the merged median BigWig is kept by default). |
| `wigtobigwig_keep_wig` | `false` | Keep the merged median `.wig` file. |

### 4.9 Assembly

Two decisions live here: **which local assembler** builds transcripts from
each sample's alignments, and **which meta-assembler** merges those into one
consensus. By default assembly runs **per chromosome** — each sample's BAM is
split by chromosome and each chromosome is assembled separately, which keeps
memory and runtime manageable on real genomes.

| Parameter | Default | What it does |
|-----------|---------|--------------|
| `skip_assembly` | `false` | Skip the whole assembly + metassembly stage (run only QC, alignment, coverage). |
| `assembly_by_chr` | `true` | Split BAMs per chromosome and assemble/metassemble per chromosome, then concatenate. Turn off to assemble whole-BAM (needs much more memory). |
| `assembly_exclude_chromosomes` | `[]` | Chromosomes to keep in the alignments but exclude from the per-chr assembly fan-out (e.g. `["chrX"]`). Useful to drop problematic contigs. |
| `assembler` | `aletsch` | The per-sample assembler. Options: `aletsch`, `stringtie3`. |
| `metassembler` | `beaver` | The meta-assembler. Options: `beaver`, `transmeta`. |
| `transmeta_use_annotation` | `true` | When TransMeta is the meta-assembler, guide it with the reference annotation. |
| `aletsch_fallback_insert_size` | – | If Aletsch cannot estimate the insert-size mean (zero profile), use this value instead. Handy for weird libraries; `null` means "trust Aletsch". |
| `aletsch_fallback_insert_std` | `10` | Standard deviation paired with the fallback insert-size mean. |
| `aletsch_keep_bam` | `false` | Keep the Aletsch output BAM (extra disk). |
| `beaver_keep_aletsch_gtf` | `false` | Keep the per-sample Aletsch GTFs after Beaver runs (extra disk). |
| `stringtie3_enable_nascent_assembly` | `false` | StringTie 3: run *de novo* assembly of nascent (unspliced) RNA (`-N`). Mutually exclusive with the next option. |
| `stringtie3_include_nascent_rna` | `false` | StringTie 3: include nascent RNA in the final assembly (`--nasc`). |

### 4.10 Polishing

Polishing is where the pipeline decides, per transcript, whether it is
real. It extracts every intron, predicts the spliceosome machinery
(intronIC), and classifies introns using SpliceAI + MaxEnt scores and
coverage. Transcripts carrying retained introns, read-throughs, weak
retained introns, or other artifacts are stripped from the final set and
published separately so nothing is silently lost.

| Parameter | Default | What it does |
|-----------|---------|--------------|
| `xloci_upstream_flank` | `100` | Bases of exonic context to keep upstream of each intron when extracting intron sequences. |
| `xloci_downstream_flank` | `100` | Bases of exonic context to keep downstream of each intron. |
| `isotools_classify_spliceai_min_ss_signal` | `0.02` | Minimum SpliceAI signal for a splice site to count as "significant". |
| `isotools_classify_rt_frequency_threshold` | `0.5` | Frequency threshold for classifying a transcript as read-through. |
| `isotools_classify_intron_frequency_threshold` | `0.5` | Frequency threshold for classifying intron retention. |
| `isotools_classify_maxent_min_ss_signal` | `1.5` | Minimum MaxEnt score for a splice site to count as real. |
| `do_twopass_polish` | `false` | Re-review the transcripts that were discarded for **intron retention** in the first pass, this time ignoring UTR retentions (which are often real). The second-pass `--toga` reference is the union of the reference annotation and the XORF ORF models. Rescued transcripts get their own XORF round; those HQ ORFs are merged with the first-pass HQ, truncation discards from both XORF runs are unioned (`10_final/truncations`), then NMD runs on the combined HQ. Costs a second classification round plus a second XORF. |

### 4.10a ANNEVO annotation (optional)

By default the reference annotation (`--annotation`) is the only gene model
xasm compares against during polishing. With ANNEVO you can add a second,
*ab initio* gene model: POLISH merges it with the reference annotation into one
coordinate-sorted BED and feeds that single file to `iso-orphan` (`--ref`),
`iso-fusion` (`--ref`), `iso-classify` (`--toga`) and the two-pass `--toga`
union. Introns supported by either source stay categorized.

Two ways to supply it:

- `--annevo_annotation FILE` — an ANNEVO GFF3 (or any GTF/GFF/BED) you already
  have. Converted to BED12 with `gxf2bed`.
- `--annevo_predict` — run ANNEVO on the genome and use its output. Requires
  `--annevo_lineage`. Ignored when `--annevo_annotation` is given. ANNEVO needs
  a container profile (`docker`/`apptainer`/`singularity`); it is not available
  with `-profile conda`.

| Parameter | Default | What it does |
|-----------|---------|--------------|
| `annevo_annotation` | `null` | Precomputed ANNEVO annotation (`.bed`/`.gtf`/`.gff`/`.gff3`, gzipped input is not supported). Wins over `annevo_predict`. |
| `annevo_predict` | `false` | Run ANNEVO on the genome and use the prediction as second annotation. |
| `annevo_lineage` | `null` | ANNEVO lineage: `Mammalia`, `Insecta`, `Aves`, `Actinopteri`, `Magnoliopsida`, `Fungi`. **Required** with `annevo_predict`. |
| `annevo_scatter` | `chromosome` | How the genome is split: `none`, `chromosome` (one task per sequence), `weighted` (length-balanced bins). |
| `annevo_bins` | `8` | Number of bins for `--annevo_scatter weighted`. |
| `annevo_overlap_pred` | `null` | Allow overlapping gene predictions. Defaults to lineage-specific ANNEVO behaviour (`Mammalia`, `Actinopteri`: on). |

> ANNEVO output is published under `08_annevo`: the `*.annevo.gff3` (only with
> `--annevo_predict`) plus the BED12 POLISH merges into the reference (passthrough
> for `.bed` input). The merged reference itself lives in the work directory
> (`<prefix>.reference.sorted.bed`).

### 4.10b TIBERIUS annotation (optional)

TIBERIUS works exactly like ANNEVO above — a further *ab initio* gene model that
POLISH merges with the reference annotation into the single coordinate-sorted
BED fed to `iso-orphan`, `iso-fusion` and `iso-classify` — with one difference:
the two TIBERIUS sources may be **combined**. ANNEVO and TIBERIUS may also be
active at the same time; every active source contributes its BED to the merge.

Two ways to supply it (both at once is allowed):

- `--tiberius_annotation FILE` — a TIBERIUS GTF/GFF3 (or any GTF/GFF/BED) you
  already have. Converted to BED12 with `gxf2bed`.
- `--tiberius_predict` — run TIBERIUS on the genome and use its output. Requires
  `--tiberius_model_cfg` (a Hiller alias from `/opt/tiberius/models.tsv`, e.g.
  `mammalia_nosoftmasking_v2`). TIBERIUS needs a container profile
  (`docker`/`apptainer`/`singularity`); like ANNEVO it is not available with
  `-profile conda`.

| Parameter | Default | What it does |
|-----------|---------|--------------|
| `tiberius_annotation` | `null` | Precomputed TIBERIUS annotation (`.bed`/`.gtf`/`.gff`/`.gff3`, gzipped input is not supported). May be combined with `tiberius_predict`. |
| `tiberius_predict` | `false` | Run TIBERIUS on the genome and use the prediction as third annotation. |
| `tiberius_model_cfg` | `null` | TIBERIUS model alias from `/opt/tiberius/models.tsv`. **Required** with `tiberius_predict`. |
| `tiberius_scatter` | `chromosome` | How the genome is split: `none`, `chromosome` (one task per sequence), `weighted` (length-balanced bins). |
| `tiberius_bins` | `8` | Number of bins for `--tiberius_scatter weighted`. |
| `tiberius_batch_size` | `null` | Batch size forwarded as TIBERIUS `--batch_size`. |
| `tiberius_seq_len` | `null` | Sequence length forwarded as TIBERIUS `--seq_len`. |
| `tiberius_protseq` | `false` | Also emit the TIBERIUS protein FASTA (`*.prot`). |
| `tiberius_codingseq` | `false` | Also emit the TIBERIUS coding-sequence FASTA (`*.cds`). |
| `tiberius_extra_args` | `''` | Extra args passed through to `tiberius.py`. |

> TIBERIUS output is published under `08_tiberius`: the merged `*.gtf`/`*.gff3`
> (plus `*.prot`/`*.cds` when requested) and the BED12(s) POLISH merges into the
> reference (passthrough for `.bed` input).

### 4.11 BigBed conversion

The final transcripts are published as BED files. If you want UCSC
genome-browser tracks, the pipeline also converts every final BED into a
`.bb` (BigBed) file. BigBed needs an `autosql` definition when the BED has
more than the standard fields.

| Parameter | Default | What it does |
|-----------|---------|--------------|
| `autosql` | – | Path to an `.as` (AutoSQL) file describing extra BED columns. Only needed if you intend to load the tracks into the UCSC browser. |
| `skip_bb_conversion` | `false` | Skip BigBed conversion entirely. The pipeline then reports `10_final/*.bed` as its final output instead of `11_bbs/*.bb`. |

### 4.12 Email notifications (advanced)

On completion, xasm can email you a per-sample summary report (trimmed read
counts, decontamination and alignment rates, assembled transcript counts,
plus the final predictions). Requires a working SMTP account — or your
cluster's `mailx` command.

| Parameter | Default | What it does |
|-----------|---------|--------------|
| `sent_email` | `false` | Master switch: send the completion email or not. |
| `email_to` | – | Recipient address(es). |
| `email_on_fail` | – | Extra address(es) to notify if the pipeline fails. |
| `plaintext_email` | `false` | Send a plain-text email instead of an HTML report. |
| `use_mailx` | `false` | Use the system `mailx` command instead of SMTP (all `smtp_*` options are then ignored). |
| `interactive` | `false` | Also write a browsable interactive HTML report. |
| `smtp_server` | `smtp.gmail.com` | SMTP server hostname. |
| `smtp_port` | `465` | SMTP port. |
| `smtp_user` | – | SMTP username (sender). |
| `smtp_password` | – | SMTP password (sender). Keep it out of git! |
| `smtp_security` | `ssl` | Connection security: `none`, `tls`, or `ssl`. |

### 4.13 Boilerplate (rarely touched)

| Parameter | Default | What it does |
|-----------|---------|--------------|
| `publish_dir_mode` | `copy` | How results are published (`copy`, `symlink`, `link`). |
| `publish_all` | `false` | Publish everything a process produces instead of only the main outputs. |
| `help` | `false` | Show the command-line help and exit. |
| `help_full` | `false` | Show the full help including hidden parameters. |
| `show_hidden` | `false` | Show hidden parameters in help output. |
| `version` | `false` | Print the pipeline version and exit. |
| `trace_report_suffix` | timestamp | Suffix for the pipeline report files in `pipeline_info/`. |
| `config_profile_name` | – | Display name of the profile (cosmetic). |
| `config_profile_description` | – | Description of the profile (cosmetic). |

---

## 5. Profiles

Profiles are chosen with `-profile <name>` (comma-separate to combine, e.g.
`-profile apptainer,slurm`).

| Profile | What it does |
|---------|--------------|
| `local` | Run everything on the local machine (default). |
| `slurm` | Submit jobs to a SLURM cluster (queue `batch`/`long`, job arrays for the heavy steps). |
| `apptainer` | Use Apptainer containers (recommended for clusters). |
| `singularity` | Use Singularity containers. |
| `docker` | Use Docker containers (recommended for laptops/workstations). |
| `podman`, `shifter`, `charliecloud` | Other container engines. |
| `conda` / `mamba` | Use Conda environments instead of containers. |
| `wave` | Build environments on the fly via Wave. |
| `gpu` | Add GPU flags to the container run (only useful if a process requests an accelerator). |
| `gitpod` | Limited resources (4 CPUs, 8 GB) for Gitpod-style environments. |
| `debug` | More verbose Nextflow diagnostics. |
| `test`, `test-sb`, `test-tm`, `test-rustar`, `test-bqc-star`, `test-polish` | Smoke-test profiles that run the bundled test data (see section 2); `test-rustar` is a smoke-checked CBQ+ruSTAR path, `test-bqc-star` is the CBQ+`bqc`+decode+STAR path checked against the same samplesheet numbers as `test` (see [4.6a](#46a-cbq-reads-bqtools--bqc)). |

---

## 6. Running on a SLURM cluster

The helper script `assets/scripts/xasm.sh` runs one pipeline per line of a
manifest file, each as its own SLURM job array task. Edit the four path
variables at the top of the script (container cache, working dir, pipeline
dir, and the manifest path), then:

```bash
sbatch --array=1-<N> xasm.sh
```

where `<N>` is the number of lines in your manifest. Manifest format is
tab-separated, no header, absolute paths:

```
<input_dir>  <genome>  <annotation>  <splice_scores_dir>  <repeats>  <prefix>
```

Each task writes its own `params.json` and `results/` under a
`<prefix>_xasm/` subdirectory, and Nextflow itself submits all compute jobs
to SLURM as child jobs. The SLURM-specific tuning (queues, array sizes,
resource tiers) lives in `nextflow.config` under the `slurm` profile.

---

## 7. Checkpoints: resuming from the middle

The full pipeline is the expensive part (QC → alignment → assembly). If you
already ran it once, you can re-run only the tail steps:

**`--from polish`** — start at the polishing step, skipping everything up to
and including metassembly.

```
nextflow run main.nf -params-file params.json -profile apptainer \
    --from polish --polish_path /path/to/metassembly.gtf
```

Required: `polish_path` (your metassembly GTF/BED), `genome`, `annotation`.
`repeats` and `splice_scores_dir` are strongly recommended.

**`--from bigbed`** — skip everything and just convert existing final BED
files to BigBed.

```
nextflow run main.nf -params-file params.json -profile apptainer \
    --from bigbed --all_bed_path /path/to/beds
```

Required: `all_bed_path` (directory of final `.bed` files), `genome`.
Useful when you ran with `skip_bb_conversion = true` before and changed your
mind.

---

## 8. Outputs

Everything lands under `output_dir` (default `./results`):

```
results/
├── 00_prepare/
│   ├── deacon_index/          built or downloaded Deacon index
│   ├── strand/                strand inference support: reference CDS FASTA,
│   │                          reusable Salmon index (bqc_salmon_index/) and
│   │                          per-sample bqc sniff strand reports
│   └── wget/                  downloaded index files
├── 01_deacon_filter/          decontaminated reads (*fastq, symlinks)
├── 02_star_index/             STAR genome index
├── 03_star_2pass/             final alignments (*bam, *bai)
├── 04_junctions/              merged, filtered junction set (*bed)
├── 05_coverage/               merged median coverage (*bigwig)
├── 06_beaver/  (or 06_transmeta/)
│                              per-chromosome assemblies + the final
│                              concatenated metassembly (<prefix>.gtf)
├── 07_remove_dirt/            cleaned metassembly (*gtf)
├── 08_gxf2bed/                BED versions of the annotation + assembly (*bed)
├── 08_annevo/                 ANNEVO annotation (*annevo.gff3 + *.bed),
│                              only when an ANNEVO source is given
├── 08_tiberius/               TIBERIUS annotation (*.gtf/*.gff3 + *.bed,
│                              plus *.prot/*.cds when requested),
│                              only when a TIBERIUS source is given
├── 09_polish/
│   ├── fusions/               fusion candidate transcripts (*bed)
│   ├── orphans/               iso-orphan classification (*bed, *tsv)
│   ├── intron_sequences/      extracted intron sequences (*tsv)
│   ├── iic/                   intronIC spliceosome predictions (*tsv)
│   ├── classify/              intron classifications (*tsv)
│   ├── retentions/            retention descriptors
│   ├── xorf/                  ORF calls and truncation/duplicate discards (*bed)
│   └── twopass/               second-pass re-classification (only with
│                              --do_twopass_polish)
├── 10_final/                  ← FINAL transcripts
│   ├── <prefix>.bed           the polished, clean transcriptome
│   ├── retentions/            transcripts discarded as retained introns
│   ├── strong_rts/            read-through discards (strong signal)
│   ├── weak_rts/              read-through discards (weak signal)
│   ├── artifacts/             other artifact discards
│   ├── truncations/           XORF 3′UTR truncation discards (only when present)
│   ├── nmd/                   NMD-positive transcripts (only when present)
│   ├── scraps/                orphan/scrap transcripts
│   └── fusions/               fusion transcripts
├── 11_bbs/                    ← FINAL tracks (BigBed, if not skipped)
│   └── <prefix>.bb            browser-ready versions of all 10_final BEDs
├── samplesheets/
│   └── samplesheet.csv        one row per sample with QC stats
└── pipeline_info/
    └── execution_timeline/report/trace + pipeline DAG (html/txt)
```

**Which files are "the answer"?** The final transcriptome is
`10_final/<prefix>.bed` (plus the per-category discard BEDs if you want to
know what was filtered out), or the BigBed tracks in `11_bbs/` if you kept
conversion enabled. On completion the pipeline prints the paths of these
files to the console.

### The samplesheet

`samplesheets/samplesheet.csv` holds one line per sample with the QC story
of that sample, in order:

| Column | Meaning |
|--------|---------|
| 1 | Sample id (fastq prefix) |
| 2–3 | Input fastq file names (mate 2 empty for single-end) |
| 4 | Reads left after trimming |
| 5 | Percent of reads kept after trimming |
| 6 | Percent of sequences retained after Deacon decontamination |
| 7 | Percent of reads uniquely mapped by STAR |
| 8 | Final BAM size (in 100 MB units) |
| 9 | Number of transcripts assembled for that sample |
| 10 | Inferred strandedness (`forward`/`reverse`/`unstranded`; see [4.6a](#46a-cbq-reads-bqtools--bqc)) |

---

## 9. Tuning compute and memory

Per-step CPU/memory/time and the retry strategy are defined in
`nextflow.config` (`process` section, labels like `process_low` /
`process_medium` / `process_high`). The sensible defaults scale by retry
attempt; only touch them if a specific step on your data keeps failing or is
over-provisioned. If the same step always dies with an out-of-memory error,
find its label in the config and raise the corresponding tier.
