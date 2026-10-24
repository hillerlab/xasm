<p align="center">
  <p align="center">
    <img width=100 align="center" src="./assets/figures/logo.png" >
  </p>

<p align="center">
  <picture>
    <source
      media="(prefers-color-scheme: dark)"
      srcset="./assets/figures/hillerlab-dark.png"
    >
    <source
      media="(prefers-color-scheme: light)"
      srcset="./assets/figures/hillerlab-light.png"
    >
    <img
      width="200"
      alt="Hiller Lab"
      src="./assets/figures/hillerlab-light.png"
    >
  </picture>
</p>

  <span>
    <h1 align="center">
        xasm
    </h1>
  </span>

  <p align="center">
    <a href="https://github.com/hillerlab/xasm" reference="_blank">
      <img alt="GitHub License" src="https://img.shields.io/github/license/hillerlab/xasm?color=blue">
    </a>
  </p>

  <p align="center">
    <samp>
        <span> meta-assemble transcriptomes from bulk RNA-seq datasets at scale  </span>
        <br>
        <span> The Hiller Lab at the Senckenberg Research Institute </span>
        <br>
        <br>
        <a href="https://github.com/alejandrogzi/xasm/blob/main/assets/docs/usage.md">usage</a> .
        <a href="https://nbisweden.github.io/workshop-RNAseq/2011/lab_assembly.html">metassembly</a> .
        <a href="https://github.com/hillerlab/xasm/blob/main/assets/pipeline/xasm.mermaid">pipeline</a> .
        <a href="https://hillerlab.com/">us</a> 
    </samp>
  </p>

</p>

---

<div align="center">

<pre style="font-size: 18px;">

[ QC + trim ] ──► [ decontaminate ] ──► [ STAR-2pass align ]

S1  ───────●●──────────●●●──────────────────●●────────────
S2  ──────────────●●●────────●●─────────────●●●───────────
S3  ────●●──────────────────●●●●────────────────●─────────
S4  ─────────●●●────────────●●───────────●●●──────────────
S5  ─────────────●●────────────────●●●──────────●●────────

▼

═══════════════●●════════════●●●●════════════●●═══════════

└────────────────── metassembled transcriptome ──────────────┘
</pre>

</div>

---

## Usage

> [!NOTE]
> Requirements: Nextflow ≥ 25.04.6, Docker or Apptainer, Java.

```bash
git clone https://github.com/hillerlab/xasm.git
cd xasm
```

Edit `params.json` (set `genome`, `assembly_prefix`), then:

```bash
# Docker
nextflow run main.nf -params-file params.json -profile docker

# Apptainer / Singularity
nextflow run main.nf -params-file params.json -profile apptainer
```

End-to-end tests use the committed 256 KB fixture and validate stable counts
through trimming, decontamination, alignment, chromosome-split assembly,
metassembly, and polishing:

```bash
# Default Aletsch + Beaver smoke test
assets/ci/e2e/run.sh test

# StringTie 3 + Beaver, TransMeta, ruSTAR + CBQ, bqc + STAR, or the complete matrix
assets/ci/e2e/run.sh test-sb
assets/ci/e2e/run.sh test-tm
assets/ci/e2e/run.sh test-rustar
assets/ci/e2e/run.sh test-bqc-star
assets/ci/e2e/run.sh all
```

Set `TEST_ENGINE=apptainer` to change container engines. The profiles write to
separate directories under `test_results/`; `test-polish` remains the small
polishing-only fixture. Activate the compatible environment first with
`mamba activate nextflow`, or set `NEXTFLOW_BIN` explicitly.

Retention-only review can be enabled for normal and polishing-checkpoint runs
with `--do_twopass_polish true`. It reuses the first-pass intronIC evidence,
reclassifies retention discards, ignores UTR retentions in the second pass, and
adds rescued transcripts back to the clean HQ set.

> [!NOTE]
> You can also specify these options directly in `params.json`.

A helper sh script is provided to run the pipeline on a SLURM cluster. See details below.

<details>
<summary>Click to expand</summary>


Edit the path variables at the top of `assets/hpc/do_xasm.sh` (cache dir, container image, manifest path), then submit:

```bash
sbatch --array=1-<N> do_xasm.sh
```

Each array task spawns one Nextflow head job that submits all compute as child SLURM jobs.

STAR_ALIGN_1PASS, STAR_ALIGN_2PASS, FASTP, DEACON_FILTER, and BEDGRAPHTOBIGWIG run as SLURM job arrays. Partition routing, array sizes, and resource tiers are documented inline in `nextflow.config` — edit there to match your cluster.

</details>

---

## Output

```
results/
├── 00_prepare/
│   ├── deacon_index/          built or downloaded Deacon index
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
├── 11_bbs/                    ← FINAL tracks (BigBed, if not skipped), one per 10_final class
│   ├── <prefix>.flnc.bb       the polished, clean transcriptome
│   ├── <prefix>.retentions.bb
│   ├── <prefix>.strong_rts.bb
│   ├── <prefix>.weak_rts.bb
│   ├── <prefix>.artifacts.bb
│   ├── <prefix>.truncations.bb (only when present)
│   ├── <prefix>.nmd.bb        (only when present)
│   ├── <prefix>.scraps.bb
│   └── <prefix>.fusions.bb
├── samplesheets/
│   └── samplesheet.csv        one row per sample with QC stats
└── pipeline_info/
    └── execution_timeline/report/trace + pipeline DAG (html/txt)
```

---

## Where to edit

| File | What |
|------|------|
| `params.json` | Genome paths, alignment settings, checkpoints — per run |
| `nextflow.config` | Compute resources, profiles, container, SLURM — rarely |
