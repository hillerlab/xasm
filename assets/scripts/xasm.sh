#!/bin/bash
# ─────────────────────────────────────────────────────────────────────────────
# SLURM submission script for xasm
#
# Each array task runs one soft-masking job. Nextflow itself is the
# "main job" — it submits all compute work as child SLURM jobs and only needs
# a small memory footprint.
#
# MANIFEST FILE FORMAT (species_list)
# ─────────────────────────────────────
# A tab-separated file with one run per line, no header:
#
#   <input_dir>  <genome> <annotation> <splice_scores_dir> <repeats> <prefix>
#
# Example:
#   /path/to/input_dir  hg38.2bit  hg38.bed  splice_scores_dir  hg38
#
# Paths must be absolute. 
#
# USAGE
# ─────
# Edit the four path variables below, then submit with:
#   sbatch --array=1-<N> xasm.sh
# where <N> is the number of lines in your manifest file.
# ─────────────────────────────────────────────────────────────────────────────

#SBATCH --job-name=XASM
#SBATCH --array=1-10        # set upper bound to number of lines in species_list
#SBATCH -t 2-0
#SBATCH --output=/path/to/logs/%A.%a.out  # MODIFY THIS!
#SBATCH --error=/path/to/logs/%A.%a.err   # MODIFY THIS!
#SBATCH --mem=20G          # memory for the Nextflow process itself (not compute jobs)
#SBATCH -p long            # partition name
#SBATCH -q long            # queue name

# ── Load required modules (adjust to your cluster's module system) ────────────
module load nextflow
module load openjdk

# ── Environment ───────────────────────────────────────────────────────────────
export SLURM_SKIP_EPILOG=1

# Directory where Apptainer caches pulled container images
export NXF_APPTAINER_CACHEDIR=/scratch/$USER/xasm/apptainer

# Optional: pre-build a named SIF to avoid the auto-derived cache filename.
# Build once with:
#   apptainer build $NXF_APPTAINER_CACHEDIR/xasm.sif \
#       ghcr.io/hillerlab/xasm:latest
# Then uncomment:
# export NXF_CONTAINER_IMAGE=$NXF_APPTAINER_CACHEDIR/xasm.sif

# Give Nextflow's JVM enough heap for large runs (thousands of jobs)
export NXF_OPTS="-Xms4g -Xmx16g"

# ── Paths — edit these ────────────────────────────────────────────────────────
species_list="/path/to/manifest.tsv"   # tab-separated manifest (see format above)
working_dir="/path/to/output"          # one subdirectory per genome will be created here
pipeline_dir="/path/to/xasm"       # cloned pipeline repo

# ── Parse manifest line for this array task ───────────────────────────────────
row=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$species_list")
input_dir=$(echo "$row" | cut -f1)
genome=$(echo "$row"   | cut -f2)
annotation=$(echo "$row" | cut -f3)
splice_scores_dir=$(echo "$row" | cut -f4)
repeats=$(echo "$row" | cut -f5)
prefix=$(echo "$row" | cut -f6)


if [[ -z "$input_dir" || -z "$genome" || -z "$annotation" || -z "$splice_scores_dir" || -z "$repeats" || -z "$prefix" ]]; then
    echo "ERROR: could not parse line ${SLURM_ARRAY_TASK_ID} of ${species_list}" >&2
    exit 1
fi

# ── Per-pair working directory ─────────────────────────────────────────────────
run_dir="${working_dir}/${prefix}_xasm"
mkdir -p "${run_dir}/logs"

# ── Write params.json for this pair ───────────────────────────────────────────
# Scientific parameters go here; infrastructure stays in nextflow.config.
cat > "${run_dir}/params.json" <<EOF
{
    "//1": "── Input / Output + Mandatory inputs ──────────────────────────────────────────────────────",
    "input_dir":           ${input_dir},
    "output_dir":          ${run_dir}/results,
    "genome":              ${genome},
    "annotation":          ${annotation},
    "splice_scores_dir":   ${splice_scores_dir},
    "repeats":             ${repeats},

    "//1a": "── Checkpoints (--from alignment|polish|bigbed) ─────────────────────────────────",
    "bam_dir":             null,
    "polish_path":         null,

    "//1b": "── Checkpoint: start from bigbed step (skip metassembly + polishing) ─────────────────────",
    "all_bed_path":        null,

    "//2": "── Linting [optional] ──────────────────────────────────",
    "fq_skip_linting_at_start":          false,
    "fq_skip_linting_after_trimming":    true,
    "fq_skip_linting_after_deacon":      true,

    "//3": "── Decontamination [optional] ───────────────────────────────────────────",
    "deacon_index_path":                   null,
    "deacon_download_index":               false,
    "deacon_make_single_index":            true,
    "deacon_single_index_use_background":  true,

    "//4": "── Alignment [optional] ──────────────────────────────────",
    "star_index_path":              null,
    "star_ignore_gtf_for_index":    false,
    "star_ignore_gtf_for_mapping":  true,

    "//5": "── Assembly [optional] ──────────────────────────────────",
    "skip_assembly":  false.

    "//6: ── Bigbed conversion [optional] ───────────────────────────────────────────",
    "autosql": null,
    "skip_bb_conversion": false
}
EOF

# cd into run_dir so each run's .nextflow.log is saved there
cd "$run_dir"

nextflow run "${pipeline_dir}/main.nf" \
    -params-file "${run_dir}/params.json" \
    -profile     apptainer,slurm \
    -w           "${run_dir}/work"
