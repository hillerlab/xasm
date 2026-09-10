#!/usr/bin/env nextflow

/*
Copyright (c) 2026 The Hiller Lab at the Senckenberg Gessellschaft für Naturforschung
Distributed under the terms of the Apache License, Version 2.0.
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    xasm

    meta-assemble bulk RNA-seq datasets at scale
    Authors: Alejandro Gonzales-Irribarren, Michael Hiller

    GitHub:  https://github.com/hillerlab/xasm
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    HELP
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

if (params.help) {
    log.info """
    xasm v${workflow.manifest.version}
    Meta-assemble bulk RNA-seq datasets at scale

    Authors: ${workflow.manifest.author}
    Github:  ${workflow.manifest.homePage}

    Usage (full run):
        nextflow run main.nf \\
            --input_dir   PATH    Path to input directory (fastqs)
            --genome      PATH    Path to genome file (fasta/2bit)
            --annotation  PATH    Path to annotation file (gtf/gff/bed) [defines STAR + iso-orphan inputs]
            --splice_scores_dir  PATH      Path to splicing scores dir
            --repeats     PATH      Path to repeats.{bed/gff/gtf}
            --output_dir  PATH    Output directory [default: ./results]
            -profile         apptainer,slurm

    Pass all parameters from a JSON file (replaces old --params_from_file):
        nextflow run main.nf -params-file my_params.json

    Required parameters (full run) :
        --input_dir           PATH    Path to input directory (fastqs)
        --genome              PATH    Path to genome file (fasta/2bit)
        --annotation          PATH    Path to annotation file (gtf/gff/bed) [defines STAR + iso-orphan inputs]
        --output_dir          PATH    Output directory [default: ./results]

    Optional parameters (common):
        --splice_scores_dir   PATH      Path to splicing scores dir 
        --repeats             PATH      Path to repeats.{bed/gff/gtf}
        --from                STRING    Checkpoint to resume from [options: polish]
        --polish_path         PATH      Path to polished assembly [required if --from polish]
        --do_twopass_polish   BOOLEAN   Re-review retention discards with --ignore-utr [default: false]
        --xorf_call_orfs      BOOLEAN   Run XORF ORF calling on first-pass HQ transcripts [default: false]
        --xorf_esm            BOOLEAN   Add CPU ESMFold2-Fast pLDDT scores to BLAST [default: false]
        --xorf_esmfold_local_weights PATH Local Hugging Face hub cache (ESMFold2-Fast + ESMC-6B); skips download [default: null]
        --xorf_custom_database PATH     Custom protein database for ORF calling (.dmnd/.dmnd.gz replaces the default; .fa/.fasta appended to SwissProt) [default: null; rest of XORF options in params.json]

    Optional parameters (ANNEVO annotation):
        --annevo_annotation   PATH      Precomputed ANNEVO annotation (.bed/.gtf/.gff/.gff3); merged with
                                        the reference annotation and used by iso-orphan, iso-fusion and
                                        iso-classify [default: null]
        --annevo_predict      BOOLEAN   Run ANNEVO on the genome and use its output as second annotation.
                                        Ignored when --annevo_annotation is given [default: false]
        --annevo_lineage      STRING    ANNEVO lineage [options: Mammalia, Insecta, Aves, Actinopteri,
                                        Magnoliopsida, Fungi] [required with --annevo_predict]
        --annevo_scatter      STRING    ANNEVO scatter mode [options: none, chromosome, weighted]
                                        [default: chromosome]
        --annevo_bins         INT       Number of bins for --annevo_scatter weighted [default: 8]
        --annevo_overlap_pred BOOLEAN   Allow overlapping gene predictions [default: lineage-specific]

    Optional parameters (TIBERIUS annotation):
        --tiberius_annotation   PATH      Precomputed TIBERIUS annotation (.bed/.gtf/.gff/.gff3); merged with
                                          the reference annotation and used by iso-orphan, iso-fusion and
                                          iso-classify [default: null]
        --tiberius_predict      BOOLEAN   Run TIBERIUS on the genome and use its output as third annotation.
                                          May be combined with --tiberius_annotation [default: false]
        --tiberius_model_cfg    STRING    TIBERIUS model alias from /opt/tiberius/models.tsv
                                          [required with --tiberius_predict]
        --tiberius_scatter      STRING    TIBERIUS scatter mode [options: none, chromosome, weighted]
                                          [default: chromosome]
        --tiberius_bins         INT       Number of bins for --tiberius_scatter weighted [default: 8]
        --tiberius_batch_size   INT       TIBERIUS --batch_size [default: null]
        --tiberius_seq_len      INT       TIBERIUS --seq_len [default: null]
        --tiberius_protseq      BOOLEAN   Also emit TIBERIUS protein FASTA [default: false]
        --tiberius_codingseq    BOOLEAN   Also emit TIBERIUS coding-sequence FASTA [default: false]
        --tiberius_extra_args   STRING    Extra args passed through to tiberius.py [default: '']

    Optional parameters (CBQ reads):
        --aligner                          STRING    Aligner [options: STAR, ruSTAR] [default: STAR]
                                                     STAR reads fastq only (CBQ is decoded after deacon);
                                                     ruSTAR reads fastq and cbq
        --bqtools_encode_fastqs            BOOLEAN   Encode fastq inputs to .cbq up front; QC then runs
                                                     via bqc instead of fastp. With --aligner STAR the
                                                     path is encode → bqc → deacon → decode → STAR
                                                     [default: false]
        --bqtools_encode_before_alignment  BOOLEAN   Keep fastp+deacon on fastq and encode to .cbq only
                                                     for the aligner (ruSTAR only) [default: false]
        Native .cbq files in --input_dir are detected automatically and need neither flag.

        --infer_strandedness               BOOLEAN   Infer each sample's strandedness from its raw CBQ
                                                     before trimming (bqc sniff strand) [default: true]
        --strand_transcriptome             PATH      Prebuilt transcriptome FASTA for the strand index;
                                                     skips XLOCI CDS extraction [default: null]
        --strand_salmon_index              PATH      Prebuilt Salmon 2.x index directory; skips CDS
                                                     extraction + index build [default: null]

    Profiles:
        local       Run on local machine (default)
        slurm       Submit jobs to SLURM cluster
        conda       Use conda environments
        apptainer   Use Apptainer containers
        singularity Use Singularity containers
        docker      Use Docker containers
        arm         Force linux/amd64 containers on ARM hosts (use as: -profile docker,arm)
        test        Run with bundled test data

    Use --help to show this message.
    """.stripIndent()
    System.exit(0)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { METASSEMBLE } from './subworkflows/metassembly/main'
include { POLISH } from './subworkflows/polish/main'
include { ANNEVO_ANNOTATION } from './subworkflows/annevo/main'
include { TIBERIUS_ANNOTATION } from './subworkflows/tiberius/main'
include { EMAIL_RESULTS } from './modules/custom/email/main'
include { TWOBIT_TO_FA } from './modules/custom/ucsc/twobittofa/main'
include { CHROMSIZE } from './modules/custom/chromsize/main'
include { BED2GTF } from './modules/custom/bed2gtf/main'
include { GXF2BED } from './modules/custom/gxf2bed/main'
include { GUNZIP as GUNZIP_FASTA } from './modules/custom/gunzip/main'
include { GUNZIP as GUNZIP_GTF } from './modules/custom/gunzip/main'
include { SORT_BED } from './modules/custom/sort/main'

include { BEDTOBIGBED } from './modules/custom/bigtools/bedtobigbed/main'
include { BEDTOBIGBED as BEDTOBIGBED_HQ } from './modules/custom/bigtools/bedtobigbed/main'
include { BEDTOBIGBED as BEDTOBIGBED_RETENTION } from './modules/custom/bigtools/bedtobigbed/main'
include { BEDTOBIGBED as BEDTOBIGBED_TRUNCATIONS } from './modules/custom/bigtools/bedtobigbed/main'
include { BEDTOBIGBED as BEDTOBIGBED_STRONG_RTS } from './modules/custom/bigtools/bedtobigbed/main'
include { BEDTOBIGBED as BEDTOBIGBED_WEAK_RTS } from './modules/custom/bigtools/bedtobigbed/main'
include { BEDTOBIGBED as BEDTOBIGBED_ARTIFACTS } from './modules/custom/bigtools/bedtobigbed/main'
include { BEDTOBIGBED as BEDTOBIGBED_FUSIONS } from './modules/custom/bigtools/bedtobigbed/main'
include { BEDTOBIGBED as BEDTOBIGBED_SCRAPS } from './modules/custom/bigtools/bedtobigbed/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// ANNEVO annotation options. Shared by every entry point that reaches POLISH.
//
def validateAnnevo() {
    def errors = []
    def lineages = ['Mammalia', 'Insecta', 'Aves', 'Actinopteri', 'Magnoliopsida', 'Fungi']

    if (params.annevo_predict && !params.annevo_lineage) {
        errors << "  --annevo_predict requires --annevo_lineage [options: ${lineages.join(', ')}]"
    }

    if (params.annevo_annotation) {
        def source = file(params.annevo_annotation)
        if (!source.exists()) {
            errors << "  --annevo_annotation '${params.annevo_annotation}' does not exist"
        } else if (!(source.name ==~ /.*\.(bed|gtf|gff|gff3)$/)) {
            errors << "  --annevo_annotation must be .bed/.gtf/.gff/.gff3 (gzip is not supported)"
        }
    }

    return errors
}

//
// TIBERIUS annotation options. Shared by every entry point that reaches POLISH.
// Unlike ANNEVO, --tiberius_annotation and --tiberius_predict may be combined.
//
def validateTiberius() {
    def errors = []
    def scatter_modes = ['none', 'chromosome', 'weighted']

    if (params.tiberius_predict && !params.tiberius_model_cfg) {
        errors << "  --tiberius_predict requires --tiberius_model_cfg (hiller alias from /opt/tiberius/models.tsv)"
    }

    if (params.tiberius_annotation) {
        def source = file(params.tiberius_annotation)
        if (!source.exists()) {
            errors << "  --tiberius_annotation '${params.tiberius_annotation}' does not exist"
        } else if (!(source.name ==~ /.*\.(bed|gtf|gff|gff3)$/)) {
            errors << "  --tiberius_annotation must be .bed/.gtf/.gff/.gff3 (gzip is not supported)"
        }
    }

    if (!(params.tiberius_scatter in scatter_modes)) {
        errors << "  --tiberius_scatter must be one of: ${scatter_modes.join(', ')} (got '${params.tiberius_scatter}')"
    }

    if (params.tiberius_scatter == 'weighted' && params.tiberius_bins != null && (params.tiberius_bins as int) < 1) {
        errors << "  --tiberius_scatter weighted requires --tiberius_bins >= 1"
    }

    return errors
}

//
// XORF ESM options. Mirrors the standalone xorf validation for local weights.
//
def validateXorf() {
    def errors = []

    if (params.xorf_esm && params.xorf_esmfold_local_weights) {
        def weights = file(params.xorf_esmfold_local_weights)
        if (!weights.exists()) {
            errors << "  --xorf_esmfold_local_weights does not exist: '${params.xorf_esmfold_local_weights}'"
        } else if (!file("${weights}/models--biohub--ESMFold2-Fast").exists() || !file("${weights}/models--biohub--ESMC-6B").exists()) {
            errors << "  --xorf_esmfold_local_weights is not a Hugging Face hub cache (need models--biohub--ESMFold2-Fast and models--biohub--ESMC-6B)"
        }
    }

    return errors
}

def validateRun() {
    def errors = []
    if (!params.input_dir)   errors << "  --input_dir is required"
    if (!params.genome)      errors << "  --genome is required"
    if (!params.annotation)  errors << "  --annotation is required"

    if (!(params.aligner in ['STAR', 'ruSTAR'])) {
        errors << "  --aligner must be 'STAR' or 'ruSTAR' (got '${params.aligner}')"
    }

    if (params.bqtools_encode_fastqs && params.bqtools_encode_before_alignment) {
        errors << "  --bqtools_encode_fastqs and --bqtools_encode_before_alignment are mutually exclusive"
    }

    if (params.do_twopass_polish && !params.xorf_call_orfs) {
        errors << "  --do_twopass_polish requires --xorf_call_orfs true (twopass needs ORF calls from XORF)"
    }

    // rustar-aligner 0.2.0 rejects --outWigStrand Unstranded and only produces stranded
    // signal (str1 + str2); COVERAGE expects the single unstranded bedGraph STAR emits.
    if (params.aligner == 'ruSTAR' && params.star_make_coverage) {
        errors << "  --aligner ruSTAR does not support coverage tracks yet; set --star_make_coverage false"
    }

    // The input scan replaces the checkIfExists that the two read globs in
    // METASSEMBLE cannot use (a run legitimately supplies only one format).
    def has_cbq = false
    if (params.input_dir) {
        def dir = file(params.input_dir)
        def has_fastq = dir.list().any { it ==~ /.*[12]\.f.*q\.gz$/ }
        has_cbq   = dir.list().any { it.endsWith('.cbq') }

        if (!has_fastq && !has_cbq) {
            errors << "  --input_dir '${params.input_dir}' contains no *{1,2}.f*q.gz or *.cbq reads"
        }

        // STAR cannot read CBQ natively. encode_fastqs / native .cbq are decoded
        // after deacon; encode_before_alignment would encode only to decode again.
        if (params.aligner == 'STAR' && params.bqtools_encode_before_alignment) {
            errors << "  --bqtools_encode_before_alignment requires --aligner ruSTAR (STAR cannot read CBQ; encoding just before STAR would be immediately decoded)"
        }
    }

    // Strandedness inference needs CBQ at the pre-trim point; a FASTQ-only run
    // (fastp path) or encode_before_alignment has none and stays 'unstranded'.
    if (params.infer_strandedness) {
        if (!params.bqtools_encode_fastqs && !has_cbq) {
            log.warn "[validateRun] --infer_strandedness: no CBQ reads at the trim point (fastp path); strandedness stays 'unstranded'. Enable with --bqtools_encode_fastqs or native .cbq inputs."
        }
        if (params.bqtools_encode_before_alignment) {
            log.warn "[validateRun] --bqtools_encode_before_alignment encodes after trimming; those samples cannot be sniffed and stay 'unstranded'."
        }
    } else if (params.strand_salmon_index || params.strand_transcriptome) {
        log.warn "[validateRun] --strand_salmon_index/--strand_transcriptome are ignored because --infer_strandedness is false"
    }

    if (params.strand_salmon_index && !file(params.strand_salmon_index).exists()) {
        errors << "  --strand_salmon_index does not exist: '${params.strand_salmon_index}'"
    }
    if (params.strand_transcriptome && !file(params.strand_transcriptome).exists()) {
        errors << "  --strand_transcriptome does not exist: '${params.strand_transcriptome}'"
    }

    errors += validateAnnevo()
    errors += validateTiberius()
    errors += validateXorf()

    if (errors) {
        log.error "Parameter validation failed:\n${errors.join('\n')}"
        System.exit(1)
    }
}

def validateFromPolishing() {
    def errors = []
    if (!params.polish_path)   errors << "  --polish_path is required"
    if (!params.genome)      errors << "  --genome is required"
    if (!params.annotation)  errors << "  --annotation is required"

    if (params.do_twopass_polish && !params.xorf_call_orfs) {
        errors << "  --do_twopass_polish requires --xorf_call_orfs true (twopass needs ORF calls from XORF)"
    }

    errors += validateAnnevo()
    errors += validateTiberius()
    errors += validateXorf()

    if (errors) {
        log.error "Parameter validation failed:\n${errors.join('\n')}"
        System.exit(1)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    LOCAL SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


workflow PIPELINE_COMPLETION {

    take:
    email
    email_on_fail
    plaintext_email
    outdir
    use_mailx
    ch_samplesheet

    main:

    if (params.sent_email) {
        EMAIL_RESULTS (
            email,
            email_on_fail,
            plaintext_email,
            outdir,
            use_mailx,
            ch_samplesheet
        )
    }

    workflow.onError {
        log.error "ERROR: Pipeline failed. Please refer to github issues: https://github.com/alejandrogzi/metassembly/issues"
    }

}

// ── Default: full pipeline ─────────────────────────────────────────────────
workflow FULL_RUN {

    log.info """
    > xasm v${workflow.manifest.version}
    > Meta-assemble bulk RNA-seq datasets at scale
    > The Hiller Lab at the Senckenberg Research Institute

    Authors: ${workflow.manifest.author}
    Github:  ${workflow.manifest.homePage}

      Input     : ${params.input_dir}
      Genome    : ${params.genome}
      Annotation: ${params.annotation}
      ANNEVO    : ${params.annevo_annotation ?: (params.annevo_predict ? "predict (${params.annevo_lineage})" : 'off')}
      TIBERIUS  : ${params.tiberius_annotation && params.tiberius_predict ? "file + predict (${params.tiberius_model_cfg})" : (params.tiberius_annotation ?: (params.tiberius_predict ? "predict (${params.tiberius_model_cfg})" : 'off'))}
      Outdir    : ${params.output_dir}
      Prefix    : ${params.prefix}
      Profile   : ${workflow.profile}
    """.stripIndent()


    validateRun()

    METASSEMBLE (
        params.input_dir,
        params.genome,
        params.annotation,
        params.star_index_path,
        params.star_ignore_gtf_for_index,
        params.deacon_index_path,
        params.deacon_download_index,
        params.deacon_make_single_index,
        params.star_make_coverage,
        params.skip_assembly,
        params.splice_scores_dir,
        params.output_dir
    )

    if (!params.skip_assembly) {
        // NOTE: METASSEMBLE.out.genome is emitted as soon as the 2bit/gz conversion
        // finishes, so ANNEVO and TIBERIUS run in parallel with alignment + assembly.
        ANNEVO_ANNOTATION(METASSEMBLE.out.genome)
        TIBERIUS_ANNOTATION(METASSEMBLE.out.genome)

        POLISH (
            METASSEMBLE.out.metassembly,
            METASSEMBLE.out.genome,
            METASSEMBLE.out.annotation_bed,
            ANNEVO_ANNOTATION.out.bed,
            TIBERIUS_ANNOTATION.out.bed,
            params.repeats,
            params.splice_scores_dir,
            params.do_twopass_polish,
            params.prefix
        )
    }

    if (!params.skip_bb_conversion) {
        ch_autosql = params.autosql ? file(params.autosql, checkIfExists: true) : Channel.of([])

        BEDTOBIGBED_HQ(
            POLISH.out.hq,
            METASSEMBLE.out.chrom_sizes,
            ch_autosql
        )
        BEDTOBIGBED_RETENTION(
            POLISH.out.retentions,
            METASSEMBLE.out.chrom_sizes,
            ch_autosql
        )
        BEDTOBIGBED_TRUNCATIONS(
            POLISH.out.truncations,
            METASSEMBLE.out.chrom_sizes,
            ch_autosql
        )
        BEDTOBIGBED_STRONG_RTS(
            POLISH.out.strong_rts,
            METASSEMBLE.out.chrom_sizes,
            ch_autosql
        )
        BEDTOBIGBED_WEAK_RTS(
            POLISH.out.weak_rts,
            METASSEMBLE.out.chrom_sizes,
            ch_autosql
        )
        BEDTOBIGBED_ARTIFACTS(
            POLISH.out.artifacts,
            METASSEMBLE.out.chrom_sizes,
            ch_autosql
        )
        BEDTOBIGBED_FUSIONS(
            POLISH.out.fusions,
            METASSEMBLE.out.chrom_sizes,
            ch_autosql
        )
        BEDTOBIGBED_SCRAPS(
            POLISH.out.scraps,
            METASSEMBLE.out.chrom_sizes,
            ch_autosql
        )
    }

    PIPELINE_COMPLETION (
        params.email_to,
        params.email_on_fail,
        params.plaintext_email,
        params.output_dir,
        params.use_mailx,
        METASSEMBLE.out.samplesheet
    )
}

// ── Checkpoint: start from polishing step (skip metassembly) ─────────────────────
workflow FROM_POLISHING {
    validateFromPolishing()

    ch_versions = Channel.empty()

    def genome_file = file(params.genome, checkIfExists: true)
    def genome_path = genome_file.toString()

    ch_chrom_sizes = CHROMSIZE([[:], genome_file]).chromsize.map { it[1] }

    // INFO: if fasta is .2bit or .gz, convert or uncompress it
    if (genome_path.endsWith(".2bit")) {
        ch_fasta = TWOBIT_TO_FA([[:], genome_file]).fasta.map { it[1] }
        ch_versions = ch_versions.mix(TWOBIT_TO_FA.out.versions)
    } else if (genome_path.endsWith(".gz")) {
        ch_fasta = GUNZIP_FASTA([[:], genome_file]).gunzip.map { it[1] }
        ch_versions = ch_versions.mix(GUNZIP_FASTA.out.versions)
    } else {
        ch_fasta = Channel.value(genome_file)
    }

    // INFO: preparing annotation 
    if (params.annotation.endsWith('.gz') || params.annotation.endsWith('.gtf') || params.annotation.endsWith('.gff')) {
      Channel.value(file(params.annotation, checkIfExists: true))
        .map { it -> [ [ id: it.baseName ], it ] }
        .set { ch_gtf }

      if (params.annotation.endsWith('.gz')) {
         GUNZIP_GTF(ch_gtf)

         ch_gtf = GUNZIP_GTF.out.gunzip
         ch_versions = ch_versions.mix(GUNZIP_GTF.out.versions)
      }

      // INFO: converting to bed
      GXF2BED(ch_gtf)
      ch_bed = GXF2BED.out.bed

      ch_versions = ch_versions.mix(GXF2BED.out.versions)
    } else if (params.annotation.endsWith('.bed')) {
      // INFO: converting to gtf
      Channel.value(file(params.annotation, checkIfExists: true))
        .map { it -> [ [ id: it.baseName ], it ] }
        .set { ch_bed }

      BED2GTF(
        ch_bed,
        Channel.of([[], []])
      )

      ch_gtf = BED2GTF.out.gtf

      ch_versions = ch_versions.mix(BED2GTF.out.versions)
    } else {
      ch_gtf = Channel.of([:])
      ch_bed = Channel.of([[], []])
    }

    ANNEVO_ANNOTATION(ch_fasta)
    TIBERIUS_ANNOTATION(ch_fasta)

    POLISH (
        Channel.fromPath(params.polish_path, checkIfExists: true)
            .map { it -> [ [ id: it.baseName ], it ] },
        ch_fasta,
        ch_bed,
        ANNEVO_ANNOTATION.out.bed,
        TIBERIUS_ANNOTATION.out.bed,
        params.repeats,
        params.splice_scores_dir,
        params.do_twopass_polish,
        params.prefix
    )

    if (!params.skip_bb_conversion) {
        ch_autosql = params.autosql ? file(params.autosql, checkIfExists: true) : Channel.of([])

        BEDTOBIGBED_HQ(
            POLISH.out.hq,
            ch_chrom_sizes,
            ch_autosql
        )
        BEDTOBIGBED_RETENTION(
            POLISH.out.retentions,
            ch_chrom_sizes,
            ch_autosql
        )
        BEDTOBIGBED_TRUNCATIONS(
            POLISH.out.truncations,
            ch_chrom_sizes,
            ch_autosql
        )
        BEDTOBIGBED_STRONG_RTS(
            POLISH.out.strong_rts,
            ch_chrom_sizes,
            ch_autosql
        )
        BEDTOBIGBED_WEAK_RTS(
            POLISH.out.weak_rts,
            ch_chrom_sizes,
            ch_autosql
        )
        BEDTOBIGBED_ARTIFACTS(
            POLISH.out.artifacts,
            ch_chrom_sizes,
            ch_autosql
        )
        BEDTOBIGBED_FUSIONS(
            POLISH.out.fusions,
            ch_chrom_sizes,
            ch_autosql
        )
        BEDTOBIGBED_SCRAPS(
            POLISH.out.scraps,
            ch_chrom_sizes,
            ch_autosql
        )
    }
}

// ── Checkpoint: start from bigbed step (skip metassembly + polishing) ─────────────────────
workflow FROM_BIGBED {
    // validateFromBigBed()

    ch_versions = Channel.empty()

    def genome_file = file(params.genome, checkIfExists: true)
    def genome_path = genome_file.toString()

    ch_chrom_sizes = CHROMSIZE([[:], genome_file]).chromsize.map { it[1] }
    ch_autosql = params.autosql ? file(params.autosql, checkIfExists: true) : Channel.of([])
    Channel.fromPath("${params.all_bed_path}/*.bed", checkIfExists: true)
        .map { it -> [ [ id: it.baseName ], it ] }
        .set { ch_beds }

    SORT_BED(ch_beds)
    BEDTOBIGBED(
        SORT_BED.out.sorted,
        ch_chrom_sizes,
        ch_autosql.first()
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow XASM {
    if (params.annevo_annotation && params.annevo_predict) {
        log.warn "--annevo_annotation is set — ignoring --annevo_predict"
    }

    if (params.from == "polish") {
        // ── Checkpoint: start from polishing step (skip metassembly) ─────────────────────
        log.info "Resuming from ${params.from} checkpoint — skipping meta-assembly"
        FROM_POLISHING()
    } else if (params.from == "bigbed") {
        // ── Checkpoint: start from bigbed step (skip metassembly + polishing) ─────────────────────
        log.info "Resuming from ${params.from} checkpoint — skipping meta-assembly + polishing"
        FROM_BIGBED()
    }
    else {
        // ── Default: full pipeline ─────────────────────────────────────────────────
        log.info "Starting full pipeline!"
        FULL_RUN()
    }
}

workflow {
    XASM()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION HANDLER
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (workflow.success) {
        def results_dir = new File(params.output_dir as String, '10_final')
        def final_beds = results_dir.exists() ? (results_dir.listFiles()?.findAll { it.name.endsWith('.bed') } ?: []) : []
        if (final_beds && params.skip_bb_conversion) {
            log.info "Pipeline completed successfully!"
            log.info "Final predictions: ${final_beds.collect { it.toString() }.join(', ')}"
        } else {
            if ((params.from == "bigbed") || (!params.skip_bb_conversion)) {
                def bigbed_dir = new File(params.output_dir as String, '11_bbs')
                def final_bbs = bigbed_dir.exists() ? (bigbed_dir.listFiles()?.findAll { it.name.endsWith('.bb') } ?: []) : []
                log.info "Pipeline completed successfully!"
                if (final_bbs) {
                    log.info "Final predictions: ${final_bbs.collect { it.toString() }.join(', ')}"
                }
            } else {
              log.warn "Pipeline reported success but final bed file was not produced - check that all steps ran"
            }
        }
        log.info "Run time   : ${workflow.duration}"
    } else {
        log.error "Pipeline FAILED — ${workflow.errorMessage}"
    }
}

workflow.onError {
    log.error "Pipeline error: ${workflow.errorMessage}"
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
