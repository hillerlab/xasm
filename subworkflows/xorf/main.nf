/*
Copyright (c) 2026 The Hiller Lab at the Senckenberg Gessellschaft für Naturforschung
Distributed under the terms of the Apache License, Version 2.0.
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { XORF } from '../../modules/xorf/src/subworkflows/xorf/main.nf'
include { WGET as WGET_SAMBA_WEIGHTS } from '../../modules/xorf/src/modules/wget/main.nf'
include { WGET as WGET_PROTEIN_DATABASE } from '../../modules/xorf/src/modules/wget/main.nf'
include { UNTAR } from '../../modules/xorf/src/modules/untar/main.nf'
include { GUNZIP as GUNZIP_DATABASE } from '../../modules/xorf/src/modules/gunzip/main.nf'
include { FASTA_MERGE } from '../../modules/xorf/src/modules/diamond/merge/main.nf'
include { DIAMOND_MAKEDB } from '../../modules/xorf/src/modules/diamond/makedb/main.nf'
include { DOWNLOAD_ESMFOLD_WEIGHTS } from '../../modules/xorf/src/modules/esmfold/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    LOCAL SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow XORF_RUN {
    take:
        regions                 // channel: [ meta, bed ] first-pass HQ transcripts
        sequence                // channel: [ meta, fasta ] genome (2bit ok)
        out_subdir              // val: 09_polish subdir for the XORF_PIPELINE_INFO storeDir [default: 'xorf']

    main:
        // Optional take arg: callers omitting it land on the first-pass location.
        def xorf_out_dir = "${params.output_dir}/09_polish/${out_subdir ?: 'xorf'}"

        // Bridge: the vendored BLAST module reads global params.esm for its HF env
        // exports, but xasm exposes xorf_esm. Mirror the flag so the env exports
        // and the --esm ext.args flag (nextflow.config) always agree. Idempotent
        // across both XORF_RUN invocations.
        params.esm = params.xorf_esm ?: false
        /*
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            SAMBA WEIGHTS
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */

        ch_samba_weights = Channel.empty()
        if (params.xorf_samba_local_weights) {
            ch_samba_weights = Channel.value(
              file(params.xorf_samba_local_weights, checkIfExists: true)
            ).map { path -> [ [id : path.baseName ], path ] }
        } else {
            WGET_SAMBA_WEIGHTS(
              Channel.value(
                params.xorf_samba_weights
              ).map { url -> [ [id : url.tokenize('/')[-1]], url ] }
            )
            ch_samba_weights = WGET_SAMBA_WEIGHTS.out.outfile
        }

        /*
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            ESMFOLD WEIGHTS
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */

        // Mirrors the standalone xorf ESM weights preparation (xorf src/main.nf).
        // The channel carries a bare path: the HF hub cache dir when ESM is on, a
        // dummy file otherwise (BLAST stages it via `each`, so it must not be empty).
        // ponytail: downloaded once per XORF_RUN invocation (same as samba weights);
        // hoist to POLISH if twopass+ESM double downloads ever matter.
        ch_esm_weights = Channel.empty()
        ch_esm_versions = Channel.empty()
        if (params.xorf_esm) {
            if (params.xorf_esmfold_local_weights) {
                ch_esm_weights = Channel.value(
                  file(params.xorf_esmfold_local_weights, checkIfExists: true)
                )
            } else {
                DOWNLOAD_ESMFOLD_WEIGHTS(
                  file("${projectDir}/modules/xorf/modules/orf/src/esmfold2_fast.py", checkIfExists: true)
                )
                ch_esm_weights = DOWNLOAD_ESMFOLD_WEIGHTS.out.cache
                ch_esm_versions = DOWNLOAD_ESMFOLD_WEIGHTS.out.versions
            }
        } else {
            ch_esm_weights = Channel.value(file("${projectDir}/params.json"))
        }

        /*
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            PROTEIN DATABASE
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */

        // Mirrors the standalone xorf database preparation (xorf src/main.nf).
        // BLAST declares `each path(database)`, so ch_database carries bare files.
        ch_database = Channel.empty()
        ch_database_versions = Channel.empty()
        if (params.xorf_custom_database) {
            if (params.xorf_custom_database.endsWith('.dmnd')) {
                ch_database = Channel.fromPath(params.xorf_custom_database, checkIfExists: true)
            } else if (params.xorf_custom_database.endsWith('.dmnd.gz')) {
                GUNZIP_DATABASE(
                    Channel.value(
                        [ [id: params.xorf_custom_database.tokenize('/')[-1]], params.xorf_custom_database ]
                    )
                )
                GUNZIP_DATABASE.out.gunzip
                    .map { meta, it -> it }
                    .set { ch_database }
            } else if (params.xorf_custom_database.endsWith('.fa') || params.xorf_custom_database.endsWith('.fasta')
                        || params.xorf_custom_database.endsWith('.fa.gz') || params.xorf_custom_database.endsWith('.fasta.gz')) {
                WGET_PROTEIN_DATABASE(
                    Channel.value('https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz')
                    .map { it -> [ [ id: 'uniprot_sprot.fasta.gz' ], it ] }
                )
                FASTA_MERGE(
                    Channel.fromPath(params.xorf_custom_database, checkIfExists: true)
                        .map { it -> [ [ id: it.baseName ], it ] }
                        .combine(WGET_PROTEIN_DATABASE.out.outfile.map { meta, raw -> raw })
                )
                DIAMOND_MAKEDB(
                    FASTA_MERGE.out.fasta.map { meta, fasta -> [ [ id: 'merged_database' ], fasta ] },
                    [],
                    [],
                    []
                )
                DIAMOND_MAKEDB.out.db
                    .map { meta, it -> it }
                    .set { ch_database }
                ch_database_versions = FASTA_MERGE.out.versions.mix(DIAMOND_MAKEDB.out.versions)
            } else {
                error """
                ERROR: custom_database extension not recognized.
                Please provide a custom database in one of these formats:
                  .dmnd/.dmnd.gz   -> replaces the default database entirely
                  .fa/.fasta/.fa.gz/.fasta.gz -> appended to the default SwissProt database and reindexed
                """.stripIndent()
            }
        } else {
            WGET_PROTEIN_DATABASE(
                Channel.value('https://zenodo.org/records/21399231/files/diamond_db.tar.gz?download=1')
                    .map { it -> [ [ id: 'uniprot_sprot.tar.gz' ], it ] }
            )
            UNTAR(WGET_PROTEIN_DATABASE.out.outfile)
            ch_database = UNTAR.out.contents.map { meta, it -> it }
        }

        /*
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            MAIN WORKFLOW
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */

        // meta.id carries the source transcript id (`<t>_clean`); XORF keeps it as the
        // group key and prefixes the renamed ORF ids with it. `chr` is only a grouping
        // hash: CHUNKER uses it as its prefix and the subworkflow groups by
        // `${meta.id}@${meta.chr}`. Using meta.id keeps every XORF-derived filename
        // deterministic (required by the e2e golden checks).
        XORF(
            regions.map { meta, f -> [ meta + [ chr: meta.id ], f ] },
            sequence.map { meta, fa -> fa },   // XORF takes a bare-file sequence channel
            ch_database,
            xorf_out_dir,
            params.xorf_chunk_size ?: 20,
            ch_samba_weights,
            ch_esm_weights,
            params.xorf_predict_keep_raw,
            params.xorf_selenocysteine_sites,                          
            // Keep this positional order in sync with modules/xorf/src/subworkflows/xorf/main.nf.
            // v0.1.2 added skip_transaid after skip_netstart and strip_truncations after
            // do_polishing; v0.2.0 keeps the same 19-input signature.
            params.xorf_skip_netstart,
            params.xorf_skip_transaid,
            params.xorf_rename_deactivate,
            params.xorf_do_polishing,
            params.xorf_strip_truncations,
            params.xorf_skip_joined_concat,
            false, null, null,             // run_only_on / run_only_mode / run_only_target — masking not ported
            ch_database_versions.mix(ch_esm_versions)
        )

    emit:
        files        = XORF.out.files        // [ meta, bed, tsv ] stripped HQ when xorf_do_polishing
        truncations  = XORF.out.truncations  // [ meta, bed ] 3'UTR truncation discards
        duplicates   = XORF.out.duplicates   // [ meta, bed ] #DU discards
        counts       = XORF.out.counts
        versions     = XORF.out.versions
}
