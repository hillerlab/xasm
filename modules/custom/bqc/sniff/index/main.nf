/*
Copyright (c) 2026 The Hiller Lab at the Senckenberg Gessellschaft für Naturforschung
Distributed under the terms of the Apache License, Version 2.0.
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    BQC_SNIFF_INDEX — Build a reusable Salmon 2.x transcriptome index.

    Runs once per pipeline run on the reference CDS FASTA from XLOCI_CDS; the
    resulting directory is staged into every BQC_SNIFF_STRAND task. bqc prints
    the index path to stdout, so the directory is pinned with
    --output-directory/--prefix instead.

    NOTE: container-only. The bqc image must be built with the `sniff-strand`
    feature (`cargo install bqc --features sniff-strand`) — the default image
    build compiles the subcommand out.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process BQC_SNIFF_INDEX {
    tag "$meta.id"
    label 'custom_process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        '' :
        'ghcr.io/hillerlab/bqc:latest' }"

    input:
    tuple val(meta), path(sequence)

    output:
    tuple val(meta), path("${prefix}"), emit: index
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: 'bqc_salmon_index'
    def args = task.ext.args ?: ''
    """
    bqc \\
        sniff index \\
        --sequence $sequence \\
        --output-directory . \\
        --prefix ${prefix} \\
        -T $task.cpus \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bqc: \$(bqc --version | sed 's/bqc //g')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: 'bqc_salmon_index'
    """
    mkdir -p ${prefix}
    touch ${prefix}/info.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bqc: \$(bqc --version | sed 's/bqc //g')
    END_VERSIONS
    """
}
