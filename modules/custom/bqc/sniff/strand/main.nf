/*
Copyright (c) 2026 The Hiller Lab at the Senckenberg Gessellschaft für Naturforschung
Distributed under the terms of the Apache License, Version 2.0.
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    BQC_SNIFF_STRAND — Infer RNA-seq strandedness by mapping CBQ reads against
    a Salmon 2.x transcriptome index.

    Runs on the raw CBQ in parallel with BQC trimming; the JSON report is turned
    into the sample's `strandedness` upstream in METASSEMBLE. bqc classifies as
    forward/reverse/unstranded/undetermined and exits 0 on all four (undetermined
    is a valid answer, not a failure), so only real errors fail the task.

    NOTE: container-only. The bqc image must be built with the `sniff-strand`
    feature (`cargo install bqc --features sniff-strand`).
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process BQC_SNIFF_STRAND {
    tag "$meta.id"
    label 'custom_process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        '' :
        'ghcr.io/hillerlab/bqc:latest' }"

    input:
    tuple val(meta), path(cbq), path(index)

    output:
    tuple val(meta), path("*.strand.json"), emit: report
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    bqc \\
        sniff strand \\
        $cbq \\
        --index $index \\
        --format json \\
        -o ${prefix}.strand.json \\
        -T $task.cpus \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bqc: \$(bqc --version | sed 's/bqc //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo '{"result":{"strandedness":"stub"}}' > ${prefix}.strand.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bqc: \$(bqc --version | sed 's/bqc //g')
    END_VERSIONS
    """
}
