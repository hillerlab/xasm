/*
Copyright (c) 2026 The Hiller Lab at the Senckenberg Gessellschaft für Naturforschung
Distributed under the terms of the Apache License, Version 2.0.
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    STRINGTIE3 — Transcript assembly and quantification for RNA-seq
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process STRINGTIE3 {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/stringtie:3.0.3--h29c0135_0' :
        'quay.io/biocontainers/stringtie:3.0.3--h29c0135_0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.transcripts.gtf"), emit: transcript_gtf
    tuple val(meta), env(LINE_COUNT)         , emit: assembled_transcripts
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args ?: ''
    def prefix    = task.ext.prefix ?: "${meta.id}"

    def strandedness = ''
    if (meta.strandedness == 'forward') {
        strandedness = '--fr'
    } else if (meta.strandedness == 'reverse') {
        strandedness = '--rf'
    }
    """
    stringtie \\
        $bam \\
        $strandedness \\
        -o ${prefix}.transcripts.gtf \\
        -p $task.cpus \\
        $args

    # StringTie exits 0 and may leave no transcript rows. grep -w exits 1 in
    # that case, and pipefail would fail the task. awk prints 0.
    touch ${prefix}.transcripts.gtf
    LINE_COUNT=\$(awk '\$3 == "transcript" { n++ } END { print n + 0 }' ${prefix}.transcripts.gtf)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stringtie: \$(stringtie --version 2>&1 | sed 's/^.*version //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.transcripts.gtf

    LINE_COUNT=0

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stringtie: \$(stringtie --version 2>&1 | sed 's/^.*version //; s/ .*\$//')
    END_VERSIONS
    """
}
