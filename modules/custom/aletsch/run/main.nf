// Copyright (c) 2025 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

process ALETSCH {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/aletsch:1.1.3--hdbdd923_0' :
        'biocontainers/aletsch:1.1.3--hdbdd923_0' }"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.gtf")       , emit: gtf
    tuple val(meta), path("*profile")    , emit: profile
    tuple val(meta), env(LINE_COUNT)     , emit: assembled_transcripts
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Aletsch's third sample-info field is sequencing protocol, not strandedness.
    def library_type = meta.single_end ? 'single_end' : 'paired_end'
    def fallback_insert_size = params.aletsch_fallback_insert_size ?: ''
    def fallback_insert_std = params.aletsch_fallback_insert_std ?: 10
    def chromosome = meta.chr ? "-l ${meta.chr}" : ''

    """
    # Create necessary directories
    mkdir -p ${prefix}_profile
    mkdir -p ${prefix}_gtf

    # Create sample info file
    echo "${bam}\t${bai}\t${library_type}" > ${prefix}.info

    # Run Aletsch profile generation
    aletsch \\
        --profile \\
        -i ${prefix}.info \\
        -p ${prefix}_profile \\
        $args

    # Compact per-chromosome inputs may not meet Aletsch's internal preview
    # sample minimum. Use an explicit known insert size only if profiling left
    # its mean at zero.
    if [ -n "$fallback_insert_size" ]; then
        for profile in ${prefix}_profile/*.profile; do
            awk -v mean="$fallback_insert_size" -v std="$fallback_insert_std" '
                \$1 == "insertsize_ave" && \$2 == 0 { \$2 = mean }
                \$1 == "insertsize_std" && \$2 == 0 { \$2 = std }
                { print }
            ' "\$profile" > "\$profile.tmp"
            mv "\$profile.tmp" "\$profile"
        done
    fi

    # Run Aletsch assembly
    aletsch \\
        -i ${prefix}.info \\
        -o ${prefix}_gtf/${prefix}.gtf \\
        -p ${prefix}_profile \\
        -d ${prefix}_gtf \\
        $chromosome \\
        $args

    # Move output to current directory (Aletsch can exit 0 without writing a GTF
    # on chromosomes with (almost) no signal; keep an empty GTF flowing so
    # RENAME/BEAVER see LINE_COUNT=0 instead of failing on a missing file)
    if [ -f "${prefix}_gtf/${prefix}.gtf" ]; then
        mv ${prefix}_gtf/${prefix}.gtf ${prefix}.gtf
    else
        echo "[WARN] ${meta.id}: aletsch produced no ${prefix}.gtf; emitting empty GTF" >&2
        touch ${prefix}.gtf
    fi
    mv ${prefix}_profile ${prefix}.profile

    LINE_COUNT=\$(awk '\$3 == "transcript" { n++ } END { print n + 0 }' ${prefix}.gtf)

    rm -rf ${prefix}_gtf/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        aletsch: \$(aletsch --version 2>&1 | sed 's/^.*aletsch //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.gtf

    LINE_COUNT=0

    mkdir -p ${prefix}.profile
    touch ${prefix}.profile/0.profile

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        aletsch: \$(aletsch --version 2>&1 | sed 's/^.*aletsch //; s/ .*\$//')
    END_VERSIONS
    """
}
