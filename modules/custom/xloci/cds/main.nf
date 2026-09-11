/*
Copyright (c) 2026 The Hiller Lab at the Senckenberg Gessellschaft für Naturforschung
Distributed under the terms of the Apache License, Version 2.0.
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    XLOCI_CDS — Extract reference CDS sequences as a transcriptome FASTA.

    Feeds BQC_SNIFF_INDEX, which builds the reusable Salmon index used by
    BQC_SNIFF_STRAND. xloci is format-agnostic for both inputs (2bit/FASTA
    genome, BED/GTF/GFF regions), so the raw reference annotation is passed
    straight through.

    NOTE: container-only; the xloci image is pulled from ghcr.io.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process XLOCI_CDS {
    tag "$meta.id"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        '' :
        'ghcr.io/alejandrogzi/xloci:latest' }"

    input:
    tuple val(meta), path(genome), path(regions)

    output:
    tuple val(meta), path("*.fa") , emit: fasta
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    xloci \\
        -f cds \\
        -s $genome \\
        -r $regions \\
        -o . \\
        -p ${prefix} \\
        --unmask \\
        -t $task.cpus \\
        --ignore-errors \\
        $args

    # A CDS-less annotation (or a BED without CDS coordinates) leaves nothing to
    # write; fail here rather than handing an empty FASTA to the index builder.
    if [ ! -s ${prefix}.fa ]; then
        echo "ERROR: xloci extracted no CDS sequences from '${regions}'; strandedness inference needs a CDS-bearing annotation" >&2
        exit 1
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        xloci: \$( xloci --version | head -n 1 | sed 's/xloci //g' | sed 's/ (.*//g' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        xloci: \$( xloci --version | head -n 1 | sed 's/xloci //g' | sed 's/ (.*//g' )
    END_VERSIONS
    """
}
