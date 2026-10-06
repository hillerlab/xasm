process BEAVER {
    tag "$meta.id"
    label 'process_long_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        '' : 
        'ghcr.io/alejandrogzi/beaver:latest' }"

    input:
    tuple val(meta), path(gtfs)

    output:
    tuple val(meta), path("beaver_output/*gtf"), optional: true, emit: gtf
    tuple val(meta), path("beaver_output/*csv"), optional: true, emit: csv
    path "*.txt"                               , emit: gtf_list
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "beaver_output"
    """
    for gtf in ${gtfs}; do
        echo "\$gtf" >> ${prefix}.gtf_list.txt
    done

    # Create output directory
    mkdir -p beaver_output

    # Assemble only when there is something to assemble: a per-chromosome run with
    # no transcripts (Aletsch assembled nothing) must skip beaver rather than die
    # on the unconditional mv below.
    if awk -F'\\t' '\$3 == "transcript" { f = 1 } END { exit !f }' ${gtfs}; then
        # Run Beaver (a non-zero exit still fails the task; only exit-0 with
        # missing outputs is tolerated below)
        beaver \\
            ${prefix}.gtf_list.txt \\
            ${prefix} \\
            -t ${task.cpus} \\
            $args

        # Beaver can exit 0 yet emit nothing on sparse contigs (e.g. chrUn
        # scaffolds where every candidate is filtered). Only move what exists;
        # both outputs are optional:true so an empty beaver_output simply means
        # "no metassembly for this chromosome" downstream.
        if [ -f "${prefix}.gtf" ]; then
            mv ${prefix}.gtf beaver_output/
        else
            echo "[WARN] ${meta.id}: beaver produced no ${prefix}.gtf (all candidates filtered?); skipping" >&2
        fi
        if [ -f "${prefix}_feature.csv" ]; then
            mv ${prefix}_feature.csv beaver_output/
        else
            echo "[WARN] ${meta.id}: beaver produced no ${prefix}_feature.csv; skipping" >&2
        fi
    else
        echo "[WARN] ${meta.id}: no assembled transcripts in input (likely too few reads on this chromosome); skipping beaver" >&2
    fi

    if [ ${params.beaver_keep_aletsch_gtf} == false ]; then
        for gtf in ${gtfs}; do
            if [ -L "\$gtf" ]; then
                realpath=\$(readlink -f "\$gtf")
                rm -f "\$gtf"
                if [ -n "\$realpath" ]; then
                    rm -f "\$realpath"
                fi
            else
                rm -f "\$gtf"
            fi
        done
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        beaver: \$(beaver --version 2>&1 | sed 's/^.*beaver //; s/ .*\$//' || echo "0.0.1")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "beaver_output"
    """
    mkdir -p beaver_output
    touch ${prefix}.gtf_list.txt
    touch beaver_output/${prefix}.gtf
    touch beaver_output/${prefix}.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        beaver: 0.0.1
    END_VERSIONS
    """
}
