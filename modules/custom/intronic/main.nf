process INTRONIC {
    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        '' :
        'ghcr.io/alejandrogzi/isox-py:latest' }"

    input:
    tuple val(meta), path(introns)

    output:
    tuple val(meta), path("*.meta.iic")      , optional: true, emit: iic
    path  "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args          = task.ext.args   ?: ''
    def prefix        = task.ext.prefix ?: "${meta.id}"
    """
    if [ ! -s "$introns" ]; then
        touch ${prefix}.meta.iic
    else
        # xloci repeats a genomic intron once per isoform. -q loads every row,
        # and iso-classify keeps one label per coordinate.
        awk -F '\\t' 'NF && !seen[\$1]++' "$introns" > ${prefix}.unique.iic

        intronIC classify \\
          -q ${prefix}.unique.iic \\
          -n ${prefix} \\
          $args \\
          -p ${task.cpus}

        if ! compgen -G "*.meta.iic" > /dev/null; then
            touch ${prefix}.meta.iic
        fi
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        intronIC: \$( intronIC --version | head -n 1 | sed 's/intronIC //g' | sed 's/ (.*//g' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.meta.iic

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        intronIC: \$( intronIC --version | head -n 1 | sed 's/intronIC //g' | sed 's/ (.*//g' )
    END_VERSIONS
    """
} 
