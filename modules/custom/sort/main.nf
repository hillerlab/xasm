process SORT_BED {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-22105f082fdf56207f1dcc5b6da71a14394e28d7:387d955c0a2cdb831ec519d636e4ffd7062d6ae1-0':
        'quay.io/biocontainers/mulled-v2-22105f082fdf56207f1dcc5b6da71a14394e28d7:387d955c0a2cdb831ec519d636e4ffd7062d6ae1-0' }"

    input:
    tuple val(meta), path(bed)

    output:
    tuple val(meta), path("*.sorted.bed"),   optional: true, emit: sorted
    path "versions.yml",                      emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # bigtools compares chromosome names as raw bytes.
    LC_ALL=C sort -k1,1 -k2,2n -k3,3n ${bed} > ${prefix}.sorted.bed

    if [[ ! -s ${prefix}.sorted.bed ]]; then
        rm ${prefix}.sorted.bed
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sort: \$(sort --version | sed 's/sort (GNU coreutils) //; s/ Copyright.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.sorted.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sort: \$(sort --version | sed 's/sort (GNU coreutils) //; s/ Copyright.*\$//')
    END_VERSIONS
    """
}
