process CONCAT_FASTA {
    tag "${meta.id}"
    label 'process_single'

    conda "conda-forge::python=3.11"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11' :
        'biocontainers/python:3.11' }"

    input:
    tuple val(meta), path(files, stageAs: "seqs/*")

    output:
    tuple val(meta), path("*.fasta"), emit: fasta
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix       = task.ext.prefix ?: "${meta.id}"
    def args     = task.ext.args ?: ''
    """
    group_mask_chunks.py \\
        ${args} \\
        --input seqs/ \\
        --output ${prefix}.fasta \\
        --verbose

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //')
        group_mask_chunks: \$(group_mask_chunks.py --version 2>&1 | sed 's/^.* //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //')
        group_mask_chunks: \$(group_mask_chunks.py --version 2>&1 | sed 's/^.* //')
    END_VERSIONS
    """
}
