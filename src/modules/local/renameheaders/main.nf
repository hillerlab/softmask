process RENAMEHEADERS {
    tag "${meta.id}"
    label 'process_single'

    conda "conda-forge::python=3.11"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11' :
        'biocontainers/python:3.11' }"

    input:
    tuple val(meta), path(genome)
    val  assembly_prefix          // optional; pass [] (or '') to fall back to meta.id

    output:
    tuple val(meta), path("${prefix}.fasta{,.gz}"), emit: fasta
    tuple val(meta), path("${prefix}.rename.tsv") , emit: dict
    path "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // `prefix` is declared without `def` so it is visible to the output block.
    prefix       = task.ext.prefix ?: "${meta.id}"
    def args     = task.ext.args ?: ''
    def asm      = (assembly_prefix?.toString()?.trim()) ?: prefix
    // Output is `<prefix>.fasta`, or `<prefix>.fasta.gz` when `-G`/`--gzip` is
    // in ext.args (the script appends `.gz` itself). The brace-glob above
    // captures either form, so we never have to parse args here.
    """
    renameheaders.py \\
        ${genome} \\
        --output ${prefix}.fasta \\
        --dict ${prefix}.rename.tsv \\
        --prefix ${asm} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //')
        renameheaders: \$(renameheaders.py --version 2>&1 | sed 's/^.* //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.fasta
    touch ${prefix}.rename.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //')
        renameheaders: \$(renameheaders.py --version 2>&1 | sed 's/^.* //')
    END_VERSIONS
    """
}
