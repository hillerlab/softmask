process REPEATMASKER {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/repeatmasker:4.1.5--pl5321hdfd78af_0':
        'quay.io/biocontainers/repeatmasker:4.1.5--pl5321hdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)
    tuple val(meta1), path(lib)

    output:
    tuple val(meta), path("${prefix}.masked")   , emit: masked
    tuple val(meta), path("${prefix}.out")      , emit: out
    tuple val(meta), path("${prefix}.tbl")      , emit: tbl, optional: true
    tuple val(meta), path("${prefix}.gff")      , emit: gff, optional: true
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args     ?: ''
    prefix      = task.ext.prefix   ?: "${meta.id}"
    def lib_arg = lib               ? "-lib $lib"   : ''

    def out_fasta    = fasta.getBaseName(fasta.name.endsWith('.gz') ? 1 : 0)
    def fasta_gz_cmd = fasta.name.endsWith('.gz') ? "gunzip -c ${fasta} > ${out_fasta}" : ""

    """
    ${fasta_gz_cmd}
    RepeatMasker \\
        $lib_arg \\
        -pa ${task.cpus} \\
        -dir ${prefix} \\
        --xsmall \\
        ${args} \\
        ${out_fasta}

    if [[ -f $prefix/${out_fasta}.masked ]]; then
        mv $prefix/${out_fasta}.masked ${prefix}.masked
    else
        # RepeatMasker may omit .masked when it finds no repetitive sequence.
        # Preserve the original chunk so downstream concatenation cannot alter
        # the assembly coordinate system.
        cp ${out_fasta} ${prefix}.masked
    fi

    mv $prefix/${out_fasta}.out ${prefix}.out

    if [[ -f $prefix/${out_fasta}.tbl ]]; then
        mv $prefix/${out_fasta}.tbl ${prefix}.tbl
    fi

    if [[ -f $prefix/${out_fasta}.out.gff ]]; then
        mv $prefix/${out_fasta}.out.gff ${prefix}.gff
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        repeatmasker: \$(RepeatMasker -v | sed 's/RepeatMasker version //1')
    END_VERSIONS
    """

    stub:
    prefix          = task.ext.prefix       ?: "${meta.id}"
    def args        = task.ext.args         ?: ''
    def touch_gff   = args.contains('-gff') ? "touch ${prefix}.gff" : ''

    """
    touch ${prefix}.masked
    touch ${prefix}.out
    touch ${prefix}.tbl
    $touch_gff

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        repeatmasker: \$(RepeatMasker -v | sed 's/RepeatMasker version //1')
    END_VERSIONS
    """
}
