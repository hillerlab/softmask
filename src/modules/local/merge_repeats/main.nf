process MERGE_REPEATMASKER_REPEATS {
    tag   "$meta.id"
    label 'process_single'

    conda "conda-forge::gawk=5.3.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:5.3.0':
        'quay.io/biocontainers/gawk:5.3.0' }"

    input:
    tuple val(meta), path(files)

    output:
    tuple val(meta), path("${prefix}.repeats.tsv"), emit: merged
    path "versions.yml",                            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Build the file-of-filenames in deterministic order so the global ID
    # numbering is reproducible. -V keeps natural chunk order
    # (..._1-100000 < ..._100001-200000 < ...). `printf` is a bash builtin,
    # so listing the staged inputs here is itself ARG_MAX-safe.
    printf '%s\\n' ${files} | sort -V > ${prefix}.fofn

    merge_repeatmasker_repeats.sh ${prefix}.fofn > ${prefix}.repeats.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version 2>&1 | head -n1)
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "SW_score\tperc_div\tperc_del\tperc_ins\tquery\tbegin\tend\tq_left\tstrand\trepeat\tclass_family\tr_begin\tr_end\tr_left\tID\toverlap\tchrom\tabs_begin\tabs_end" > ${prefix}.repeats.tsv
    touch versions.yml
    """
}
