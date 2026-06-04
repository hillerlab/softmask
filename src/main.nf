#!/usr/bin/env nextflow

/*
Copyright (c) 2026 The Hiller Lab at the Senckenberg Gessellschaft für Naturforschung
Distributed under the terms of the Apache License, Version 2.0.
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    softmask

    Lighweight workflow for soft-masking genomes
    Authors: Alejandro Gonzales-Irribarren, Michael Hiller

    GitHub:  https://github.com/hillerlab/softmask
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT WORKFLOWS AND MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SOFTMASK as MASK_GENOME } from './subworkflows/local/softmask.nf'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    HELP
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

if (params.help) {
    log.info """
    softmask v${workflow.manifest.version}
    Lightweight workflow for soft-masking genomes

    Authors: ${workflow.manifest.author}
    Github:  ${workflow.manifest.homePage}

    Usage (full run):
        nextflow run main.nf \\
            --genome             PATH    Genome file (FASTA or .2bit) \\
            --assembly_prefix    STRING  Genome assembly prefix (e.g. hg38) \\
            --repeat_library     PATH    Repeat library (FASTA or .2bit) \\
            --chunks             INT     Number of chunks to split genome into \\
            --output_format      STRING  Output format (fasta, fasta.gz, 2bit) \\
            --outdir             PATH    Output directory (default: ./results) \\
            --use_container      BOOL    Use container (default: false) \\
            --help               BOOL    Show this message


    Pass all parameters from a JSON file:
        nextflow run main.nf -params-file my_params.json

    Required parameters:
        --genome             PATH    Genome file (FASTA or .2bit)
        --assembly_prefix    STRING  Genome assembly prefix (e.g. hg38)
        --repeat_library     PATH    Repeat library (FASTA or .2bit)
        --chunks             INT     Number of chunks to split genome into
        --output_format      STRING  Output format (fasta, fasta.gz, 2bit)

    Optional parameters:
        --outdir             PATH    Output directory (default: ./results)
        --use_container      BOOL    Use container (default: false)
        --help               BOOL    Show this message

    Profiles:
        local       Run on local machine (default)
        slurm       Submit jobs to SLURM cluster
        conda       Use conda environments
        apptainer   Use Apptainer containers
        singularity Use Singularity containers
        docker      Use Docker containers
        test        Run with bundled test data

    Use --help to show this message.
    """.stripIndent()
    System.exit(0)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    SOFTMASK()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ENTRY WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SOFTMASK {
    log.info """
    > softmask v${workflow.manifest.version}
    > Lighweight workflow for soft-masking genomes

    Authors: ${workflow.manifest.author}
    Github:  ${workflow.manifest.homePage}

      Genome:    ${params.genome}
      Output:    ${params.output_format}
      Outdir:    ${params.outdir}
      Profile:   ${workflow.profile}
    """.stripIndent()

    MASK_GENOME(
        params.genome,
        params.assembly_prefix,
        params.repeat_library,
        params.chunks,
        params.output_format,
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION HANDLER
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (workflow.success) {
        def masked_genome = file("${params.outdir}/05_final/${params.assembly_prefix}.masked.${params.output_format}")
        log.info "Pipeline completed successfully!"
        if (masked_genome.exists()) {
            log.info "Masked genome: ${masked_genome}"
        } else {
            log.warn "Masked genome not found: ${masked_genome}"
        }
        log.info "Run time   : ${workflow.duration}"
    } else {
        log.error "Pipeline FAILED — ${workflow.errorMessage}"
    }
}

workflow.onError {
    log.error "Pipeline error: ${workflow.errorMessage}"
}
