/*
Copyright (c) 2026 The Hiller Lab at the Senckenberg Gessellschaft für Naturforschung
Distributed under the terms of the Apache License, Version 2.0.
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SOFTMASK subworkflow
    Performs soft-masking of genome repeats

    Steps (both conditional on params flags):
    1. CONVERSION: Convert input genome to .fa
    2. RE-FORMAT HEADERS: Reformat headers to UCSC format
    3. FXSPLIT: Split into chunks + get rid of soft-masking at the same time
    4. REPEATMODELER_DATABASE: Run RepeatModeler on genome (whole-genome) if 'repeat_library = null'
    5. REPEEATMODELER_MODEL: Run RepeatModeler on genome
    6. REPEATMASKER: Run RepeatMasker on chunks determined by params.chunks
    7. MERGE_REPEATMASKER_REPEATS: Merge per-chunk RepeatMasker .out files into a single tab-separated table.
    8. CONCAT_FASTA: Concatenate masked chunks and produce single output

    Emits: soft-masked genome -> *.masked.{fasta,fasta.gz,2bit}
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { TWOBIT_TO_FA } from '../../modules/local/ucsc/twobittofa/main'
include { GUNZIP as GUNZIP_FASTA } from '../../modules/local/gunzip/main'
include { RENAMEHEADERS as REFORMAT_GENOME_HEADERS } from '../../modules/local/renameheaders/main'
include { FXSPLIT } from '../../modules/local/fxsplit/main'
include { REPEATMASKER } from '../../modules/local/repeatmasker/main'
include { REPEATMODELER_BUILDDATABASE } from '../../modules/local/repeatmodeler/builddatabase/main'
include { REPEATMODELER_MODEL } from '../../modules/local/repeatmodeler/model/main'
include { CONCAT_FASTA } from '../../modules/local/concat/main'
include { GAWK_JOIN as GAWK_JOIN_FASTA_FILES } from '../../modules/local/gawk/join/main'
include { FA_TO_TWO_BIT } from '../../modules/local/ucsc/fatotwobit/main'
include { GZIP as GZIP_FASTA } from '../../modules/local/gzip/main'
include { PUBLISH as PUBLISH_FASTA } from '../../modules/local/publish/main'
include { MERGE_REPEATMASKER_REPEATS } from '../../modules/local/merge_repeats/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SOFTMASK {
    take:
      genome           // str: path to genome
      assembly_prefix  // str: optional; pass [] (or '') to fall back to meta.id
      repeat_library   // str: path to repeat library
      chunks           // int: number of chunks to split genome into
      output_format    // str: output format (fasta, fasta.gz, 2bit)

    main:
        ch_versions = Channel.empty()

        def genome_file = file(genome, checkIfExists: true)
        def genome_path = genome_file.toString()

        // INFO: if fasta is .2bit or .gz, convert or uncompress it
        ch_fasta = Channel.empty()
        if (genome_path.endsWith(".2bit")) {
            ch_fasta = TWOBIT_TO_FA([[:], genome_file]).fasta.map { it[1] }
            ch_versions = ch_versions.mix(TWOBIT_TO_FA.out.versions)
        } else if (genome_path.endsWith(".gz")) {
            ch_fasta = GUNZIP_FASTA([[:], genome_file]).gunzip.map { it[1] }
            ch_versions = ch_versions.mix(GUNZIP_FASTA.out.versions)
        } else {
            ch_fasta = Channel.value(genome_file)
        }

        REFORMAT_GENOME_HEADERS(
            ch_fasta.map { it -> [ [id: assembly_prefix], it ] },
            assembly_prefix
        )
        ch_versions = ch_versions.mix(REFORMAT_GENOME_HEADERS.out.versions)

        ch_repeats = Channel.empty()
        if (!repeat_library) {
            REPEATMODELER_BUILDDATABASE(
                REFORMAT_GENOME_HEADERS.out.fasta,
            )

            REPEATMODELER_MODEL(
                REPEATMODELER_BUILDDATABASE.out.db,
            )

            ch_repeats = REPEATMODELER_MODEL.out.fasta
            ch_versions = ch_versions.mix(REPEATMODELER_MODEL.out.versions)
            ch_versions = ch_versions.mix(REPEATMODELER_BUILDDATABASE.out.versions)
        } else {
            Channel.value(repeat_library)
              .map { it -> [ [id: assembly_prefix + '_repeats'], it ] }
              .set { ch_repeats }
        }

        FXSPLIT(
            REFORMAT_GENOME_HEADERS.out.fasta,
        )
        FXSPLIT.out.fastx
            .flatMap { meta, chunks ->
                def files = chunks instanceof List ? chunks : [chunks]

                files.withIndex().collect { it, idx ->
                    def name = it.baseName

                    /*
                     * Remove trailing chunk coordinate block:
                     *
                     *   chr1_additionalstring_100-150_assembly_name
                     *   -> chr1_additionalstring
                     *
                     * If no coordinate block exists, keep the full name:
                     *
                     *   mock1r_test_assembly
                     *   -> mock1r_test_assembly
                     */
                    def header = name.replaceFirst(/_\d+-\d+(?:_.*)?$/, '')

                    [
                        [
                            id: name,
                            header: header,
                            chunk: idx
                        ],
                        it
                    ]
                }
            }
            .set { ch_chunks }
        ch_versions = ch_versions.mix(FXSPLIT.out.versions)

        REPEATMASKER(
            ch_chunks,
            ch_repeats,
        )
        ch_versions = ch_versions.mix(REPEATMASKER.out.versions)

        REPEATMASKER.out.masked
          .map { meta, masked_chunk -> [ meta.header, meta, masked_chunk ] }    
          .groupTuple()                                      
          .map { chr, metas, files ->
              [ [ id: chr, n_chunks: metas.size() ], files ]
          }
          .set { ch_grouped_masked_chunks }
        CONCAT_FASTA(
            ch_grouped_masked_chunks
        )
        CONCAT_FASTA.out.fasta
          .map { meta, fasta -> fasta }
          .collect()
          .map { fastas -> [ [id: assembly_prefix + '.masked' ], fastas ] }
          .set { ch_concat_fasta }
        GAWK_JOIN_FASTA_FILES(
            ch_concat_fasta,
            'fasta'
        )

        ch_versions = ch_versions.mix(CONCAT_FASTA.out.versions)
        ch_versions = ch_versions.mix(GAWK_JOIN_FASTA_FILES.out.versions)

        ch_masked_genome = Channel.empty()
        if (output_format == '2bit') {
            FA_TO_TWO_BIT(GAWK_JOIN_FASTA_FILES.out.output)
            ch_masked_genome = FA_TO_TWO_BIT.out.twobit
            ch_versions = ch_versions.mix(FA_TO_TWO_BIT.out.versions)
        } else if (output_format == 'gz') {
            GZIP_FASTA(GAWK_JOIN_FASTA_FILES.out.output)
            ch_masked_genome = GZIP_FASTA.out.gzip
            ch_versions = ch_versions.mix(GZIP_FASTA.out.versions)
        } else {
            PUBLISH_FASTA(GAWK_JOIN_FASTA_FILES.out.output)
            ch_masked_genome = GAWK_JOIN_FASTA_FILES.out.output
        }

        REPEATMASKER.out.out
          .map { meta, rp -> rp }
          .collect()
          .map { repeats -> [ [id: assembly_prefix], repeats ] }
          .set { ch_repeatmasker_repeats }
        MERGE_REPEATMASKER_REPEATS(ch_repeatmasker_repeats)

      emit:
        masked_genome = ch_masked_genome
        versions      = ch_versions
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
