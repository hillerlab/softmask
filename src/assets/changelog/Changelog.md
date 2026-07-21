# Changelog

All notable changes to the softmask pipeline are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [0.0.2] - 2026-07-21

### Added

- **Test suite** covering core pipeline components:
  - `test_group_mask_chunks.py` — unit tests for the chunk-based masking logic
  - `test_merge_repeatmasker_repeats.py` — validation of repeat merging behaviour
  - `test_pipeline_ambiguity.py` — end-to-end ambiguity masking workflow tests
- **Test configuration** (`pipeline_ambiguity.config`) to drive the ambiguity test pipeline
- **Test data** for ambiguity resolution (`ambiguity_genome.fa`, `ambiguity_repeats.fa`)
- **CI pipeline** (`.github/workflows/tests.yml`) running tests on every push

### Changed

- **Docker image builder** overhauled to improve reproducibility and reduce layer count
- **Test genome** (`test_genome.fa`) drastically reduced from ~11 800 lines to ~2 000 lines, keeping only the sequences relevant to unit testing

### Removed

- **Redundant test genome sequences** — the bulk of `test_genome.fa` was trimmed to minimise repository size while preserving all test coverage

---

## [0.0.1] - 2026-06-04

### Added

- **Initial pipeline release** with the following capabilities:
  - Genome soft-masking via RepeatMasker and RepeatModeler
  - Scaffold renaming and header normalisation
  - Genome format conversion (FASTA ↔ 2bit)
  - Chunk-based masking with configurable window size
  - Merge and reconciliation of RepeatMasker outputs
  - Gzip/gunzip handling for intermediate files
- **Nextflow DSL2 modules** for each pipeline step
- **Main workflow** (`main.nf`) and **subworkflow** (`softmask.nf`) wiring together all stages
- **Workflow chart** (`softmask.mermaid`) documenting the pipeline topology
- **Base Docker image** and build configuration
- **HPC launcher script** (`do_softmask.sh`)
- **Configuration** (`nextflow.config`) with profiles for local and cluster execution
- **Parameter schema** (`params.json`) for pipeline input validation
- **Project boilerplate** — MIT license, `.gitignore`, and this changelog
- **Branding assets** — Hiller Lab logo for documentation
