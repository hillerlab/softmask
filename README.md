<p align="center">
  <p align="center">
    <img width=200 align="center" src="./src/assets/figures/hillerlab.png" >
  </p>

  <span>
    <h1 align="center">
        softmask
    </h1>
  </span>

  <p align="center">
    <a href="https://github.com/hillerlab/softmask" reference="_blank">
      <img alt="GitHub License" src="https://img.shields.io/github/license/hillerlab/softmask?color=blue">
    </a>
  </p>

  <p align="center">
    <samp>
        <span> lightweight workflow for soft-masking genomes </span>
        <br>
        <span> The Hiller Lab at the Senckenberg Research Institute </span>
        <br>
        <br>
        <a href="https://training.galaxyproject.org/training-material/topics/genome-annotation/tutorials/repeatmasker/tutorial.html">masking</a> .
        <a href="https://github.com/hillerlab/softmask/blob/main/assets/pipeline/softmask.mermaid">pipeline</a> .
        <a href="https://hillerlab.com/">us</a> 
    </samp>
  </p>

</p>

---

<div align="center">

<pre style="font-size: 18px;">
CAGATGATGATGATGATGATGATGATGAGCTT
█████████████████████░░░░░░░░░░░
CAGATGATGATGATGATGATgatgatgagctt
└────────unique────┘└──repeat──┘
</pre>

</div>

---

> [!IMPORTANT]
> - **Masking**: This pipeline is designed to mask genomes with soft-masking not hard-masking. This framework is compatible with our whole-genome alignment chain pipeline [make_lastz_chains](https://github.com/hillerlab/make_lastz_chains)
> - **Scaffold names**: Input genome is renamed. See []() for the specific formatting rules.
> - **Inputs accepted**: `.fasta`, `.2bit`, or `.gz`.
> - **Custom library**: Repeat library can be provided under `repeat_library` in `params.json`.
> - **Container image**: We offer a pre-built container image for the whole pipeline as well as individual modules. By default the pipeline runs with [ghcr.io/hillerlab/softmask:latest](https://github.com/hillerlab/containers/pkgs/container/softmask). Additional images can be found at [containers](https://github.com/hillerlab/containers) and nextflow modules at [core](https://github.com/hillerlab/core).
> - **Input format**: Input assemblies may contain Ns and other ambiguous bases. Ambiguous-only chunks are preserved unchanged and produce no repeat annotations. Sequence previously replaced by Ns cannot be recovered; users requiring a completely re-derived soft-masked assembly should provide the original unmasked genome where available.

---

## Usage

> [!NOTE]
> Requirements: Nextflow ≥ 25.04.6, Docker or Apptainer, Java.

```bash
git clone https://github.com/hillerlab/softmask.git
cd softmask
```

Edit `params.json` (set `genome`, `assembly_prefix`), then:

```bash
# Docker
nextflow run main.nf -params-file params.json -profile docker

# Apptainer / Singularity
nextflow run main.nf -params-file params.json -profile apptainer
```

Smoke test:
```bash
nextflow run main.nf -profile test,apptainer
```

## Tests

Run the fast Python unit tests from the repository root:

```bash
python -m unittest discover -s src/tests -p 'test_*.py'
```

The end-to-end ambiguity test requires Nextflow 25.04.6 or newer and Docker.
It builds the pipeline image as `softmask:test` when that image is not already
available, then runs the pipeline once on a synthetic genome containing normal,
all-N, mixed, and IUPAC-ambiguity chunks:

```bash
SOFTMASK_RUN_INTEGRATION=1 python src/tests/test_pipeline_ambiguity.py
```

> [!NOTE]
> You can also specify these options directly in `params.json`.

A helper sh script is provided to run the pipeline on a SLURM cluster. See details below.

<details>
<summary>Click to expand</summary>


Edit the path variables at the top of `assets/hpc/do_softmask.sh` (cache dir, container image, manifest path), then submit:

```bash
sbatch --array=1-<N> do_softmask.sh
```

Each array task spawns one Nextflow head job that submits all compute as child SLURM jobs.

REPEATMASKER run as SLURM job arrays. Partition routing, array sizes, and resource tiers are documented inline in `nextflow.config` — edit there to match your cluster.

</details>

---

## Output

```
results/
├── 01_renamed/      *fasta
├── 02_database/     {assembly}.* [ repeatmodeler database ]
├── 03_model/        *.fa 
├── 04_mask/         *.masked/*out/*tbl
├── 05_final/        *.{2bit,fasta,gz}/*.repeats.tsv
└── pipeline_info/    timeline, trace, DAG
```

---

## Where to edit

| File | What |
|------|------|
| `params.json` | Genome paths, alignment settings, checkpoints — per run |
| `nextflow.config` | Compute resources, profiles, container, SLURM — rarely |
