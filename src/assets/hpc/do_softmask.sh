#!/bin/bash
# ─────────────────────────────────────────────────────────────────────────────
# SLURM submission script for softmask
#
# Each array task runs one soft-masking job. Nextflow itself is the
# "main job" — it submits all compute work as child SLURM jobs and only needs
# a small memory footprint.
#
# MANIFEST FILE FORMAT (species_list)
# ─────────────────────────────────────
# A tab-separated file with one pair per line, no header:
#
#   <genome_path>  <assembly_prefix>
#
# Example:
#   /data/genomes/hg38.2bit           Human
#   /data/genomes/mm39.fa             Mouse
#   /data/genomes/rn7.fasta.gz        Rat
#
# Genome files can be FASTA (.fa / .fasta / .fasta.gz / fa.gz) or 2bit (.2bit).
# Paths must be absolute. 
#
# USAGE
# ─────
# Edit the four path variables below, then submit with:
#   sbatch --array=1-<N> do_softmask.sh
# where <N> is the number of lines in your manifest file.
# ─────────────────────────────────────────────────────────────────────────────

#SBATCH --job-name=SOFTMASK
#SBATCH --array=1-10        # set upper bound to number of lines in species_list
#SBATCH -t 2-0
#SBATCH --output=/path/to/logs/%A.%a.out  # MODIFY THIS!
#SBATCH --error=/path/to/logs/%A.%a.err   # MODIFY THIS!
#SBATCH --mem=20G          # memory for the Nextflow process itself (not compute jobs)
#SBATCH -p long            # partition name
#SBATCH -q long            # queue name

# ── Load required modules (adjust to your cluster's module system) ────────────
module load nextflow
module load openjdk

# ── Environment ───────────────────────────────────────────────────────────────
export SLURM_SKIP_EPILOG=1

# Directory where Apptainer caches pulled container images
export NXF_APPTAINER_CACHEDIR=/scratch/$USER/softmask/apptainer

# Optional: pre-build a named SIF to avoid the auto-derived cache filename.
# Build once with:
#   apptainer build $NXF_APPTAINER_CACHEDIR/softmask.sif \
#       ghcr.io/hillerlab/softmask:latest
# Then uncomment:
# export NXF_CONTAINER_IMAGE=$NXF_APPTAINER_CACHEDIR/softmask.sif

# Give Nextflow's JVM enough heap for large runs (thousands of jobs)
export NXF_OPTS="-Xms4g -Xmx16g"

# ── Paths — edit these ────────────────────────────────────────────────────────
species_list="/path/to/manifest.tsv"   # tab-separated manifest (see format above)
working_dir="/path/to/output"          # one subdirectory per genome will be created here
pipeline_dir="/path/to/softmask"       # cloned pipeline repo

# ── Parse manifest line for this array task ───────────────────────────────────
row=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$species_list")
genome_path=$(echo "$row" | cut -f1)
assembly_prefix=$(echo "$row"   | cut -f2)


if [[ -z "$genome_path" || -z "$assembly_prefix" ]]; then
    echo "ERROR: could not parse line ${SLURM_ARRAY_TASK_ID} of ${species_list}" >&2
    exit 1
fi

# ── Per-pair working directory ─────────────────────────────────────────────────
pair_dir="${working_dir}/${genome_path}_${query_name}_mask"
mkdir -p "${pair_dir}/logs"

# ── Write params.json for this pair ───────────────────────────────────────────
# Scientific parameters go here; infrastructure stays in nextflow.config.
cat > "${pair_dir}/params.json" <<EOF
{
    "query_genome":  "${query_fa}",
    "genome": "${genome_path}",
    "assembly_prefix": "${assembly_prefix}",
    "repeat_library": null,
    "include_ltr": true,
    "engine" : "ncbi",
    "chunks": 100000,
    "output_format": "2bit",
    "outdir": "${pair_dir}/results",
    "use_container": true
}
EOF

# cd into pair_dir so each run's .nextflow.log is saved there
cd "$pair_dir"

nextflow run "${pipeline_dir}/main.nf" \
    -params-file "${pair_dir}/params.json" \
    -profile     apptainer,slurm \
    -w           "${pair_dir}/work"
