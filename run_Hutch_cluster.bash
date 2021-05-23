#!/bin/bash

# https://vaneyckt.io/posts/safer_bash_scripts_with_set_euxo_pipefail/
set -euo pipefail

slurm_scratch_dir="results/_scratch/slurm_logs/"
mkdir -p $slurm_scratch_dir

printf "Running snakemake...\n"

snakemake \
    -j 100 \
    --cluster-config cluster.yaml \
    --cluster "sbatch -c {cluster.cpus} -t {cluster.time} -J {cluster.name} --mem={cluster.memory} -o $slurm_scratch_dir/slurm-%j.out -e $slurm_scratch_dir/slurm-%j.out" \
    --latency-wait 60 \
    --use-conda \
    --rerun-incomplete \
    --keep-going

printf "Run of snakemake complete.\n"
