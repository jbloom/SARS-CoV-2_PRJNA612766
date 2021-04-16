#!/bin/bash

# https://vaneyckt.io/posts/safer_bash_scripts_with_set_euxo_pipefail/
set -euo pipefail

printf "Running snakemake"

slurm_scratch_dir="results/_scratch/slurm_logs/"
mkdir -p $slurm_scratch_dir

snakemake \
    -j 100 \
    --cluster-config cluster.yaml \
    --cluster "sbatch -c {cluster.cpus} -t {cluster.time} -J {cluster.name} --mem={cluster.memory} -o $slurm_scratch_dir/slurm-%j.out -e $slurm_scratch_dir/slurm-%j.out" \
    --latency-wait 60 \
    --use-conda \
    --rerun-incomplete

printf "Script complete."
