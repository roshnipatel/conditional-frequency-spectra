#!/usr/bin/bash

snakemake --rerun-incomplete --keep-going --jobs 100 \
        --latency-wait 60 --wait-for-files \
        --cluster-config scripts/simulations_sm_slurm_config.json \
        -s Snakefile \
        --rerun-triggers mtime \
        --cluster "sbatch -p {cluster.queue} \
                        -t {cluster.time} \
                        --ntasks-per-node={cluster.tasks} \
                        --job-name={cluster.name} \
                        -o {cluster.output} \
                        -e {cluster.error} \
                        --nodes={cluster.nodes} \
                        --mem={cluster.memory}"
