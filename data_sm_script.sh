#!/usr/bin/bash

snakemake --rerun-incomplete --keep-going --jobs 50 \
        --latency-wait 60 --wait-for-files \
        --cluster-config scripts/data_sm_slurm_config.json \
        --rerun-triggers mtime \
        -s Snakefile \
        --cluster "sbatch -p {cluster.queue} \
                        -t {cluster.time} \
                        --ntasks={cluster.tasks} \
                        --job-name={cluster.name} \
                        -o {cluster.output} \
                        --cpus-per-task {cluster.cores} \
                        -e {cluster.error} \
                        --nodes={cluster.nodes} \
                        --mem={cluster.memory}"
