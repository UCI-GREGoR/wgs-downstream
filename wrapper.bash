#!/usr/bin/env bash
snakemake -j23 --profile ../sge-profile -p --rerun-incomplete --rerun-triggers mtime --use-conda --use-singularity --cluster-config config/cluster.yaml
