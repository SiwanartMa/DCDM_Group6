#!/bin/bash -l
#SBATCH --job-name=ops-r
#SBATCH --partition=cpu

Rscript merge_files_batch.R
