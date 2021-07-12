#!/bin/bash
#
#SBATCH -p medium
#SBATCH --mem=1G
#SBATCH --job-name=sim1
#SBATCH -o /mnt/shared/projects/nhm/goswamilab/rnf_methods/reports/sim1_%A.out
#SBATCH -e /mnt/shared/projects/nhm/goswamilab/rnf_methods/reports/sim1_%A.err

conda activate r_r
Rscript --no-save ~/projects/nhm/goswamilab/rnf_methods/script_200_30.R
