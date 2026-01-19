#!/bin/bash -l
#SBATCH --partition=scu-cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --job-name=case
#SBATCH --time=48:00:00
#SBATCH --mem=200G
#SBATCH --output=case_%j.out
#SBATCH --error=case_%j.err
#SBATCH --mail-user=mde4023@med.cornell.edu
#SBATCH --mail-type=END

# Initialize conda (needed for batch jobs)
source /home/mde4023/miniconda3/etc/profile.d/conda.sh

# Activate R environment
conda activate r_env

# Set cores for parallel processing
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Navigate to working directory
cd /home/mde4023/FDR_Datasplitting/Case

# Run the R script
R CMD BATCH --no-save Finetuning_Qin2014_xgb_selectionrates.R case.Rout

# Deactivate conda
conda deactivate

echo "Job completed at $(date)"
EOF