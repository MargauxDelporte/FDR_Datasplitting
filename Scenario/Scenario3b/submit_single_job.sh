#!/bin/bash -l
#SBATCH --partition=scu-cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=50
#SBATCH --job-name=scenario3
#SBATCH --time=48:00:00
#SBATCH --mem=200G
#SBATCH --output=scenario3_%j.out
#SBATCH --error=scenario3_%j.err
#SBATCH --mail-user=mde4023@med.cornell.edu
#SBATCH --mail-type=END

# Initialize conda (needed for batch jobs)
source /home/mde4023/miniconda3/etc/profile.d/conda.sh

# Activate R environment
conda activate r_env

# Set cores for parallel processing
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Navigate to working directory
cd /home/mde4023/FDR_Datasplitting/Scenario/Scenario3b

# Run the R script
R CMD BATCH --no-save Run_scenario3.R scenario3.Rout

# Deactivate conda
conda deactivate

echo "Job completed at $(date)"
