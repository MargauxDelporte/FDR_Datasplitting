#!/bin/bash -l
#SBATCH --partition=scu-cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --job-name=scen1
#SBATCH --time=48:00:00
#SBATCH --mem=200G
#SBATCH --output=scenario1_%j.out
#SBATCH --error=scenario1_%j.err
#SBATCH --mail-user=mde4023@med.cornell.edu
#SBATCH --mail-type=END

# Initialize conda (needed for batch jobs)
source /home/mde4023/miniconda3/etc/profile.d/conda.sh

# Activate R environment
conda activate r_env

# Set cores for parallel processing
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Navigate to working directory
cd /home/mde4023/FDR_Datasplitting/Scenario/Scenario1

# Run the R script
R CMD BATCH --no-save run_scenario1.R scenario1.Rout

# Deactivate conda
conda deactivate

echo "Job completed at $(date)"
EOF