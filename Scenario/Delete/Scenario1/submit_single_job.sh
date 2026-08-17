#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=51
#SBATCH --job-name=scen1
#SBATCH --time=48:00:00
#SBATCH --mem=200G
#SBATCH --output=scenario1_%j.out
#SBATCH --error=scenario1_%j.err

# Load R inside the batch job
module purge
module load r/4.4.0

# Show which R is being used
echo "R location:"
which R
R --version

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

cd /home/margaux_delporte_uri_edu/Scenario/Scenario1

echo "Starting R job at $(date)"

R CMD BATCH --no-save run_scenario1.R scenario1.Rout

echo "Job completed at $(date)"