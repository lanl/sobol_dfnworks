#!/bin/bash -l

#SBATCH -J calc_sobol
#SBATCH -p general
#SBATCH -c 1
#SBATCH -t 10:00:00
#SBATCH --output=SLURMOUT/calc_sobol_a_%A-%a.out

# Load software and activate environment (I loaded the modules into this env already):
module add R/4.1.2
R CMD BATCH "--no-save" sobol_on_seq_design_runs.R ./logfiles/calc_sobol_$SLURM_ARRAY_TASK_ID.Rout

