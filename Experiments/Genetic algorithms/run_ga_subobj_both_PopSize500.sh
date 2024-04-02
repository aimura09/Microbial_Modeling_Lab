#!/bin/bash
#SBATCH -c 4  # number of cores
#SBATCH -t 0-22:00:00   # time in d-hh:mm:ss
#SBATCH -p general      # Partition (default for jobs < 4 hours) 
#SBATCH -q public       # QOS (default for main cpu partitions)
#SBATCH -o ./log/slurm-%j.out # File to save job's STDOUT to (%j = JobId)
#SBATCH -e ./log/slurm-%j.err # File to save job's STDERR to (%j = JobId)
#SBATCH --mail-type=ALL # Send a notification when a job starts, stops, or fails
#SBATCH --mail-user=%u@asu.edu # Send-to address (%u expands to username)

# Always purge modules to ensure a consistent environment
module purge 
# Load version of Matlab
module load matlab
# Run script
matlab -nodisplay -nosplash -r "addpath('cobratoolbox/'); initCobraToolbox; GA_subtraction_object_cmdarg sub_object_both_pop500 1e-4 0.85 500 12 1000 both; exit"
