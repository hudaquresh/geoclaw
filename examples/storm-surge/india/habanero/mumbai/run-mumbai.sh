#!/bin/sh
#
# Simple "Hello World" submit script for Slurm.
#
# Replace <ACCOUNT> with your account name before submitting.
#
#SBATCH --account=apam      # The account name for the job.
#SBATCH --job-name=IkeRun    # The job name.
#SBATCH -c 1                  # The number of cpu cores to use.
#SBATCH --time=20:00              # The time the job will take to run.
#SBATCH --mem-per-cpu=1gb        # The memory the job will use per cpu core.
#SBATCH --nodes=4        # The memory the job will use per cpu core.

make all  

# End of script
