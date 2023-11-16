#!/bin/bash

#SBATCH --job-name=individual-gramians
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=jpic@umich.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=180G
#SBATCH --time=4:00:00
#SBATCH --account=indikar0
#SBATCH --partition=standard

# The application(s) to execute along with its input arguments and options:
python individualGramian.py
# module load matlab
# matlab -nodisplay -r "driverTest4('GL', $SLURM_ARRAY_TASK_ID,100)"
