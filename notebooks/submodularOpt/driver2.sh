#!/bin/bash

#SBATCH --job-name=soss2015
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=jpic@umich.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --mem-per-cpu=50G
#SBATCH --time=48:00:00
#SBATCH --account=indikar0
#SBATCH --partition=standard
# SBATCH --array=

module load matlab
echo $1
echo $2
echo $3
matlab -nodisplay -r "addpath(genpath(pwd)); driver2($1, 100, $2, $3)"