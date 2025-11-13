#!/bin/bash

#SBATCH --time=3-01:30:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --mem=0
#SBATCH --exclusive
#SBATCH --constraint=intel


# Old Module Load
#module load gcc
#module load r
#module load udunits
#module load r-rgdal
#module load proj
#module load r gsl

# modules for Rscript run
module purge
module load r/4.4.1
module load r-rgdal gsl udunits/2.2.28-et3j662

Rscript ./SIR_simulation.R "$1" 


