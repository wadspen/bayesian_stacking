#!/bin/bash

#SBATCH --time=30:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --mem=0
#SBATCH --exclusive
#SBATCH --constraint=intel
module load gcc
module load r
module load udunits
module load r-rgdal
module load proj
module load r gsl



cd ../../FluSight-forecast-hub
git pull

cd ../gpo_flusight_24/gpo_flusight_24
ls



Rscript ./dyn_weight_forecasts.R #"$1" 
bash push_forecasts.sh
