#!/bin/bash

#SBATCH --time=1-16:04:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --mem=0
#SBATCH --exclusive
#SBATCH --constraint=intel

#Rscript get_ws_gibbs_optim.R 
#Rscript time_var_sim.R
#Rscript optimize_eta.R
Rscript make_matrices_all_weeks.R
