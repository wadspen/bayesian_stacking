#!/bin/bash

#SBATCH --time=1-16:04:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --mem=0
#SBATCH --exclusive
#SBATCH --constraint=intel

Rscript simple_iid.R 

