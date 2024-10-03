library(dplyr)
library(lubridate)
library(stringr)
library(cmdstanr)
library(tidyr)
library(ggplot2)
source("../simulations/stack_functions.R")
mod <- cmdstan_model(stan_file = '../stan_models/emp_mix_crps_time_weight.stan')


mod_loc <- "../../../forecast-hub/FluSight-forecast-hub/model-output/" #local machine
models <- list.files(mod_loc)
models <- models[models != "README.md"]
sub_dates <- substr(list.files(paste0(mod_loc, "FluSight-baseline")), 1, 10)
sub_dates <- sub_dates[-length(sub_dates)]
horizons <- -1:3
get_loc_file <- list.files(paste0(mod_loc, "FluSight-baseline/"))[4]
get_loc_forc <- read.csv(paste0(mod_loc, "FluSight-baseline/", get_loc_file))
locations <- unique(get_loc_forc$location)

#locations <- locations[sample(length(locations), 6)]

comp_forcs <- readRDS("../comp_forcs.rds")
all_flu <- read.csv("../../../forecast-hub/FluSight-forecast-hub/target-data/target-hospital-admissions.csv") #local machine
#etas <- seq(.5, 30, length.out = 30)
etas <- seq(-1, 5, length.out = 20)
#locations <- c("01", "16")
dat <- sub_dates[1]
horizon <- 0

fits_list <- list()
loc <- "US"
# ds <- c(6, 12, 18, 24)
ds <- 1:length(sub_dates)

for (k in 1:length(ds)) {
            d <- ds[k]
            comp_file <- paste0("../crps_comps2/crps_comps_", loc, ".rds")	
            #if (file.exists(comp_file) == FALSE) {next}	
            crps_comps <- readRDS(comp_file)
            
            absdiff_h_arr <- crps_comps[[1]]
            all_h_mse <- crps_comps[[2]]
            h_date <- crps_comps[[3]]
            
            h_date <- h_date %>%
              mutate(forecast_date = date(reference_date) + (horizon - 1)*7)
            
            nmods <- nrow(all_h_mse)
            wt <- rep(1/nmods, nmods)		
            weight <- matrix(NA, nrow = nmods, ncol = length(sub_dates))
            weight[,1] <- rep(1/nrow(weight), nrow(weight))
            # weight[,2] <- weight[,1]
            all_mse <- all_h_mse[,which(h_date$horizon == 1)]
            absdiff_arr <- absdiff_h_arr[,,which(h_date$horizon == 1)]	
            etad <- c()
            mean_crps <- mix_mat_crps(wt, all_mse, absdiff_arr)
            stack_crps <- c()
            
            etad[1] <- 1
            
            
            if (d == 1) {
              stack_crps[d] <- mix_mat_crps(weight[,d], all_mse[,d], 
                                            absdiff_arr[,,d])
              wts <- weight[,d]
            } else {
              ev_grid <- c()
              for (i in 1:length(etas)) {
                ev_grid[i] <- 
                  try(learning_rate(etas[i], d - 1, mse_mat = all_mse, 
                                    absdiff_arr = absdiff_arr, 
                                    mod = mod, power = 1))
                
              }
              etad[d] <- exp(etas[which.min(ev_grid)])
              
              mae <- all_mse[,1:d]
              absdiff = aperm(absdiff_arr[,,1:d])
              stan_dat <- list(
                T = d,
                num_comp = nrow(all_mse),
                eta = exp(log(etad[d])),
                alpha = rep(1, nrow(all_mse)),
                mae = mae,
                absdiff = absdiff
                , tweight = .98
                , power = 1)
              
              fits_list[[k]] <- mod$sample(data = stan_dat,
                                chains = 4, 
                                parallel_chains = 2,
                                iter_warmup = 10000,
                                iter_sampling = 50000)
              
            
              
              
              
            } 
}

saveRDS(fits_list, "us_all_season_fits")

