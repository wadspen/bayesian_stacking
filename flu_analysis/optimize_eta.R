library(dplyr)
library(lubridate)
library(stringr)
library(cmdstanr)
library(evalcast)
library(tidyr)
source("../simulations/stack_functions.R")
library(parallel)
library(doParallel)
library(doMC)
n.cores <- detectCores()
#n.cores <- 1
my.cluster <- makeCluster(n.cores, type = "PSOCK")
doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()
foreach::getDoParWorkers()
registerDoMC(cores = n.cores)
print("gets here")
mod <- cmdstan_model(stan_file = '../stan_models/emp_mix_crps_time_weight.stan')


mod_loc <- "../../FluSight-forecast-hub/model-output/"
# mod_loc <- "../../forecast-hub/FluSight-forecast-hub/model-output/" #local machine
models <- list.files(mod_loc)
models <- models[models != "README.md"]
sub_dates <- substr(list.files(paste0(mod_loc, "FluSight-baseline")), 1, 10)
sub_dates <- sub_dates[sub_dates < "2024-07-01"]
sub_dates <- sub_dates[-length(sub_dates)]
horizons <- -1:3
get_loc_file <- list.files(paste0(mod_loc, "FluSight-baseline/"))[4]
get_loc_forc <- read.csv(paste0(mod_loc, "FluSight-baseline/", get_loc_file))
locations <- unique(get_loc_forc$location)

#locations <- locations[sample(length(locations), 6)]

comp_forcs <- readRDS("comp_forcs.rds")
all_flu <- read.csv("../../FluSight-forecast-hub/target-data/target-hospital-admissions.csv")
# all_flu <- read.csv("../../forecast-hub/FluSight-forecast-hub/target-data/target-hospital-admissions.csv") #local machine
#etas <- seq(.5, 30, length.out = 30)
all_flu <- all_flu %>%
  mutate(reference_date = date, true_value = value) %>%
  dplyr::select(-date, -value)
etas <- seq(-1, 5, length.out = 20)
etas <- 1
#locations <- c("01", "16")
dat <- sub_dates[1]
horizon <- 0


stack_res <- foreach(loc = locations,
        .packages = c("cmdstanr", "stringr",
                      "lubridate", "dplyr")
        ,.errorhandling = "remove"
        ,.combine = rbind) %:% 
    foreach(d = 1:(length(sub_dates))
            ,.combine = rbind) %dopar% {
            # for (d in 1:(length(sub_dates))) {

            #look natural 45, 
            #look unatrual 46 not better after 10, 22 not better after 120ish
            # loc <- sample(locations, 1)
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
            etad[2] <- 1
            
              
              if (d == 1) {
                stack_crps[d] <- mix_mat_crps(weight[,d], all_mse[,d], 
                                              absdiff_arr[,,d])
	    	wts <- weight[,d]
              } else if (d == 2) {
                
                wts <- try(learning_rate(log(etad[d]), d-1, mse_mat = all_mse, 
                                         absdiff_arr = absdiff_arr, 
                                         mod = mod, power = 1, 
                                         alpha = 1, return_wts = TRUE))
                wts[wts < 0] <- 0
                wts <- wts/sum(wts)
                weight[,d] <- wts
                
                stack_crps[d] <- mix_mat_crps(weight[,d], all_mse[,d], 
                                              absdiff_arr[,,d])
              
              } else {
                  ev_grid <- c()
                  # for (i in 1:length(etas)) {
                  #   ev_grid[i] <- 
                  #     try(learning_rate(etas[i], d - 2, mse_mat = all_mse, 
                  #                       absdiff_arr = absdiff_arr, 
                  #                       mod = mod, power = 1, alpha = 1))
                  #   
                  # }
                  # etad[d] <- exp(etas[which.min(ev_grid)])
                  wts <- try(learning_rate(etas, d-1, mse_mat = all_mse, 
                                           absdiff_arr = absdiff_arr, 
                                           mod = mod, power = 1, 
					   alpha = 1, return_wts = TRUE))
                  wts[wts < 0] <- 0
                  wts <- wts/sum(wts)
                  weight[,d] <- wts
                  
                  stack_crps[d] <- mix_mat_crps(weight[,d], all_mse[,d], 
                                               absdiff_arr[,,d])
              } 
            
            
            
            comps <- comp_forcs %>% 
              filter(location == loc)
            
            comp_mods <- unique(comps$model)
            
            forcs <- data.frame()
            for (c in 1:length(comp_mods)) {
              comp_file <- paste0(mod_loc, comp_mods[c], "/", sub_dates[d], 
                                  "-", comp_mods[c], ".csv")
              
              forc <- read.csv(comp_file) %>%
                mutate(location = as.character(location)) %>%
                mutate(location = ifelse(nchar(location) < 2, 
                                         paste0("0", location), location)) %>%
                filter(horizon == horiz, target == "wk inc flu hosp",
                       location == loc) %>%
                left_join(all_flu, by = c("reference_date", "location"))
              
              forc$model <- comp_mods[c]
              forcs <- rbind(forcs, forc)
            }
            
            wts_df <- data.frame(model = comp_mods, wt = wts)
            forcs <- forcs %>%
              left_join(wts_df, by = "model")
            
            
            wis_wt <- forcs %>% 
              group_by(output_type_id, true_value) %>%
              summarise(value = sum(wt*value)) %>%
              ungroup() %>%
              summarise(wis = weighted_interval_score(as.numeric(output_type_id), 
                                                      value, unique(true_value)))
            
            
            
            wis_med <- forcs %>% 
              group_by(output_type_id, true_value) %>%
              summarise(value = median(value)) %>%
              ungroup() %>%
              summarise(wis = weighted_interval_score(as.numeric(output_type_id), 
                                                      value, unique(true_value)))
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
                ress <- data.frame(location = loc, week = d, 
                           forecast_date = sub_dates[d],
                           stack_crps = stack_crps[d], eq_crps = mean_crps[d],
                           stack_wis = wis_wt, med_wis = wis_med
			   eta = etad[d])
                
        

	    	ress
	        #wtsdf <- data.frame(t(as.matrix(wts)))
		#colnames(wtsdf) <- paste0("wt", 1:length(wts))

		#cbind(ress, wtsdf)
                
            }
            
           saveRDS(stack_res, "sgp_eq.rds") 
