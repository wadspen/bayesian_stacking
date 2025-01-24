library(dplyr)
library(lubridate)
library(stringr)
library(cmdstanr)
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
#mod_loc <- "../../../forecast-hub/FluSight-forecast-hub/model-output/" #local machine
models <- list.files(mod_loc)
models <- models[models != "README.md"]
sub_dates <- substr(list.files(paste0(mod_loc, "FluSight-baseline")), 1, 10)
horizons <- -1:3
get_loc_file <- list.files(paste0(mod_loc, "FluSight-baseline/"))[4]
get_loc_forc <- read.csv(paste0(mod_loc, "FluSight-baseline/", get_loc_file))
locations <- unique(get_loc_forc$location)

#locations <- locations[sample(length(locations), 6)]

comp_forcs <- readRDS("comp_forcs.rds")
all_flu <- read.csv("../../FluSight-forecast-hub/target-data/target-hospital-admissions.csv")
#all_flu <- read.csv("../../../forecast-hub/FluSight-forecast-hub/target-data/target-hospital-admissions.csv") #local machine

#locations <- "16"
loc <- "16"
dat <- sub_dates[1]
horizon <- 0
ms_et <- 3
mean_score <- c(); print(locations)
etas <- seq(.5, 20, length.out = 20)
mean_scores <- foreach(loc = locations,
                    .packages = c("cmdstanr", "stringr",
                                  "lubridate", "dplyr")
                    #,.errorhandling = "remove"
                    ,.combine = rbind) %dopar% {

        comp_file <- paste0("../crps_comps/crps_comps_", loc, ".rds")	
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
	all_mse <- all_h_mse[,which(h_date$horizon == 1)]
	absdiff_arr <- absdiff_h_arr[,,which(h_date$horizon == 1)]	
	for (d in 1:(length(sub_dates) - 1)) {
		
		if (d == 1) {
			mae <- matrix(all_mse[,1:d], nrow = nrow(all_mse))
			absdiff <- aperm(array(absdiff_arr[,,1:d], dim = c(nrow(all_mse), nrow(all_mse), 1)))
			
		} else {
			mae <- all_mse[,1:d]
			absdiff <- aperm(absdiff_arr[,,1:d])
		}
	#	pre_forcs <- which(h_date$forecast_date <= sub_dates[d])
	#	if (length(pre_forcs) == 1) {all_mse <- all_h_mse; absdiff_arr <- absdiff_h_arr}
	#	else {
	#		all_mse <- all_h_mse[,pre_forcs]
	#		absdiff_arr <- absdiff_h_arr[,,pre_forcs]
	#	}
  	  
	  for (e in 1:length(etas)) {
		wts <- learning_rate(etas[e], d - 1, mse_mat = mae[,1:(d-1)], absdiff_arr = absdiff_arr[,,1:(d-1)],
			      mod = mod, tweight = .98, power = 1, return_wts = TRUE)
	  }
  	  stan_dat <- list(
        	     T = d,
        	     num_comp = nrow(all_mse),
        	     eta = ms_et,
        	     alpha = rep(1, nrow(all_mse)),
        	     mae = mae,
        	     absdiff = absdiff,
		     power = 1
                     , tweight = .98)
	  
	   fit <- mod$sample(data = stan_dat,
	                     chains = 1,
	                     iter_warmup = 5000,
	                     iter_sampling = 5000,
	   		     init = list(list(omega = rep(1/nrow(all_mse), nrow(all_mse)))))
	  
	   #fit <- mod$variational(data = stan_dat)
	  
	   draws <- fit$draws(format = "df") %>% 
	     select(contains("omega"))
	   weight[,d+1] <- apply(draws, MARGIN = 2, FUN = mean)
		
	}

	weight_df <- weight %>% 
		t() %>%
		as.data.frame()
	colnames(weight_df) <- paste0("wt", 1:ncol(weight_df))
	one_week <- which(h_date$horizon == 1)
        all_mse <- all_h_mse[,one_week]
 	absdiff_arr <- absdiff_h_arr[,,one_week]	
	mean_crps <- mix_mat_crps(wt, all_mse, absdiff_arr)
	stack_crps <- c()
	for (i in 1:ncol(all_mse)) {
		stack_crps[i] <- mix_mat_crps(weight[,i], all_mse[,i], absdiff_arr[,,i])
	}
	
	flu_res <- cbind(data.frame(location = loc, eta = ms_et, stack_crps = stack_crps, 
			 eq_crps = mean_crps), weight_df)

	saveRDS(flu_res, paste0("./loc_scores/loc", loc, ".rds"))
	
	

}

#saveRDS(mean_scores, "mean_score_weight.rds")
