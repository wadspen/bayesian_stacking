library(dplyr)
library(lubridate)
library(stringr)
source("./stack_functions.R")


mod_loc <- "../FluSight-forecast-hub/model-output/"
models <- list.files(mod_loc)
models <- models[models != "README.md"]
sub_dates <- substr(list.files(paste0(mod_loc, "FluSight-baseline")), 1, 10)
horizons <- -1:3
get_loc_file <- list.files(paste0(mod_loc, "FluSight-baseline/"))[4]
get_loc_forc <- read.csv(paste0(mod_loc, "FluSight-baseline/", get_loc_file))
locations <- unique(get_loc_forc$location)



comp_forcs <- readRDS("comp_forcs.rds")
all_flu <- read.csv("../FluSight-forecast-hub/target-data/target-hospital-admissions.csv")


loc <- "12"
dat <- sub_dates[1]
horizon <- 0

mean_scores <- c()
m <- 1
for (loc in locations) {
        comp_file <- paste0("./crps_comps/crps_comps_", loc, ".rds")	
	if (file.exists(comp_file) == FALSE) {next}	
	crps_comps <- readRDS(comp_file)
	
	absdiff_h_arr <- crps_comps[[1]]
	all_h_mse <- crps_comps[[2]]
	h_date <- crps_comps[[3]]

	h_date <- h_date %>%
		mutate(forecast_date = date(reference_date) + (horizon - 1)*7)

	nmods <- nrow(all_h_mse)
	wt <- rep(1/nmods, nmods)		
	weights <- matrix(NA, nrow = nmods, ncol = length(sub_dates))
	weights[,1] <- rep(1/nrow(weights), nrow(weights))
	
	all_mse <- all_h_mse[,which(h_date$horizon == 1)]
	absdiff_arr <- absdiff_h_arr[,,which(h_date$horizon == 1)]	
	for (d in 1:(length(sub_dates) - 1)) {
	#	pre_forcs <- which(h_date$forecast_date <= sub_dates[d])
	#	if (length(pre_forcs) == 1) {all_mse <- all_h_mse; absdiff_arr <- absdiff_h_arr}
	#	else {
	#		all_mse <- all_h_mse[,pre_forcs]
	#		absdiff_arr <- absdiff_h_arr[,,pre_forcs]
	#	}
		
		wts <- optim(log(weights[,d] + 1), sum_crps, mse_mat = all_mse[,1:d], 
			     absdiff_arr = absdiff_arr[,,1:d])$par
		weights[,d+1] <- exp(wts)/sum(exp(wts))
	}


	one_week <- which(h_date$horizon == 1)
        all_mse <- all_h_mse[,one_week]
 	absdiff_arr <- absdiff_h_arr[,,one_week]	
	mean_crps <- mix_mat_crps(wt, all_mse, absdiff_arr)
	stack_crps <- c()
	for (i in 1:ncol(all_mse)) {
		stack_crps[i] <- mix_mat_crps(weights[,i], all_mse[,i], absdiff_arr[,,i])
	}
	mean_scores[m] <- mean(mean_crps/stack_crps, na.rm = TRUE)
	#print(mean(stack_crps < mean_crps, na.rm = TRUE))
	print(mean_scores[m])
	m <- m + 1
	

}
