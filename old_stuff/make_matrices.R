library(dplyr)
library(lubridate)
library(stringr)
source("./simulations/stack_functions.R")


mod_loc <- "../FluSight-forecast-hub/model-output/"
models <- list.files(mod_loc)
models <- models[models != "README.md"]
sub_dates <- substr(list.files(paste0(mod_loc, "FluSight-baseline")), 1, 10)
horizons <- 0:3
get_loc_file <- list.files(paste0(mod_loc, "FluSight-baseline/"))[4]
get_loc_forc <- read.csv(paste0(mod_loc, "FluSight-baseline/", get_loc_file))
locations <- unique(get_loc_forc$location)



keep_models <- c()
for (i in 1:length(models)) {
	forcs <- list.files(paste0(mod_loc, models[i]))
	forc_dates <- substr(forcs, 1, 10)
	miss_date <- setdiff(sub_dates, forc_dates)
	if (length(miss_date) == 0) {keep_models <- c(keep_models, models[i])}
	
}

#note that Iowa is location 19

all_forcs <- data.frame()
for (i in 1:length(keep_models)) {
	mod <- keep_models[i]
	files <- paste0(mod_loc, mod, "/", list.files(paste0(mod_loc, mod), pattern = "\\.csv"))
	tables <- lapply(files, read.csv)
	forcs <- do.call(rbind, tables)
	
	#for (j in 1:length(sub_dates)) {
	#	mod <- keep_models[i]
	#	dat <- sub_dates[j]
		#forcs <- read.csv(paste0(mod_loc, mod, "/", dat, "-", mod, ".csv"))
		forcs <- forcs %>%
			mutate(location = as.character(location), model = mod) %>%
			mutate(location = ifelse(nchar(location) < 2, paste0(0, location),
			       location)) %>%
			filter(horizon != -1, output_type == "quantile")
		
		all_forcs <- rbind(all_forcs, forcs); print(i)
	#	print(i); print(j)
	#}
}

all_forcs <- all_forcs %>%
	filter(str_detect(model, "FluSight") == FALSE | model == "FluSight-baseline")


comp_forcs <- all_forcs %>%
	filter(reference_date %in% sub_dates, output_type_id == .025, horizon %in% 0:3) %>%
	group_by(model, location) %>% 
	summarise(n = n()) %>%
        ungroup() %>%
	mutate(totn = max(n)) %>%
	filter(n == totn)	

#saveRDS(comp_forcs, "comp_forcs.rds")
#comp_forcs <- readRDS("comp_forcs.rds")
saveRDS(comp_forcs, "comp_forcs_new.rds")

all_flu <- read.csv("../FluSight-forecast-hub/target-data/target-hospital-admissions.csv")


loc <- "12"
dat <- sub_dates[1]
horizon <- 0


all_mean_crps <- list()
all_stack_crps <- list()
for (l in 2:length(locations)) {
	loc <- locations[l]
	draws_loc <- "../quantile_modeling/quantile_fitting/model-fits/"
	
	
	loc_forcs <- comp_forcs %>%
		filter(location == loc)
	
	mods <- unique(loc_forcs$model); print(mods)
	if (locations[l] == "10") {mods <- mods[mods != "UMass-trends_ensemble"]}
	nmods <- length(mods)
	wt <- rep(1/nmods, nmods)	
	absdiff_arr <- array(NA, dim = c(nmods, nmods, length(sub_dates)))
	all_mse <- matrix(NA, nrow = nmods, ncol = length(sub_dates))
	for (d in 1:length(sub_dates)) {
		dat <- sub_dates[d]
	loc_flu <- all_flu %>%
		filter(location == loc, date == dat)
	
	mc_dist <- matrix(NA, nrow = length(mods), ncol = length(mods))
	mse <- c()
		for (i in 1:length(mods)) {
		
			#for (j in i:length(mods)) {
		
			modi_name <- paste0(draws_loc, mods[i], "/draws/", dat, "-",
						    loc, "-", horizon, "-", mods[i], ".rds")
			modi <- readRDS(modi_name)
			
		 	mse[i] <- mean(abs(modi$dist_samp - log(loc_flu$value + 1)))	
			for (j in i:length(mods)) {
				modj_name <- paste0(draws_loc, mods[j], "/draws/", dat, "-",
						    loc, "-", horizon, "-", mods[j], ".rds")
				modj <- readRDS(modj_name)
		
				mc_dist[i,j] <- mean(abs(sample(modi$dist_samp, replace = FALSE) -
				         sample(modj$dist_samp, replace = FALSE)))
		
				mc_dist[j,i] <- mc_dist[i,j]		
			}
		}
		
		absdiff_arr[,,d] <- mc_dist
		all_mse[,d] <- mse
		print(d)
	}
	
	
	
	#vMv <- function(mat, wt) {
	#		return(wt%*%mat%*%wt)
	#}
	#
	#mix_mat_crps <- function(wt, mse_mat, absdiff_arr) {
	#			#wt <- exp(wt)/sum(exp(wt))
	#			one <- apply(mse_mat, MARGIN = 2, FUN = function(x) {sum(wt*x)})
	#			two <- apply(absdiff_arr, MARGIN = 3, FUN = vMv, wt = wt)
	#
	#			crpss <- one - (1/2)*two
	#			return(crpss)
	#}
	#
	#sum_crps <- function(wt, mse_mat, absdiff_arr) {
	#	wt <- exp(wt)/sum(exp(wt))
	#	crpss <- mix_mat_crps(wt, mse_mat, absdiff_arr)
	#	return(sum(crpss))
	#}
	
	
	
#	weights <- matrix(NA, nrow = nmods, ncol = length(sub_dates))
#	weights[,1] <- rep(1/nrow(weights), nrow(weights))
#	
#	
#	for (d in 1:(length(sub_dates) - 1)) {
#		
#		wts <- optim(log(weights[,d] + 1), sum_crps, mse_mat = all_mse[,1:d], 
#			     absdiff_arr = absdiff_arr[,,1:d])$par
#		weights[,d+1] <- exp(wts)/sum(exp(wts))
#	}
#
#
#
#	mean_crps <- mix_mat_crps(wt, all_mse, absdiff_arr)
#	stack_crps <- c()
#	for (i in 1:ncol(all_mse)) {
#		stack_crps[i] <- mix_mat_crps(weights[,i], all_mse[,i], absdiff_arr[,,i])
#	}
#	print(mean(stack_crps < mean_crps, na.rm = TRUE))
#	print(mean(mean_crps/stack_crps, na.rm = TRUE))
#	all_mean_crps[[l]] <- mean_crps
#	all_stack_crps[[l]] <- stack_crps
#
#}
