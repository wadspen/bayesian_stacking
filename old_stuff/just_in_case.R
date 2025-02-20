library(dplyr)
library(lubridate)
library(stringr)
source("./stack_functions.R")
library(cmdstanr)
library(dplyr)
library(evmix)
library(parallel)
library(doParallel)
library(doMC)
n.cores <- detectCores()
my.cluster <- makeCluster(n.cores, type = "PSOCK")
doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()
foreach::getDoParWorkers()
registerDoMC(cores = n.cores)


mod_loc <- "../FluSight-forecast-hub/model-output/"
models <- list.files(mod_loc)
models <- models[models != "README.md"]
sub_dates <- substr(list.files(paste0(mod_loc, "FluSight-baseline")), 1, 10)
horizons <- -1:3
get_loc_file <- list.files(paste0(mod_loc, "FluSight-baseline/"))[4]
get_loc_forc <- read.csv(paste0(mod_loc, "FluSight-baseline/", get_loc_file))
locations <- unique(get_loc_forc$location)


comp_forcs <- readRDS("comp_forcs.rds")

keep_comp_forcs <- comp_forcs %>%
	group_by(model) %>%
	summarise(n = n()) %>%
	filter(n == 53)
all_flu <- read.csv("../FluSight-forecast-hub/target-data/target-hospital-admissions.csv")


#loc <- "35"  #32, 33, 35
dat <- sub_dates[1]
horizon <- 0
mods <- keep_comp_forcs$model
nmods <- length(mods)
foreach(loc = locations,
        .packages = c("stringr", "lubridate", "dplyr")
        #,.errorhandling = "remove"
        ) %dopar% {

all_mean_crps <- list()
all_stack_crps <- list()
#for (l in 2:length(locations)) {
	#loc <- locations[l]
	draws_loc <- "../quantile_modeling/quantile_fitting/model-fits/"
	
	print(loc)	
	loc_forcs <- comp_forcs %>%
		filter(location == loc)
	
	mods <- unique(loc_forcs$model); print(mods)
	if (loc == "10") {mods <- mods[mods != "UMass-trends_ensemble"]}
	if (loc == "02") {mods <- mods[mods != "MIGHTE-Nsemble"]}
	if (loc == "15") {mods <- mods[mods != "CEPH-Rtrend_fluH"]}
	if (loc == "23") {mods <- mods[mods != "SigSci-TSENS"]}
	if (loc == "32") {mods <- mods[mods != "MIGHTE-Nsemble"]}
	if (loc == "33") {mods <- mods[mods != "MIGHTE-Nsemble"]}
	if (loc == "35") {mods <- mods[mods != "MIGHTE-Nsemble"]}
	
	nnmods <- length(mods)
	wt <- rep(1/nmods, nmods)	
	absdiff_arr <- array(NA, dim = c(nmods, nmods, length(sub_dates)*4))
	all_mse <- matrix(NA, nrow = nmods, ncol = length(sub_dates)*4)
	r <- 1
	h_dates <- data.frame()
	for (d in 1:length(sub_dates)) {

		dat <- sub_dates[d]
		for (horizon in 1:4) {
	loc_flu <- all_flu %>%
		filter(location == loc, date == date(dat) + (horizon - 1)*7)
	
	mc_dist <- matrix(NA, nrow = length(mods), ncol = length(mods))
	mse <- c()
	
		for (i in 1:length(mods)) {
		
			#for (j in i:length(mods)) {
			
					
				modi_name <- paste0(draws_loc, mods[i], "/draws/", dat, "-",
							    loc, "-", horizon - 1 , "-", mods[i], ".rds")
				modi <- readRDS(modi_name)
			
		 		mse[i] <- mean(abs(modi$dist_samp - log(loc_flu$value + 1)))	
				for (j in i:length(mods)) {
					modj_name <- paste0(draws_loc, mods[j], "/draws/", dat, "-",
							    loc, "-", horizon - 1, "-", mods[j], ".rds")
					modj <- readRDS(modj_name)
		
					mc_dist[i,j] <- mean(abs(sample(modi$dist_samp, replace = FALSE) -
				        	 sample(modj$dist_samp, replace = FALSE)))
		
					mc_dist[j,i] <- mc_dist[i,j]		
				}
		}		
			h_date <- data.frame(reference_date = sub_dates[d], horizon = horizon)
			h_dates <- rbind(h_dates, h_date)
			
		

			absdiff_arr[,,r] <- mc_dist
			all_mse[,r] <- mse
			r <- r + 1
		
			

		}
		print(d)
	}
	save_name <- paste0("./crps_comps/crps_comps_", loc, ".rds")
	saveRDS(list(absdiff = absdiff_arr, mse = all_mse, h_date = h_dates), save_name)
	
	

}
