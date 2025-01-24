library(dplyr)
library(lubridate)
library(stringr)
source("./simulations/stack_functions.R")
library(cmdstanr)
library(dplyr)
library(evmix)
library(scoringRules)
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


mod_loc <- "../FluSight-forecast-hub/model-output/"
models <- list.files(mod_loc)
models <- models[models != "README.md"]
sub_dates <- substr(list.files(paste0(mod_loc, "FluSight-baseline")), 1, 10)
sub_dates <- sub_dates[-length(sub_dates)]
horizons <- 0:3
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
        ,.errorhandling = "remove"
        ) %dopar% {

#for (loc in c("01")) {
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
	
	nmods <- length(mods)
	
	all_logs <- matrix(NA, nrow = nmods, ncol = length(sub_dates)*4)
	all_crps <- matrix(NA, nrow = nmods, ncol = length(sub_dates)*4)	
	r <- 1
	h_dates <- data.frame()
	for (d in 1:length(sub_dates)) {

		dat <- sub_dates[d]
		for (horizon in 1:4) {
	loc_flu <- all_flu %>%
		filter(location == loc, date == date(dat) + (horizon - 1)*7)
        mse <- c()	
	mc_dist <- c()
	logs <- c()
	crps <- c()
	
		for (i in 1:length(mods)) {
		
			#for (j in i:length(mods)) {
			
				if (length(loc_flu$value) < 1) {logs[i] <- NA; crps[i] <- NA; next}	
				modi_name <- paste0(draws_loc, mods[i], "/draws/", dat, "-",
							    loc, "-", horizon - 1 , "-", mods[i], ".rds")
				modi <- readRDS(modi_name)
				logs[i] <- logs_sample(log(loc_flu$value + 1), modi$dist_samp)	
		 		mse[i] <- mean(abs(modi$dist_samp - log(loc_flu$value + 1)))
			        mc_dist[i] <- mean(abs(sample(modi$dist_samp, replace = FALSE) -
				        	 sample(modi$dist_samp, replace = FALSE)))

				crps[i] <- mse[i] - .5*mc_dist[i]
	
		}		
			h_date <- data.frame(reference_date = sub_dates[d], horizon = horizon)
			h_dates <- rbind(h_dates, h_date)
			
		
			#print(dim(absdiff_arr))
			all_logs[,r] <- logs
			all_crps[,r] <- crps
			r <- r + 1
		
			

		}
		print(d)
	}
	
	saveRDS(list(dates = h_dates, logs = all_logs, crps = all_crps), 
		paste0("./bma_bps_stat/loc", loc, ".rds"))
}
