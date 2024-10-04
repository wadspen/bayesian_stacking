library(dplyr)
library(tidyr)
library(stringr)
library(lubridate)
source("../simulations/stack_functions.R")



mod_loc <- "../../FluSight-forecast-hub/model-output/"
#mod_loc <- "../../../forecast-hub/FluSight-forecast-hub/model-output/" #local machine
models <- list.files(mod_loc)
models <- models[models != "README.md"]
sub_dates <- substr(list.files(paste0(mod_loc, "FluSight-baseline")), 1, 10)
horizons <- 0:3
get_loc_file <- list.files(paste0(mod_loc, "FluSight-baseline/"))[4]
get_loc_forc <- read.csv(paste0(mod_loc, "FluSight-baseline/", get_loc_file))
locations <- unique(get_loc_forc$location)

#locations <- locations[sample(length(locations), 6)]

comp_forcs <- readRDS("comp_forcs.rds")
all_flu <- read.csv("../../FluSight-forecast-hub/target-data/target-hospital-admissions.csv")
#all_flu <- read.csv("../../../forecast-hub/FluSight-forecast-hub/target-data/target-hospital-admissions.csv") #local machine

#

alpha <- .98

score_files <- list.files("./bma_bps_stat/", pattern = ".rds")
et <- .01
all_crps_wts <- data.frame()
for (i in 1:length(locations)) {
	print(locations[i])
	score_rds <- readRDS(paste0("./bma_bps_stat/loc", locations[i], ".rds"))
	score_date <- score_rds$date
	score_logs <- score_rds$logs
	score_crps <- score_rds$crps

	comp_file <- paste0("../crps_comps/crps_comps_", locations[i], ".rds")
	crps_comps <- readRDS(comp_file)
	absdiff_h_arr <- crps_comps[[1]]
	all_h_mse <- crps_comps[[2]]
	h_date <- crps_comps[[3]]

	h_date <- h_date %>%
		mutate(forecast_date = date(reference_date) + (horizon - 1)*7)
        
        one_week <- which(h_date$horizon == 1)
	all_mse <- all_h_mse[,one_week]
	absdiff_arr <- absdiff_h_arr[,,one_week]	
	
	logsh1 <- score_logs[,score_date$horizon == 1]
	crpsh1 <- score_crps[,score_date$horizon == 1]
	
	C <- nrow(crpsh1)
	bma_wt <- matrix(NA, nrow = C, ncol = ncol(all_mse))
	avs_wt <- matrix(NA, nrow = C, ncol = ncol(all_mse))
	bma_wt[,1] <- rep(1/C, C)
        avs_wt[,1] <- rep(1/C, C)	


	select_ets <- c()
	for (d in 1:(ncol(crpsh1) - 1)) {
		etas <- seq(.001, 6, length.out = 35)
    		min_eta <- c()
    		for (n in 1:(d - 1)) {
      
      			avcrpss <- c()
      			for(e in 1:length(etas)) {
        			et <- etas[e]
        			wavs <- c()
        			for (m in 1:C) {
          				wavs[m] <- (1/C)*exp(-et*sum(
            				alpha^(n:1 - 1)*crpsh1[m:1:n]))
        			}
        			wavs <- wavs/sum(wavs)
				all_crps <- c()
				for (k in 1:(n + 1)) {
					all_crps[k] <- mix_mat_crps(wavs, all_mse[,k], absdiff_arr[,,k])
				}
        			avcrpss[e] <- mean(all_crps)
      			}
      			min_eta[n] <- etas[which.min(avcrpss)]
    		}
    		avs_et <- mean(min_eta)
		select_ets[d + 1] <- avs_et
		wbma <- c()
		for (m in 1:C) {
			wbma[m] <- (1/C)*prod(-exp(logsh1[m,1:d]))
		}
		bma_wt[,d + 1] <- wbma/sum(wbma)
		wavs <- c()
		for (m in 1:C) {
			wavs[m] <- (1/C)*exp(-avs_et*sum(
		   	alpha^(d:1 - 1)*crpsh1[m,1:d]))
		}
		avs_wt[,d + 1] <- wavs/sum(wavs)
	}


	stack_bma <- c()
	stack_bps <- c()
	for (j in 1:ncol(all_mse)) {
		stack_bma[j] <- mix_mat_crps(bma_wt[,j], all_mse[,j], absdiff_arr[,,j])
		stack_bps[j] <- mix_mat_crps(avs_wt[,j], all_mse[,j], absdiff_arr[,,j])
	}

stack_bma <- stack_bma[!is.na(stack_bma)]
stack_bps <- stack_bps[!is.na(stack_bps)]
print(locations[i])
crps_eta_wts <- data.frame(location = locations[i], season_week = 1:(ncol(crpsh1)),
			   bma_crps = stack_bma, 
			   bps_crps = stack_bps, eta = select_ets) 

all_crps_wts <- rbind(all_crps_wts, crps_eta_wts)
}

saveRDS(all_crps_wts, "./bma_bps_crps.rds")
