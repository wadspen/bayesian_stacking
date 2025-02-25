library(dplyr)
library(tidyr)
library(lubridate)
library(evalcast)
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

comp_forcs <- readRDS("../comp_forcs.rds")
#comp_forcs <- comp_forcs %>%
	#group_by(model) %>%
	#summarise(n = n()) %>%
	#filter(n == 53)
all_flu <- read.csv("../../FluSight-forecast-hub/target-data/target-hospital-admissions.csv")
# all_flu <- read.csv("../../forecast-hub/FluSight-forecast-hub/target-data/target-hospital-admissions.csv") #local machine
#etas <- seq(.5, 30, length.out = 30)
all_flu <- all_flu %>%
	mutate(reference_date = date, true_value = value) %>%
	dplyr::select(-date, -value)
etas <- seq(-1, 5, length.out = 20)
#locations <- c("01", "16")
dat <- sub_dates[1]
horiz <- 0

loc <- "01"

comps <- comp_forcs %>% 
	filter(location == loc)

comp_mods <- unique(comps$model)
print(sub_dates)

forcs <- data.frame()
for (c in 1:length(comp_mods)) {
	comp_file <- paste0(mod_loc, comp_mods[c], "/", sub_dates[2], 
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

wts <- MCMCprecision::rdirichlet(1, a = rep(1, length(comp_mods)))
wts <- as.vector(wts)
wts_df <- data.frame(model = comp_mods, wt = wts)
forcs <- forcs %>%
	left_join(wts_df, by = "model")


test <- forcs %>% 
	group_by(output_type_id, true_value) %>%
	summarise(value = sum(wt*value)) %>%
	ungroup() %>%
	summarise(wis = weighted_interval_score(as.numeric(output_type_id), 
						value, unique(true_value)))

