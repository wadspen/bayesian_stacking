library(dplyr)
library(lubridate)
library(stringr)
library(cmdstanr)
library(evalcast)
library(tidyr)
library(tibble)
source("../simulations/stack_functions.R")

print("gets here")
mod <- cmdstan_model(stan_file = '../stan_models/emp_mix_crps_time_weight.stan')


mod_loc <- "../../FluSight-forecast-hub/model-output/"
# mod_loc <- "../../forecast-hub/FluSight-forecast-hub/model-output/" #local machine
models <- list.files(mod_loc)
models <- models[models != "README.md"]
sub_dates <- substr(list.files(paste0(mod_loc, "FluSight-baseline")), 1, 10)
sub_dates <- sub_dates[sub_dates < "2024-07-01"]
# sub_dates <- sub_dates[-length(sub_dates)]
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
  dplyr::select(-date, -value) %>%
  mutate(true_value = log(true_value + 1))
etas <- seq(-1, 5, length.out = 20)
etas <- 1
#locations <- c("01", "16")
dat <- sub_dates[1]
horiz <- 0
loc <- "01"
d <- 2


#{

# loc <- "US"
loc <- "44" #so far 02 is the best
all_bounds <- data.frame()
unique(all_flu$location_name[all_flu$location == loc])
for (d in 2:(length(sub_dates))) {
            
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
            
            wts <- try(learning_rate(etas, d-1, mse_mat = all_mse, 
                                     absdiff_arr = absdiff_arr, 
                                     mod = mod, power = 1, 
                                     alpha = 1, return_wts = "draws"))
            
            bounds <- apply(wts, MARGIN = 2, FUN = quantile, 
                      probs = c(0.05, 0.95)) %>% 
                      t() %>% 
                      as.data.frame()
          
            colnames(bounds) <- c("lower", "upper")
            bounds$week <- sub_dates[d]
            bounds <- rownames_to_column(bounds, var = "component")
            
            all_bounds <- rbind(all_bounds, bounds)
            
}

saveRDS(all_bounds, "wt_bds_eg.rds")

library(MCMCprecision)
prior <- rdirichlet(50000, a = rep(1,10))
colnames(prior) <- paste0("omegas[", 1:10, "]")

priorq <- apply(prior, margin = 2, MARGIN = 2,
                FUN = quantile, probs = c(0.05, 0.95)) %>% 
  t() %>% 
  as.data.frame() %>% 
  mutate(week = "prior") %>% 
  rownames_to_column(var = "component") %>% 
  mutate(lower = `5%`, upper = `95%`) %>% 
  select(component, lower, upper, week)



rbind(all_bounds, priorq) %>% 
  pivot_longer(2:3, names_to = "bound") %>% 
  ggplot() +
  geom_line(aes(y = component, x = value), size = .7) + 
  facet_wrap(~week) +
  ylab("") +
  xlab("Weight") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 14),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank())

