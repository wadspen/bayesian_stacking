library(SimInf)
library(forecast)
library(tidyr)
library(ggplot2)
library(cmdstanr)
library(dplyr)
library(car)
library(stats)
#library(scoringutils)
library(scoringRules)
source("./stack_functions.R")
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

args <- commandArgs()
step <- as.numeric(args[6])

asg_mod <- cmdstan_model(stan_file = '../stan_models/asg.stan')
sir_mod <- cmdstan_model(stan_file = '../stan_models/sir.stan')
mod <- cmdstan_model(stan_file = '../stan_models/emp_mix_crps_time_weight.stan')
drawn <- 2000
warm <- 1000
reps <- 60

sir_res <- foreach(replicate = 1:reps,
                   .packages = c("cmdstanr", "stringr", "scoringRules",
                                 "lubridate", "dplyr", "tidyr")
                   ,.errorhandling = "remove"
                   ,.combine = rbind)  %dopar% {
                     N <- 800
                     I0 <- 12
                     days <- 300
                     ## Create an SIR model object.
                     model <- SIR(u0 = data.frame(S = N, I = I0, R = 0),
                                  tspan = 1:days,
                                  beta = rnorm(1, .06, .01),
                                  gamma = rnorm(1, 0.034, .009))
                     
                     ## Run the SIR model and plot the result.
                     # set.seed(22)
                     result <- run(model)
                     
                     
                     #simulate SIR data
                     I <- result@U[2,]
                     time <- result@tspan
                     wktime <- seq(5, days, by = 7)
                     wkI <- I[wktime]
                     # plot(wkI~wktime)
                     models <- c("naive", "arima", "asg", "sir")#, "mean")
                     scores <- data.frame()
                     absdiff_arr <- array(NA, dim = c(4,4,length(5:(length(wkI) - 1))))
                     anweek <- 0
                     all_mse <- matrix(NA, nrow = 4, ncol = length(5:(length(wkI) - 1)))
                     all_samps <- list()
                     
                     # Four different component models are fit in this for loop
                     # After the models are fit, expectation pieces of Monte Carlo CRPS are estimated
                     for (w in 5:(length(wkI) - 1)) {
                       anweek <- anweek + 1  
                       sdat <- wkI[1:w]
                       
                       #fit naive model
                       naive <- naive(sdat)
                       naive_m <- naive$mean[1]
                       naive_se <- sqrt(naive$model$sigma2)
                       samp_naive <- rnorm(drawn, naive_m, naive_se)
                       samp_naive[samp_naive <= 0] <- 0
                       naive_md <- mean(abs(samp_naive - wkI[w + 1]))
                       
                       #fit arima model
                       arima <- auto.arima(sdat)
                       arima_m <- forecast(arima)$mean[1]
                       arima_se <- sqrt(arima$sigma2)
                       samp_arima <- rnorm(drawn, arima_m, arima_se)
                       samp_arima[samp_arima <= 0] <- 0
                       arima_md <- mean(abs(samp_arima - wkI[w + 1]))
                       
                       
                       #fit hierarchical Bayesian model
                       m0 <- c(-4.6, -1, 15, 1.5, 1.5)
                       C0 <- diag(c(2, 1, 7, 2, 2))
                       stan_dat <- list(
                         n_weeks = length(sdat),
                         n_params = length(m0),
                         ili = sdat/N + .000001,
                         weeks = 1:length(sdat),
                         m0 = m0,
                         C0 = C0,
                         sigma_kappa = 1000
                       )
                       
                       
                       fit_asg <- asg_mod$sample(data = stan_dat, chains = 1,
                                                 iter_warmup = warm,
                                                 iter_sampling = drawn,
                                                 adapt_delta = .99)
                       
                       
                       
                       asg_draws <- fit_asg$draws(format = "df") %>% 
                         select(contains("pred_ili")) %>% 
                         as.data.frame()
                       
                       samp_asg <- asg_draws[,w + 1]*N
                       asg_md <- mean(abs(samp_asg - wkI[w + 1]))
                       
                       
                       #fit Bayesian SIR model
                       stan_dat <- list(
                         n_weeks = length(sdat),
                         weeks = 1:length(sdat),
                         ts = 1:length(sdat),
                         S0 = .9,
                         t0 = 0,
                         ili = sdat/N + .000001,
                         rho_mu = .68,
                         rho_sigma = .08,
                         beta_mu = .8,
                         beta_sigma = .3,
                         I0_mu = .005,
                         I0_sigma = .0015,
                         sigma_kappa = 1000
                       )
                       
                       fit_sir <- sir_mod$sample(data = stan_dat, chains = 1,
                                                 iter_warmup = warm,
                                                 iter_sampling = drawn,
                                                 adapt_delta = .99)
                       
                       # fit_sir <- sir_mod$variational(data = stan_dat, output_samples = drawn)
                       
                       
                       sir_draws <- fit_sir$draws(format = "df") %>% 
                         select(contains("pred_ili")) %>% 
                         as.data.frame()
                       
                       samp_sir <- sir_draws[,w + 1]*N
                       sir_md <- mean(abs(samp_sir - wkI[w + 1]))
                       
                       samp_mean <- rnorm(drawn, mean(sdat), sd(sdat))
                       mean_md <- mean(abs(samp_mean - wkI[w + 1]))
                       
                       
                       # estimate expected values of absolute difference between each model
                       all_draws <- cbind(samp_naive, samp_arima, samp_asg, samp_sir, samp_mean)
                       all_samps[[anweek]] <- all_draws
                       abs_diffs <- matrix(NA, nrow = 4, ncol = 4)
                       for (i in 1:4) {
                         for (j in i:4) {
                           mabs <- mean(abs(sample(all_draws[,i], drawn) - 
                                              sample(all_draws[,j], drawn)))
                           abs_diffs[i,j] <- mabs
                           abs_diffs[j,i] <- mabs
                         }
                       }
                       
                       all_mses <- c(naive_md, arima_md, asg_md, sir_md)#, mean_md)
                       all_mse[,anweek] <- all_mses
                       absdiff_arr[,,anweek] <- abs_diffs
                       
                     }
                       
                       
                     all_wts <- data.frame()
                     for (d in 2:ncol(all_mse)) {
                       
                       etad[1] <- 1
                       etad[2] <- 1
                       
                       for (alpha in c(1,50)) {
                         wts <- try(learning_rate(etas, d-1, mse_mat = all_mse, 
                                                  absdiff_arr = absdiff_arr, 
                                                  mod = mod, power = 1, return_wts = "draws",
                                                  tweight = .98,
                                                  alpha = alpha))
                         
                         
                         colnames(wts) <- c("naive", "arima", "asg", "sir")
                         
                         wts_long <- wts %>% 
                           pivot_longer(everything(), names_to = "model", 
                                        values_to = "draw") %>% 
                           group_by(model) %>% 
                           summarise(
                             wt = mean(draw),
                             upp = quantile(draw, probs = 0.975),
                             low = quantile(draw, probs = 0.025)
                           ) %>% 
                           mutate(week = d, alpha = alpha)
                         
                         all_wts <- rbind(all_wts, wts_long)
                         
                     }
                         
                         
                         
                    }
                       
  
                     
                     
                     all_wts <- all_wts %>% 
                       mutate(step = step, rep = replicate)
                     
                     saveRDS(all_wts, "test.rds")
                     all_wts
                     
                     
                   }

write.csv(sir_res, paste0("sir_wts/seq_", step, ".csv"))
