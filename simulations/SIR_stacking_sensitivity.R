.libPaths("~/rlibs")
library(SimInf)
library(forecast)
library(tidyr)
library(ggplot2)
library(cmdstanr)
library(dplyr)
library(car)
library(stats)
library(scoringRules)
source("./stack_functions.R")
library(parallel)
library(doParallel)
library(doMC)

# Run with: Rscript SIR_stacking_sensitivity.R <sequence_id>
# Edit these grids and reps to control the sensitivity analysis.
args <- commandArgs(trailingOnly = TRUE)
step <- if (length(args) >= 1) as.integer(args[1]) else 1L
eta_grid <- c(0.25, 1, 4)
dirichlet_grid <- c(0.25, 1, 4)
discount_grid <- c(0.8, 0.9, 1)
reps <- 100
drawn <- 2000
warm <- 2000
n.cores <- min(max(1, parallel::detectCores() - 1), reps)

my.cluster <- makeCluster(n.cores, type = "PSOCK")
doParallel::registerDoParallel(cl = my.cluster)
registerDoMC(cores = n.cores)

asg_mod <- cmdstan_model(stan_file = "../stan_models/asg.stan")
sir_mod <- cmdstan_model(stan_file = "../stan_models/sir.stan")
stack_mod <- cmdstan_model(stan_file = "../stan_models/emp_mix_crps_time_weight.stan")

sensitivity_res <- foreach(replicate = seq_len(reps),
  .packages = c("cmdstanr", "stringr", "scoringRules", "lubridate", "dplyr", "tidyr"),
  .errorhandling = "stop", .combine = rbind) %dopar% {
    N <- 800
    I0 <- 12
    days <- 300
    model <- SIR(u0 = data.frame(S = N, I = I0, R = 0), tspan = 1:days,
                 beta = rnorm(1, .06, .01), gamma = rnorm(1, 0.034, .009))
    result <- run(model)
    I <- result@U[2,]
    wktime <- seq(5, days, by = 7)
    wkI <- I[wktime]
    absdiff_arr <- array(NA_real_, dim = c(4, 4, length(5:(length(wkI) - 1))))
    all_mse <- matrix(NA_real_, nrow = 4, ncol = length(5:(length(wkI) - 1)))
    all_samps <- vector("list", length = ncol(all_mse))

    for (idx in seq_along(5:(length(wkI) - 1))) {
      w <- idx + 4
      sdat <- wkI[1:w]
      naive_fit <- naive(sdat)
      samp_naive <- pmax(0, rnorm(drawn, naive_fit$mean[1], sqrt(naive_fit$model$sigma2)))
      arima_fit <- auto.arima(sdat)
      arima_fc <- forecast(arima_fit)
      samp_arima <- pmax(0, rnorm(drawn, arima_fc$mean[1], sqrt(arima_fit$sigma2)))

      asg_dat <- list(n_weeks = length(sdat), n_params = 5,
        ili = sdat/N + .000001, weeks = seq_along(sdat),
        m0 = c(-4.6, -1, 15, 1.5, 1.5), C0 = diag(c(2, 1, 7, 2, 2)), sigma_kappa = 1000)
      fit_asg <- asg_mod$sample(data = asg_dat, chains = 1, iter_warmup = warm,
                                iter_sampling = drawn, adapt_delta = .99)
      asg_draws <- fit_asg$draws(format = "df") %>% select(contains("pred_ili"))
      samp_asg <- asg_draws[[w + 1]] * N

      sir_dat <- list(n_weeks = length(sdat), weeks = seq_along(sdat), ts = seq_along(sdat),
        S0 = .9, t0 = 0, ili = sdat/N + .000001, rho_mu = .68, rho_sigma = .08,
        beta_mu = .8, beta_sigma = .3, I0_mu = .005, I0_sigma = .0015, sigma_kappa = 1000)
      fit_sir <- sir_mod$sample(data = sir_dat, chains = 1, iter_warmup = warm,
                                iter_sampling = drawn, adapt_delta = .99)
      sir_draws <- fit_sir$draws(format = "df") %>% select(contains("pred_ili"))
      samp_sir <- sir_draws[[w + 1]] * N

      samples <- list(samp_naive, samp_arima, samp_asg, samp_sir)
      all_samps[[idx]] <- do.call(cbind, samples)
      all_mse[, idx] <- vapply(samples, function(s) crps_sample(wkI[w + 1], s), numeric(1))
      for (i in 1:4) for (j in i:4) {
        dij <- mean(abs(sample(samples[[i]], drawn, replace = TRUE) -
                        sample(samples[[j]], drawn, replace = TRUE)))
        absdiff_arr[i, j, idx] <- dij
        absdiff_arr[j, i, idx] <- dij
      }
    }

    grid <- expand.grid(eta = eta_grid, dirichlet_alpha = dirichlet_grid,
                        discount = discount_grid)
    result_rows <- vector("list", nrow(grid))
    for (g in seq_len(nrow(grid))) {
      weights <- matrix(NA_real_, nrow = 4, ncol = ncol(all_mse))
      weights[, 1] <- rep(1/4, 4)
      for (d in 2:ncol(all_mse)) {
        draws <- learning_rate(eta = grid$eta[g], i = d - 1,
          mse_mat = all_mse, absdiff_arr = absdiff_arr, mod = stack_mod,
          tweight = grid$discount[g], power = 1, return_wts = "draws",
          alpha = grid$dirichlet_alpha[g])
        wt <- colMeans(draws)
        wt <- pmax(wt, 0)
        weights[, d] <- wt / sum(wt)
      }
      crps <- vapply(seq_len(ncol(all_mse)), function(d)
        mix_mat_crps(weights[, d], all_mse[, d], absdiff_arr[, , d]), numeric(1))
      logs <- vapply(seq_len(ncol(all_mse)), function(d) {
        ens <- vapply(seq_len(drawn), function(k) {
          component <- sample.int(4, 1, prob = weights[, d])
          sample(all_samps[[d]][, component], 1)
        }, numeric(1))
        logs_sample(wkI[d + 5], ens)
      }, numeric(1))
      result_rows[[g]] <- data.frame(seq = step, rep = replicate,
        eta = grid$eta[g], dirichlet_alpha = grid$dirichlet_alpha[g],
        discount = grid$discount[g], time = seq_along(crps), crps = crps, logs = logs)
    }
    do.call(rbind, result_rows)
  }

stopCluster(my.cluster)
out_file <- paste0("sir_res/sensitivity_seq_", step, ".rds")
saveRDS(sensitivity_res, out_file)
write.csv(sensitivity_res, sub("\\.rds$", ".csv", out_file), row.names = FALSE)
