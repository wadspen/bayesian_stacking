library(SimInf)
library(forecast)
library(ggplot2)
library(cmdstanr)
library(dplyr)
library(car)
library(stats)
library(scoringutils)
source("./stack_functions.R")

asg_mod <- cmdstan_model(stan_file = '../stan_models/asg.stan')
sir_mod <- cmdstan_model(stan_file = '../stan_models/sir.stan')
mod <- cmdstan_model(stan_file = '../stan_models/emp_mix_crps_time_weight.stan')
drawn <- 1000
N <- 1000
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


I <- result@U[2,]
time <- result@tspan
wktime <- seq(5, days, by = 7)
wkI <- I[wktime]
plot(wkI~wktime)
models <- c("naive", "arima", "asg", "sir")#, "mean")
scores <- data.frame()
absdiff_arr <- array(NA, dim = c(4,4,length(5:(length(wkI) - 1))))
anweek <- 0
all_mse <- matrix(NA, nrow = 4, ncol = length(5:(length(wkI) - 1)))
for (w in 5:(length(wkI) - 1)) {
  anweek <- anweek + 1  
  sdat <- wkI[1:w]
  naive <- naive(sdat)
  naive_m <- naive$mean[1]
  naive_se <- sqrt(naive$model$sigma2)
  samp_naive <- rnorm(drawn, naive_m, naive_se)
  samp_naive[samp_naive <= 0] <- 0
  naive_md <- mean(abs(samp_naive - wkI[w + 1]))
  
  
  
  arima <- auto.arima(sdat)
  arima_m <- forecast(arima)$mean[1]
  arima_se <- sqrt(arima$sigma2)
  samp_arima <- rnorm(drawn, arima_m, arima_se)
  samp_arima[samp_arima <= 0] <- 0
  arima_md <- mean(abs(samp_arima - wkI[w + 1]))
  
  m0 <- c(-4.6, -1, 15, 1.5, 1.5)
  C0 <- diag(c(2, 1, 7, 2, 2))
  stan_dat <- list(
    n_weeks = length(sdat),
    n_params = length(m0),
    ili = sdat/N,
    weeks = 1:length(sdat),
    m0 = m0,
    C0 = C0,
    sigma_kappa = 1000
  )
  
  fit_asg <- asg_mod$sample(data = stan_dat, chains = 1,
                        iter_warmup = drawn,
                        iter_sampling = drawn,
                        adapt_delta = .99)
  
  # fit_asg <- asg_mod$variational(data = stan_dat, output_samples = drawn)
  
  asg_draws <- fit_asg$draws(format = "df") %>% 
    select(contains("pred_ili")) %>% 
    as.data.frame()
  
  samp_asg <- asg_draws[,w + 1]*N
  asg_md <- mean(abs(samp_asg - wkI[w + 1]))
  
  
  
  stan_dat <- list(
    n_weeks = length(sdat),
    weeks = 1:length(sdat),
    ts = 1:length(sdat),
    S0 = .9,
    t0 = 0,
    ili = sdat/N,
    rho_mu = .68,
    rho_sigma = .08,
    beta_mu = .8,
    beta_sigma = .3,
    I0_mu = .005,
    I0_sigma = .0015,
    sigma_kappa = 1000
  )
  
  fit_sir <- sir_mod$sample(data = stan_dat, chains = 1,
                        iter_warmup = drawn,
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
  
  all_draws <- cbind(samp_naive, samp_arima, samp_asg, samp_sir, samp_mean)
  
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
  

  
  # samp_ens <- sample(c(samp_naive, samp_arima, samp_asg, samp_sir), drawn)
  crpss <- c(crps_sample(wkI[w + 1], samp_naive),
             crps_sample(wkI[w + 1], samp_arima),
             crps_sample(wkI[w + 1], samp_asg),
             crps_sample(wkI[w + 1], samp_sir))#,
             # crps_sample(wkI[w + 1], samp_mean))
  
  
  logss <- c(logs_sample(wkI[w + 1], samp_naive),
             logs_sample(wkI[w + 1], samp_arima),
             logs_sample(wkI[w + 1], samp_asg),
             logs_sample(wkI[w + 1], samp_sir))
  
  # logss <- c(exp(-logs_sample(wkI[w + 1], samp_naive)),
  #            exp(-logs_sample(wkI[w + 1], samp_arima)),
  #            exp(-logs_sample(wkI[w + 1], samp_asg)),
  #            exp(-logs_sample(wkI[w + 1], samp_sir)))
  
  print(which.min(crpss))
  score <- data.frame(week = w, model = models, crps = crpss, logs = logss)
  scores <- rbind(scores, score)
  
  

}


crpsh1 <- scores %>% 
  select(-logs) %>%
  pivot_wider(names_from = "model", values_from = "crps") %>% 
  select(-week) %>% 
  as.matrix() %>% 
  t()

logsh1 <- scores %>% 
  select(-crps) %>% 
  pivot_wider(names_from = "model", values_from = "logs") %>% 
  select(-week) %>% 
  as.matrix() %>% 
  t()

C <- nrow(crpsh1)
bma_wt <- matrix(NA, nrow = C, ncol = ncol(all_mse))
avs_wt <- matrix(NA, nrow = C, ncol = ncol(all_mse))
bma_wt[,1] <- rep(1/C, C)
avs_wt[,1] <- rep(1/C, C)	

alpha <- .98
d <- ncol(crpsh1)
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
          alpha^(n:1 - 1)*crpsh1[m,1:n]))
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
    # wbma[m] <- (1/C)*prod(logsh1[m, 1:d]^(alpha^(d:1 - 1)))
    # wbma[m] <- (1/C)*prod(exp(-logsh1[m,1:d])^(alpha^(d:1 - 1)))
    wbma[m] <- (1/C)*exp(sum(
      alpha^(d:1 - 1)*-logsh1[m,1:d]))
  }
  bma_wt[,d + 1] <- wbma/sum(wbma)
  wavs <- c()
  for (m in 1:C) {
    wavs[m] <- (1/C)*exp(-avs_et*sum(
      alpha^(d:1 - 1)*crpsh1[m,1:d]))
  }
  avs_wt[,d + 1] <- wavs/sum(wavs)
}

scores %>% 
  group_by(model) %>% 
  summarise(n())

scores %>% 
  group_by(model) %>% 
  # filter(week == 9) %>%
  summarise(mean(crps))

scores %>% 
  # filter(model != "sir") %>% 
  ggplot(aes(x = week, y = crps, colour = model)) + 
  geom_line()

scores %>% 
  arrange(week, crps) %>% 
  mutate(ind = 1) %>% 
  group_by(week) %>% 
  mutate(rank = cumsum(ind)) %>% 
  group_by(model) %>% 
  summarise(median(rank))




eq_wt <- rep(1, 4)/4
dim(absdiff_arr)
etas <- seq(.1, 5, length.out = 20)
etas <- seq(-3, 5, length.out = 20)
# etas <- c(100, 150)
# etas <- c(3, 7)
# etas <- .5
etad <- c()
weight <- matrix(NA, nrow = 4, ncol = dim(absdiff_arr)[3])
weight[,1] <- eq_wt
stack_crps <- c()
pp_stack_crps <- c()
for (d in 32:ncol(all_mse)) {
  
  
  
  etad[1] <- 1
  
  
  if (d == 1) {
    stack_crps[d] <- mix_mat_crps(weight[,d], all_mse[,d], 
                                  absdiff_arr[,,d])
    wts <- weight[,d]
  } else {
    ev_grid <- c()
    for (i in 1:length(etas)) {
      ev_grid[i] <- 
        try(learning_rate(etas[i], d-1, mse_mat = all_mse, 
                          absdiff_arr = absdiff_arr, 
                          mod = mod, power = 1, tweight = .98,
                          alpha = 1))
      
    }
    etad[d] <- exp(etas[which.min(ev_grid)])
    wts <- try(learning_rate(log(etad[d]), d, mse_mat = all_mse, 
                             absdiff_arr = absdiff_arr, 
                             mod = mod, power = 1, return_wts = "draws",
                             tweight = .98,
                             alpha = 1))
    
    
    pp_crps <- c()
    for (i in 1:nrow(wts)) {
      wt <- wts[i,]
      wt[wt < 0] <- 0
      wt <- wt/sum(wt)
      wt <- unlist(as.vector(wt))
      pp_crps[i] <- mix_mat_crps(wt, all_mse[,d], 
                                    absdiff_arr[,,d])
    }
    
    pp_stack_crps[d] <- mean(pp_crps)
    
    wts <- apply(wts, MARGIN = 2, FUN = mean)
    wts[wts < 0] <- 0
    wts <- wts/sum(wts)
    weight[,d] <- wts
    
    stack_crps[d] <- mix_mat_crps(weight[,d], all_mse[,d], 
                                  absdiff_arr[,,d])
  }
  
}
stack_crps
weight
mean_crps <- mix_mat_crps(eq_wt, all_mse, absdiff_arr)

avs_crps <- c()
for (i in 1:ncol(all_mse)) {
  avs_crps[i] <- mix_mat_crps(avs_wt[,i], all_mse[,i], absdiff_arr[,,i])
}

bma_crps <- c()
for (i in 1:ncol(all_mse)) {
  bma_crps[i] <- mix_mat_crps(bma_wt[,i], all_mse[,i], absdiff_arr[,,i])
}





