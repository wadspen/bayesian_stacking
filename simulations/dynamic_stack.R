# library(TruncatedNormal)
library(tmvtnorm)
library(cmdstanr)
library(distr)
library(dplyr)
library(insight)
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


gibbmod <- cmdstan_model(stan_file = '../stan_models/simple_mix_norm_crps_T.stan')

start <- 1
methods <- c("bma", "avs", "sgp", "eqw")
tw <- .65
true_wt <- c(tw, 1-tw)
log_tr_wt <- log(true_wt) + 7
tmus <- c(3, 6.5)
tsigmas <- rep(1, length(tmus))



mus <- c(0,2,4,6,8,10)
C <- length(mus)
sigmas <- rep(1, C)
tweight <- .98
ms_et <- 15

T <- 100
M <- length(true_wt)
t <- 1:T

reps <- 500


stacks <- foreach(replicate = 1:reps,
                  .packages = c("cmdstanr", "dplyr", "tidyr")
                  ,.errorhandling = "remove"
                  ,.combine = rbind) %dopar% {


  zeta <- matrix(NA, nrow = T, ncol = M)
  Sigma <- diag(rep(.1,M))
  
  
  zeta[1,] <- MASS::mvrnorm(1, log_tr_wt, Sigma = Sigma)
  
  for (i in 2:T) {
    zeta[i,] <- MASS::mvrnorm(1, zeta[i-1,] + .02, Sigma)
  }
  
  omega <- matrix(NA, nrow = T, ncol = M)
  for (i in 1:T) {
    omega[i,] <- exp(zeta[i,])/sum(exp(zeta[i,]))
  }
  
  
  y <- c()
  for (i in t) {
    
    mean <- sample(tmus, 1, prob = omega[i,])
    
    y[i] <- rnorm(1, mean, tsigmas)
  }
  
  
  all_alphas <- array(NA, dim = c(C,C,T))
  all_betas <- matrix(NA, nrow = T, ncol = C)
  
  for (i in t) {
    mean_diff <- outer(mus, mus, "-")
    var_sum <- outer(sigmas^2, sigmas^2, "+")
    param_func <- array(c(mean_diff, var_sum), dim = c(C, C, 2))
    
    all_alphas[,,i] <- apply(param_func, MARGIN = c(1,2), FUN = alphaik)
    all_betas[i,] <- sapply(y[i], FUN = betai, mu = mus, sigma = sigmas)
  }
  
  
  
  
  pb <- c()
  pa <- c()
  ps <- c()
  pm <- c()
  cb <- c()
  ca <- c()
  cs <- c()
  cm <- c()
  lb <- c()
  la <- c()
  ls <- c()
  lm <- c()
  
  
  for (d in start:(nrow(all_betas) - 1)) {
    
  #########################################
  ###################BMA###################
  #########################################
    wpmp <- sapply(mus, FUN = function(x) {
      prod(dnorm(x = y[1:d], x))
    })
    wpmp <- wpmp/sum(wpmp)
    
    
  
  #########################################
  ###################AVS###################
  #########################################
  
    etas <- seq(.001, 3, length.out = 30)
    min_eta <- c()
    for (n in 1:(d - 1)) {
      ylfo <- y[1:n]
      
      avcrpss <- c()
      for(i in 1:length(etas)) {
        et <- etas[i]
        wavs <- c()
        for (m in 1:C) {
          wavs[m] <- (1/C)*exp(-et*sum(
            tweight^(n:1 - 1)*scoringRules::crps(ylfo,
                               family = "norm", mean = mus[m], sd = 1)))
        }
        wavs <- wavs/sum(wavs)
        avcrpss[i] <- mean(all_crps(y[n + 1], mus, sigmas, ws = wavs))
      }
      min_eta[n] <- etas[which.min(avcrpss)]
    }
    avs_et <- mean(min_eta)
  
    wavs <- c()
    for (m in 1:C) {
      wavs[m] <- (1/C)*exp(-avs_et*sum(
        tweight^(d:1 - 1)*scoringRules::crps(y[1:d],
                           family = "norm", mean = mus[m], sd = 1)))
    }
    wavs <- wavs/sum(wavs)
  
    
  
  #########################################
  ###################SGP###################
  #########################################

   if (d == 1) {
	   betas <- matrix(all_betas[1:d,], nrow = 1)
	   alphas <- array(all_alphas[,,1:d], dim = c(C,C,1))
   } else {
	   betas <- all_betas[1:d,]
	   alphas <- all_alphas[,,1:d]
   }

   stan_dat <- list(
      T = d,
      y = y[1:d],
      num_comp = nrow(all_alphas),
      eta = ms_et,
      alphas = aperm(alphas),
      betas = betas,
      alpha = rep(1, length(sigmas)),
      wts = rep(1/length(mus), length(mus)),
      tweight = tweight
    )
    
    
    
    # fit <- gibbmod$variational(data = stan_dat)
    fit <- gibbmod$sample(data = stan_dat, chains = 1, 
                          iter_warmup = 5000,
                          iter_sampling = 5000, 
                          init = list(list(omega = rep(1/C, C))))
    
    draws <- fit$draws(variables = "omega", format = "df") %>%
      dplyr::select(contains("omega"))
    
  
    wmsgp <- apply(draws, MARGIN = 2, FUN = mean)
    wmsgp <- wmsgp/sum(wmsgp)
  
    pb[d] <- pmixnorm(y[d + 1], mus, sigmas, wpmp)
    pa[d] <- pmixnorm(y[d + 1], mus, sigmas, wavs)
    ps[d] <- pmixnorm(y[d + 1], mus, sigmas, wmsgp)
    pm[d] <- pmixnorm(y[d + 1], mus, sigmas, rep(1/C, C))
    
    cb[d] <- scoringRules::crps(y[d + 1], family = "mixnorm",
                                m = matrix(mus, nrow = 1),
                                s = matrix(sigmas, nrow = 1),
                                w = matrix(wpmp, nrow = 1))
    
    
    ca[d] <- scoringRules::crps(y[d + 1], family = "mixnorm",
                                m = matrix(mus, nrow = 1),
                                s = matrix(sigmas, nrow = 1),
                                w = matrix(wavs, nrow = 1))
     
    cs[d] <- scoringRules::crps(y[d + 1], family = "mixnorm",
                                m = matrix(mus, nrow = 1),
                                s = matrix(sigmas, nrow = 1),
                                w = matrix(wmsgp, nrow = 1))
     
    cm[d] <- scoringRules::crps(y[d + 1], family = "mixnorm",
                                m = matrix(mus, nrow = 1),
                                s = matrix(sigmas, nrow = 1),
                                w = matrix(rep(1/C, C), nrow = 1))
    
    
    lb[d] <- scoringRules::logs(y[d + 1], family = "mixnorm",
                                m = matrix(mus, nrow = 1),
                                s = matrix(sigmas, nrow = 1),
                                w = matrix(wpmp, nrow = 1))
    
    
    la[d] <- scoringRules::logs(y[d + 1], family = "mixnorm",
                                m = matrix(mus, nrow = 1),
                                s = matrix(sigmas, nrow = 1),
                                w = matrix(wavs, nrow = 1))
    
    ls[d] <- scoringRules::logs(y[d + 1], family = "mixnorm",
                                m = matrix(mus, nrow = 1),
                                s = matrix(sigmas, nrow = 1),
                                w = matrix(wmsgp, nrow = 1))
    
    lm[d] <- scoringRules::logs(y[d + 1], family = "mixnorm",
                                m = matrix(mus, nrow = 1),
                                s = matrix(sigmas, nrow = 1),
                                w = matrix(rep(1/C, C), nrow = 1))
    
    
    
    
    
    
  }
  
  
  #mcrps <- c(mean(cb, na.rm = TRUE), mean(ca, na.rm = TRUE), 
  #          mean(cs, na.rm = TRUE), mean(cm, na.rm = TRUE))
  #
  #mlogs <- c(mean(lb, na.rm = TRUE), mean(la, na.rm = TRUE), 
  #          mean(ls, na.rm = TRUE), mean(lm, na.rm = TRUE))
  #
  #uwd1s <- c(unit_wass_dist(ecdf(pb)),
  #              unit_wass_dist(ecdf(pa)),
  #              unit_wass_dist(ecdf(ps)),
  #              unit_wass_dist(ecdf(pm)))

  crpss <- c(cb, ca, cs, cm)
  logss <- c(lb, la, ls, lm)
  pits <- c(pb, pa, ps, pm)
  times <- 1:length(cb)
  methodsl <- rep(methods, each = length(cb))
  
  
  data.frame(replicate = replicate, method = methodsl,
             eta = avs_et, time = times, mcrps = crpss, mlogs = logss, pit = pits)
  
  
}


saveRDS(stacks, "dynamic_stacks_all_weeks.rds")












