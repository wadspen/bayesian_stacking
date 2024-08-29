source("./simulations/stack_functions.R")
library(dplyr)
library(tidyr)
library(cmdstanr)
library(parallel)
library(doParallel)
library(doMC)
n.cores <- detectCores()
my.cluster <- makeCluster(n.cores, type = "PSOCK")
doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()
foreach::getDoParWorkers()
registerDoMC(cores = n.cores)




gibbmod <- cmdstan_model(stan_file = './stan_models/simple_mix_norm_crps.stan')

methods <- c("bma", "avs", "msgp")
samp_sizes <- c(10, 20, 50, 100, 200)
#models
mus <- c(0,2,4,6,8,10)
C <- length(mus)
sigmas <- rep(1, C)


distance <- foreach(replicate = 1:reps,
                    .packages = c("cmdstanr", "evmix", "distfromq", "dplyr", "tidyr",
                                  "janitor")
                    ,.errorhandling = "remove"
                    ,.combine = rbind) %:%
                foreach(N2 = samp_sizes, .combine = rbind) %dopar% {

      #data
      #N2 <- 10
      N <- ceiling(N2/2)
      tmus <- c(3,5.5)
      tw <- .7
      all_y <- c()
      for (i in 1:N2) {
        smu <- sample(tmus, 1, prob = c(tw, 1-tw))
        all_y[i] <- rnorm(1, smu)
      }
      
      #train
      train_ind <- sample(N2, N)
      y <- all_y[train_ind]
      y_test <- all_y[-train_ind]
      
      
      ##############
      #####BMA######
      ##############
      wpmp <- sapply(mus, FUN = function(x) {
         prod(dnorm(x = y, x))
        })
      wpmp <- wpmp/sum(wpmp)
      
      
      ##############
      #####AVS######
      ##############
      etas <- seq(.01, 3, length.out = 30)
      min_eta <- c()
      for (n in 1:N) {
        yloo <- y[-n]
        avcrpss <- c()
        for(i in 1:length(etas)) {
          et <- etas[i]
          wavs <- c()
          for (m in 1:C) {
            wavs[m] <- (1/C)*exp(-et*sum(
              scoringRules::crps(yloo,
                                 family = "norm", mean = mus[m], sd = 1)))
          }
          wavs <- wavs/sum(wavs)
          avcrpss[i] <- mean(all_crps(y, mus, sigmas, ws = wavs))
        }
        min_eta[n] <- etas[which.min(avcrpss)]
      }
      avs_et <- mean(min_eta)
      
      
      wavs <- c()
      for (m in 1:C) {
        wavs[m] <- (1/C)*exp(-avs_et*sum(
            scoringRules::crps(yloo,
                               family = "norm", mean = mus[m], sd = 1)))
      }
      wavs <- wavs/sum(wavs)
      
      
      
      ##############
      #####MSGP#####
      ##############
      etas <- seq(.01, 10, length.out = 30)
      min_eta <- c()
      for (n in 1:N) {
        yloo <- y[-n]
        mscrpss <- c()
        
        all_alphas <- array(NA, dim = c(C,C,length(yloo)))
        all_betas <- matrix(NA, nrow = length(yloo), ncol = C)
        for (l in 1:length(yloo)) {
          
          mean_diff <- outer(mus, mus, "-")
          var_sum <- outer(sigmas^2, sigmas^2, "+")
          param_func <- array(c(mean_diff, var_sum), dim = c(C, C, 2))
          
          all_alphas[,,l] <- apply(param_func, MARGIN = c(1,2), FUN = alphaik)
          
          all_betas[l,] <- sapply(yloo[l], FUN = betai, mu = mus, sigma = sigmas)
          
        }
        
        
        etas <- seq(.01, 20, length.out = 30)
        mscrpss <- c()
        for (i in 1:length(etas)) {
          et <- etas[i]
          stan_dat <- list(
            N = length(yloo),
            y = yloo,
            num_comp = nrow(all_alphas[,,1]),
            eta = et,
            alphas = aperm(all_alphas[,,1]),
            betas = all_betas[,],
            alpha = rep(1, C),
            mus = mus,
            sigmas = sigmas,
            sigma_s = rep(3,5)
          )
          
          fit <- gibbmod$variational(data = stan_dat)
          # fit <- gibbmod$sample(data = stan_dat, chains = 1)
          draws <- fit$draws(variables = "omega", format = "df") %>%
            select(contains("omega"))
          
          
          wmsgp <- apply(draws, MARGIN = 2, FUN = mean)
          wmsgp <- wmsgp/sum(wmsgp)
        
        
        mscrpss[i] <- mean(all_crps(y, mus, sigmas, ws = wmsgp)) + .0002*et
        }
        min_eta[n] <- etas[which.min(mscrpss)]
      }
      ms_et <- mean(min_eta)
      
      
      all_alphas <- array(NA, dim = c(C,C,N))
      all_betas <- matrix(NA, nrow = N, ncol = C)
      for (i in 1:N) {
        
        mean_diff <- outer(mus, mus, "-")
        var_sum <- outer(sigmas^2, sigmas^2, "+")
        param_func <- array(c(mean_diff, var_sum), dim = c(C, C, 2))
        
        all_alphas[,,i] <- apply(param_func, MARGIN = c(1,2), FUN = alphaik)
        
        all_betas[i,] <- sapply(y[i], FUN = betai, mu = mus, sigma = sigmas)
        
      }
      
      
      
      stan_dat <- list(
        N = N,
        y = y,
        num_comp = nrow(all_alphas[,,1]),
        eta = et,
        alphas = aperm(all_alphas[,,1]),
        betas = all_betas[,],
        alpha = rep(1, C)
      )
      
        # fit <- gibbmod$variational(data = stan_dat)
      fit <- gibbmod$sample(data = stan_dat, chains = 1)
      draws <- fit$draws(variables = "omega", format = "df") %>%
        select(contains("omega"))
        
        
      wmsgp <- apply(draws, MARGIN = 2, FUN = mean)
      wmsgp <- wmsgp/sum(wmsgp)
        
      
      
      #########################
      ########Evaluate#########
      #########################
      
      mcrps <- c(mean(all_crps(y_test, mus, sigmas, ws = wpmp)),
                 mean(all_crps(y_test, mus, sigmas, ws = wavs)),
                 mean(all_crps(y_test, mus, sigmas, ws = wmsgp)))
      
      pbma <- c()
      pavs <- c()
      pmsgp <- c()
      for (i in 1:length(y_test)) {
        pbma[i] <- pmixnorm(y_test[i], mus, sigmas, wpmp)
        pavs[i] <- pmixnorm(y_test[i], mus, sigmas, wavs)
        pmsgp[i] <- pmixnorm(y_test[i], mus, sigmas, wmsgp)
      }
      
      uwd1 <- c(unit_wass_dist(ecdf(pbma)),
                unit_wass_dist(ecdf(pavs)),
                unit_wass_dist(ecdf(pmsgp)))
      
      min_etas <- c(NA, avs_et, ms_et)
      
      data.frame(replicate = replicate, N = N2, method = methods, 
                 eta = min_etas, mcrps = mcrps, uwd1 = uwd1)

}


