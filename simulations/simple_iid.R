source("./stack_functions.R")
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




gibbmod <- cmdstan_model(stan_file = '../stan_models/simple_mix_norm_crps.stan')

methods <- c("eqw", "bma", "avs", "msgp")
samp_sizes <- c(10, 20, 50, 100, 200, 400, 800)
#models
mus <- c(0,2,4,6,8,10)
C <- length(mus)
sigmas <- rep(1, C)
reps <- 500
total_samp <- 1000
stacks <- foreach(replicate = 1:reps,
                    .packages = c("cmdstanr", "dplyr", "tidyr")
                    ,.errorhandling = "remove"
                    ,.combine = rbind) %:%
                foreach(N = samp_sizes, .combine = rbind) %dopar% {

      #data
      #N2 <- 10
      #N <- ceiling(N2/2)
      N2 <- total_samp + N
      tmus <- c(3,6.5)
      tw <- .65
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
      ######EQW#####
      ##############
      weq <- rep(1/length(mus), length(mus))


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
      etas <- seq(.001, 1.5, length.out = 30)
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
      #etas <- seq(.5, 20, length.out = 15)
      #etas <- 15
      #min_eta <- c()
      #for (n in 1:N) {
      #  yloo <- y[-n]
      #  mscrpss <- c()
      #  
      #  all_alphas <- array(NA, dim = c(C,C,length(yloo)))
      #  all_betas <- matrix(NA, nrow = length(yloo), ncol = C)
      #  for (l in 1:length(yloo)) {
      #    
      #    mean_diff <- outer(mus, mus, "-")
      #    var_sum <- outer(sigmas^2, sigmas^2, "+")
      #    param_func <- array(c(mean_diff, var_sum), dim = c(C, C, 2))
      #    
      #    all_alphas[,,l] <- apply(param_func, MARGIN = c(1,2), FUN = alphaik)
      #    
      #    all_betas[l,] <- sapply(yloo[l], FUN = betai, mu = mus, sigma = sigmas)
      #    
      #  }
      #  
      #  
      #  #etas <- seq(.5, 20, length.out = 15)
      #  mscrpss <- c()
      #  for (i in 1:length(etas)) {
      #    et <- etas[i]
      #    stan_dat <- list(
      #      N = length(yloo),
      #      y = yloo,
      #      num_comp = nrow(all_alphas[,,1]),
      #      eta = et,
      #      alphas = aperm(all_alphas[,,1]),
      #      betas = all_betas[,],
      #      alpha = rep(1, C),
      #      mus = mus,
      #      sigmas = sigmas,
      #      sigma_s = rep(3,5)
      #    )
      #    
      #    #fit <- gibbmod$variational(data = stan_dat)
      #    fit <- gibbmod$sample(data = stan_dat, chains = 1, 
      #                      iter_warmup = 500,
      #                      iter_sampling = 1000, 
      #                      init = list(list(omega = rep(1/C, C))))
      #    draws <- fit$draws(variables = "omega", format = "df") %>%
      #      dplyr::select(contains("omega"))
      #    
      #    
      #    wmsgp <- apply(draws, MARGIN = 2, FUN = mean)
      #    wmsgp <- wmsgp/sum(wmsgp)
      #  
      #  
      #  mscrpss[i] <- mean(all_crps(y, mus, sigmas, ws = wmsgp)) + .0002*et
      #  }
      #  min_eta[n] <- etas[which.min(mscrpss)]
      #}
      #ms_et <- mean(min_eta)
      ms_et <- 15 
      
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
      fit <- gibbmod$sample(data = stan_dat, chains = 1, 
                            iter_warmup = 10000,
                            iter_sampling = 50000, 
                            init = list(list(omega = rep(1/C, C))))
      draws <- fit$draws(variables = "omega", format = "df") %>%
        dplyr::select(contains("omega"))
        
        
      wmsgp <- apply(draws, MARGIN = 2, FUN = mean)
      wmsgp <- wmsgp/sum(wmsgp)
        
      
      
      #########################
      ########Evaluate#########
      #########################
      
      mcrps <- c(mean(all_crps(y_test, mus, sigmas, ws = weq)),
		 mean(all_crps(y_test, mus, sigmas, ws = wpmp)),
                 mean(all_crps(y_test, mus, sigmas, ws = wavs)),
                 mean(all_crps(y_test, mus, sigmas, ws = wmsgp)))
      
      mlogs <- c(mean(all_logs(y_test, mus, sigmas, ws = weq)),
		 mean(all_logs(y_test, mus, sigmas, ws = wpmp)),
                 mean(all_logs(y_test, mus, sigmas, ws = wavs)),
                 mean(all_logs(y_test, mus, sigmas, ws = wmsgp))) 
      
      peq <- c()
      pbma <- c()
      pavs <- c()
      pmsgp <- c()
      for (i in 1:length(y_test)) {
	peq[i] <- pmixnorm(y_test[i], mus, sigmas, weq)
        pbma[i] <- pmixnorm(y_test[i], mus, sigmas, wpmp)
        pavs[i] <- pmixnorm(y_test[i], mus, sigmas, wavs)
        pmsgp[i] <- pmixnorm(y_test[i], mus, sigmas, wmsgp)
      }
      
      uwd1 <- c(unit_wass_dist(ecdf(weq)),
		unit_wass_dist(ecdf(pbma)),
                unit_wass_dist(ecdf(pavs)),
                unit_wass_dist(ecdf(pmsgp)))
      
      min_etas <- c(NA, NA, avs_et, ms_et)
      weights <- as.data.frame(rbind(weq, wpmp, wavs, wmsgp))
      colnames(weights) <- paste0("comp", 1:length(mus)) 
      res <- data.frame(replicate = replicate, N = N, method = methods, 
                 eta = min_etas, mcrps = mcrps, mlogs = mlogs, uwd1 = uwd1)

      res <- cbind(res, weights)
      res

}

saveRDS(stacks, "iid_stacks_res.rds")
