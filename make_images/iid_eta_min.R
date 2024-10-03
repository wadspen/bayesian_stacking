source("./simulations/stack_functions.R")
library(dplyr)
library(tidyr)
library(cmdstanr)





gibbmod <- cmdstan_model(stan_file = './stan_models/simple_mix_norm_crps.stan')

methods <- c("bma", "avs", "msgp")
samp_sizes <- c(10, 20, 50, 100, 200, 400, 800)
#models
mus <- c(0,2,4,6,8,10)
C <- length(mus)
sigmas <- rep(1, C)
reps <- 500
total_samp <- 1000
    N <- 10
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
    #####AVS######
    ##############
    etas <- seq(.01, 1.5, length.out = 30)
    min_eta <- c()
      yloo <- y
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
      plot(avcrpss ~ etas, type = "l")
      
      data.frame(y = avcrpss, x = etas) %>% 
        filter(etas > .01) %>% 
        ggplot() +
        geom_line(aes(x = x, y = y), size = 2) +
        xlab(expression(eta)) +
        ylab(expression(CRPS[eta])) +
        theme_bw() +
        theme(axis.text.y=element_text(size=17),
              axis.text.x=element_text(size = 17),
              axis.title=element_text(size=24),
              strip.text.y = element_text(size = 14,),
              strip.text.x = element_text(size = 12),
              legend.position = "none")
    
    
    
    ##############
    #####MSGP#####
    ##############
    etas <- seq(.5, 20, length.out = 15)
    # etas <- 15
    min_eta <- c()
     yloo <- y
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


     #etas <- seq(.5, 20, length.out = 15)
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

       #fit <- gibbmod$variational(data = stan_dat)
       fit <- gibbmod$sample(data = stan_dat, chains = 1,
                         iter_warmup = 5000,
                         iter_sampling = 5000,
                         init = list(list(omega = rep(1/C, C))))
       draws <- fit$draws(variables = "omega", format = "df") %>%
         dplyr::select(contains("omega"))


       wmsgp <- apply(draws, MARGIN = 2, FUN = mean)
       wmsgp <- wmsgp/sum(wmsgp)


     mscrpss[i] <- mean(all_crps(y, mus, sigmas, ws = wmsgp)) + .0002*et
     }
     
     
     data.frame(y = mscrpss, x = etas) %>% 
       filter(etas > .01) %>% 
       ggplot() +
       geom_line(aes(x = x, y = y), size = 2) +
       xlab(expression(eta)) +
       # ylab(expression(CRPS[eta])) +
       ylab("") +
       theme_bw() +
       theme(axis.text.y=element_text(size=17),
             axis.text.x=element_text(size = 17),
             axis.title=element_text(size=24),
             strip.text.y = element_text(size = 14,),
             strip.text.x = element_text(size = 12),
             legend.position = "none")
    
    