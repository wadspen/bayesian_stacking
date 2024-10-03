# library(TruncatedNormal)
library(tmvtnorm)
library(cmdstanr)
library(distr)
library(dplyr)
library(insight)
library(stringr)
library(ggplot2)
library(tidyverse)
source("./simulations/stack_functions.R")



gibbmod <- cmdstan_model(stan_file = './stan_models/simple_mix_norm_crps_T.stan')

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
          omega_plot <- cbind(omega, 1:nrow(omega))
          colnames(omega_plot) <- c(paste0("comp", 1:2), "t")
          # saveRDS(omega_plot, "./make_images/true_omega_df.rds")
          
          
          
                    
                    
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
                    
             # fit_weights <- matrix(NA, nrow = nrow(omega), ncol = C)       
                    for (d in start:(nrow(all_betas) - 1)) {
                      
                      
                      
                      
                      
                     
                      
                      
                      
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
                      
                      fit_weights[d,] <- wmsgp

                  }


             # saveRDS(fit_weights, "./make_images/fit_omega_mat.rds")
             
      # test <- fit_weights
        test <- test %>% 
          as.data.frame()
        
        colnames(test) <- paste0("comp", 1:C)
      
        
        test <- cbind(test, t = 1:nrow(test))
        
        
        
        
        omega_plot %>% 
          as.data.frame() %>% 
          pivot_longer(1:2, names_to = "comp", values_to = "weight") %>% 
          ggplot() +
          geom_line(aes(x = t, y = weight, linetype = comp), size = 1.1) +
          ylab("") +
          scale_linetype_discrete(labels = c(3, 6.5)) +
          labs(linetype = expression(paste(mu, ":  "))) +
          theme_bw() +
          theme(axis.text.y=element_text(size=17),
                axis.text.x=element_text(size = 17),
                axis.title=element_text(size=23),
                strip.text.y = element_text(size = 12,),
                strip.text.x = element_text(size = 16),
                legend.title = element_text(size = 22),
                legend.text = element_text(size = 17),
                legend.position = c(.7,.91),
                legend.direction='horizontal',
                legend.background = element_rect(size = 0.5, colour = 1))
        
        test%>% 
          pivot_longer(1:C, names_to = "comp", 
                       values_to = "weight") %>% 
          group_by(comp) %>% 
          ggplot() +
          geom_line(aes(x = t, y = weight, 
                        colour = comp, linetype = comp), size = 1.1) +
          scale_colour_discrete(labels = c(0,2,4,6,8,10)) +
          scale_linetype_discrete(labels = c(0,2,4,6,8,10)) +
          ylab("") +
          labs(colour = expression(mu),
               linetype = expression(mu)) +
          theme_bw() +
          theme(axis.text.y=element_text(size=15),
                axis.text.x=element_text(size = 19),
                axis.title=element_text(size=23),
                strip.text.y = element_text(size = 12,),
                strip.text.x = element_text(size = 16),
                legend.title = element_text(size = 22),
                legend.text = element_text(size = 18),
                legend.position = c(.81,.78),
                legend.background = element_rect(size = 0.5, colour = 1))
        
        
        
        omega_plot %>% 
          as.data.frame() %>% 
          mutate(ev = 3*comp1 + 6.5*comp2) %>% 
          ggplot() +
          geom_line(data = test_ev, aes(x = t, y = ev)) +
          geom_line(aes(x = t, y = ev))
        
        test_ev <- test %>% 
          mutate(ev = comp1*0 + comp2*2 + comp3*4 + comp4*6 + comp5*8 +
                   comp6*10)
        
        
        
        
        
        
        
                     
