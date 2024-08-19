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
my.cluster <- makeCluster(n.cores, type = "PSOCK")
doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()
foreach::getDoParWorkers()
registerDoMC(cores = n.cores)

gibbmod <- cmdstan_model(stan_file = './time_mix_norm_crps.stan')
gibbmod <- cmdstan_model(stan_file = './time_mix_norm_crps_lambda.stan')
gibbmod <- cmdstan_model(stan_file = './simple_mix_norm_crps_T.stan')


true_wt <- c(.1, .1, ., .35)
true_wt <- c(.35, .1, .15, .4)
true_wt <- c(.2, .8)
log_tr_wt <- log(true_wt) + 7
log_wt <- log(1/4)
replicates <- 500
# wds <- foreach(rep = 1:replicates,
#                     .packages = c("cmdstanr", "stringr",
#                                   "lubridate", "dplyr", "tmvtnorm")
#                     ,.errorhandling = "remove"
#                     ,.combine = rbind) %dopar% {

T <- 100
M <- length(true_wt)
t <- 1:T
zeta <- matrix(NA, nrow = T, ncol = M)
Sigma <- diag(rep(.1,M))
# Sigma <- matrix(.99, nrow = M, ncol = M)
# Sigma <- Sigma + diag(rep(1, M))

zeta[1,] <- MASS::mvrnorm(1, log_tr_wt, Sigma = Sigma)
# zeta[1,] <- rtmvnorm(1, rep(.1, M), sigma = Sigma,
#           lower = rep(0,M))
for (i in 2:T) {
  zeta[i,] <- MASS::mvrnorm(1, zeta[i-1,] + .02, Sigma)
  # zeta[i,] <- rtmvnorm(1, zeta[i-1,], sigma = Sigma,
  #                       lower = rep(0,M))
}

omega <- matrix(NA, nrow = N, ncol = M)
for (i in 1:N) {
  omega[i,] <- exp(zeta[i,])/sum(exp(zeta[i,]))
  # omega[i,] <- zeta[i,]/sum(zeta[i,])
}
# z <- c()
# z[1] <- log(1/4)
# for (i in 2:T) {
#   z[i] <- rnorm(1, z[i-1], 1)
# }
# omega <- matrix(NA, ncol = 2, nrow = T)
# omega[,1] <- exp(z)/(1 + exp(z))
# omega[,2] <- 1 - omega[,1]
# log(true_wt) + 7


plot(omega[,1] ~ t, type = "l", ylim = c(0,1))
lines(omega[,2] ~ t, col = "red")
lines(omega[,3] ~ t, col = "dodgerblue")
lines(omega[,4] ~ t, col = "purple")

# a <- 1:7
a <- c(1,3,5,7,9,11,13)
MS <- length(a)
X <- MASS::mvrnorm(N, a, diag(rep(1,MS)))

# Xl <- apply(X, MARGIN = 2, FUN = function(x) {x*a})
Xl <- X
yt <- c()
# for (i in 1:nrow(Xl)) {
#   yt[i] <- sample(Xl[i,], 1, prob = omega[i,]) #+ rnorm(1)
# }
sigma <- rep(1, 4)
# yt[1] <- rnorm(1,4)
for (i in 1:nrow(Xl)) {
  # mix <- UnivarMixingDistribution(Norm(a[1], sigma[1]),
  #                                 Norm(a[2], sigma[2]),
  #                                 Norm(a[3], sigma[3]),
  #                                 Norm(a[4], sigma[4]),
  #                                 mixCoeff = omega[i,])
  
  mean <- sample(c(4,8), 1, prob = omega[i,])
  
  yt[i] <- rnorm(1, mean, 1)
}
# yt <- rnorm(50, 2)
# yt <- c(yt, rnorm(50,6))
# yt <- sample(yt, 100, replace = FALSE)
# yt <- yt + rnorm(N, 0, .4)
# yt <- apply(Xl*omega, MARGIN = 1, FUN = sum) + rnorm(N,0,.5)
plot(yt ~ t, type = "l")
# plot(yt ~ X[,4])
# points(yt ~ X[,2], col = "red")
# points(yt ~ X[,3], col = "dodgerblue")
# points(yt ~ X[,4], col = "purple")
#p1 <- c()
#p2 <- c()
#p3 <- c()
#p4 <- c()
#for (i in 1:N) {
#  p1[i] <- pnorm(yt[i], X[i,1])
#  p2[i] <- pnorm(yt[i], X[i,2])
#  p3[i] <- pnorm(yt[i], X[i,3])
#  p4[i] <- pnorm(yt[i], X[i,4])
#}
#hist(c(p1, p2, p3, p4))


all_alphas <- array(NA, dim = c(MS,MS,N))
all_betas <- matrix(NA, nrow = N, ncol = MS)
# sigma <- rep(1, M)
sigma <- rep(1, MS)
for (i in 1:N) {
  
  # mean_diff <- outer(Xl[i,], Xl[i,], "-")
  mean_diff <- outer(a, a, "-")
  var_sum <- outer(sigma^2, sigma^2, "+")
  param_func <- array(c(mean_diff, var_sum), dim = c(MS, MS, 2))
  
  all_alphas[,,i] <- apply(param_func, MARGIN = c(1,2), FUN = alphaik)
  # all_betas[i,] <- sapply(yt[i], FUN = betai, mu = Xl[i,], sigma = sigma)
  all_betas[i,] <- sapply(yt[i], FUN = betai, mu = a, sigma = sigma)
  
}


pmixnorm <- function(y, mus, sigmas, wts) {
  return(sum(wts*pnorm(y,mus,sigmas)))
}



alpha <- .98
et <- .1
lambda <- 3

ps <- c()
pm <- c()
pb <- c()
cs <- c()
cm <- c()
cb <- c()
ls <- c()
lm <- c()
lb <- c()
uwds <- matrix(NA, nrow = nrow(all_betas), ncol = 3)
wts <- data.frame()
for (d in 40:(nrow(all_betas) - 1)) {
  # et <-1
  
    stan_dat <- list(
      T = d,
      y = yt[1:d],
      num_comp = nrow(all_alphas),
      eta = et,
      alphas = aperm(all_alphas[,,1:d]),
      betas = all_betas[1:d,],
      # alpha = rep(1,MS),
      alpha = rep(1, MS),
      sigma_s = rep(1, MS - 1),
      lambda = lambda,
      wts = rep(1/length(a), length(a)),
      tweight = .99
    )
    
  
    
    fit <- gibbmod$variational(data = stan_dat)
    # fit <- gibbmod$sample(data = stan_dat,
    #                       chains = 1,
    #                       iter_warmup = 2000,
    #                       iter_sampling = 2000)
    
    # draws <- fit$draws(variables = "omegaT1", format = "df") %>%
    #   select(contains("omega"))
    
    draws <- fit$draws(variables = "omega", format = "df") %>%
      select(contains("omega"))
    
    # rps <- c()
    # for (i in 1:nrow(draws)) {
    #   wt <- draws[i,]
    #   wt <- wt/sum(wt)
    #   rps[i] <- pmixnorm(yt[d + 1], Xl[d + 1,], sigma, wt)
    # }
    
    wt <- apply(draws, MARGIN = 2, FUN = mean)
    wt <- wt/sum(wt)
    wts <- rbind(wts, data.frame(t(wt)))
    # wt <- c(omega[d,], 0, 0)
    # wt <- c(true_wt,0,0)
    
    
    
    
    uwts <- c()
    for (m in 1:length(a)) {
      uwts[m] <- (1/length(a))*exp(-et*d*mean(
        alpha^(d:2 - 1)*
          scoringRules::crps(yt[2:d],
                                     family = "norm", mean = a[m], sd = 1)))*
        prod(dnorm(yt[2:d], a[m]))
    }
    # uwt <- uwts/sum(uwts)
    bma <- c(mixtools::dmvnorm(yt[2:d], a[1]),#Xl[2:d,1]),
             mixtools::dmvnorm(yt[2:d], a[2]),#Xl[2:d,2]),
             mixtools::dmvnorm(yt[2:d], a[3]),#Xl[2:d,3]),
             mixtools::dmvnorm(yt[2:d], a[4]),#Xl[2:d,4]),
             mixtools::dmvnorm(yt[2:d], a[5]),#Xl[2:d,5]),
             mixtools::dmvnorm(yt[2:d], a[6]),
             mixtools::dmvnorm(yt[2:d], a[7]))#,
             # mixtools::dmvnorm(yt[2:d], a[8]),
             # mixtools::dmvnorm(yt[2:d], a[9]),
             # mixtools::dmvnorm(yt[2:d], a[10]))#Xl[2:d,6]))
    bwt <- bma/sum(bma)
    # bwt <- rep(1/7,7)
    mwt <- rep(1/MS,MS)
    
    ps[d] <- pmixnorm(yt[d + 1], a, sigma, wt)
    # ps[d] <- mean(rps, na.rm = TRUE)
    pm[d] <- pmixnorm(yt[d + 1], a, sigma, rep(1/MS, MS))
    pb[d] <- pmixnorm(yt[d + 1], a, sigma, bwt)
    
    cs[d] <- scoringRules::crps(yt[d + 1], family = "mixnorm",
                       m = matrix(a, nrow = 1),
                       s = matrix(sigma, nrow = 1),
                       w = matrix(wt, nrow = 1))
    # 
    cm[d] <- scoringRules::crps(yt[d + 1], family = "mixnorm",
                                m = matrix(a, nrow = 1),
                                s = matrix(sigma, nrow = 1),
                                w = matrix(mwt, nrow = 1))
    # 
    cb[d] <- scoringRules::crps(yt[d + 1], family = "mixnorm",
                                m = matrix(a, nrow = 1),
                                s = matrix(sigma, nrow = 1),
                                w = matrix(bwt, nrow = 1))
  
  
    ls[d] <- scoringRules::logs(yt[d + 1], family = "mixnorm",
                                m = matrix(a, nrow = 1),
                                s = matrix(sigma, nrow = 1),
                                w = matrix(wt, nrow = 1))
  
    lm[d] <- scoringRules::logs(yt[d + 1], family = "mixnorm",
                                m = matrix(a, nrow = 1),
                                s = matrix(sigma, nrow = 1),
                                w = matrix(mwt, nrow = 1))
  
    lb[d] <- scoringRules::logs(yt[d + 1], family = "mixnorm",
                                m = matrix(a, nrow = 1),
                                s = matrix(sigma, nrow = 1),
                                w = matrix(bwt, nrow = 1))
    
    et <- sum(ls[d], na.rm = TRUE)/sum(cs[d], na.rm = TRUE)
    uwds[d,] <- c(unit_wass_dist(ecdf(pm)),
              unit_wass_dist(ecdf(pb)),
              unit_wass_dist(ecdf(ps)))
    
    print_color(d, color = "blue")
    print("")
    # print_color(c(mean(cm, na.rm = TRUE),mean(cb, na.rm = TRUE),mean(cs, na.rm = TRUE)), 
    #             color = "green")
    
    print_color(c(cm[d], cb[d], cs[d]), 
                color = "green")
    print("")
    print_color(uwds[d,], 
                color = "green")
    
  
  
}

uwd1 <- c(unit_wass_dist(ecdf(pm)),
          unit_wass_dist(ecdf(pb)),
          unit_wass_dist(ecdf(ps)))

uwd2 <- c(unit_wass_dist(ecdf(pm), d = 2),
          unit_wass_dist(ecdf(pb), d = 2),
          unit_wass_dist(ecdf(ps), d = 2))

ks_stat <- c(ks.test(pm, y = "punif")$statistic,
             ks.test(pb, y = "punif")$statistic,
             ks.test(ps, y = "punif")$statistic)

logs <- c(mean(lm, na.rm = TRUE),
          mean(lb, na.rm = TRUE),
          mean(ls, na.rm = TRUE))

crps <- c(mean(cm, na.rm = TRUE),
          mean(cb, na.rm = TRUE),
          mean(cs, na.rm = TRUE))
uwd1/min(uwd1); uwd2/min(uwd2); ks_stat/min(ks_stat); logs/min(logs); crps/min(crps)

 data.frame(rep = rep, method = c("ms", "bma", "bs"), uwd1, uwd2, ks_stat, 
           logs, crps)

# }
# 
# saveRDS(wds, "wds.rds")


# M <- 1000
# m <- 1
# diff <- c()
# repeat{
#   diff[m] <- unit_wass_dist(ecdf(sample(pm,95,replace = FALSE)), d = 2)/
#                 unit_wass_dist(ecdf(sample(ps,95,replace = FALSE)), d = 2)
# 
#   m <- m + 1
#   if (m > M) {break}
# }
# 
# unit_wass_dist(ecdf(pm[1:100]), d = 1)/
#   unit_wass_dist(ecdf(ps[1:100]), d = 1)
# 
# ks.test(pm[1:100], y = "punif")$statistic/ks.test(ps[1:100], y = "punif")$statistic

















