# library(TruncatedNormal)
library(tmvtnorm)
library(cmdstanr)
library(dplyr)

gibbmod <- cmdstan_model(stan_file = './time_mix_norm_crps.stan')

N <- 100
M <- 4
t <- 1:N
zeta <- matrix(NA, nrow = N, ncol = M)
Sigma <- diag(rep(1,M))
# Sigma <- matrix(.99, nrow = M, ncol = M)
# Sigma <- Sigma + diag(rep(1, M))

zeta[1,] <- MASS::mvrnorm(1, rep(0,M), Sigma = Sigma)
# zeta[1,] <- rtmvnorm(1, rep(.1, M), sigma = Sigma,
#          lower = rep(0,M))
for (i in 2:N) {
  zeta[i,] <- MASS::mvrnorm(1, zeta[i-1,] + .02, Sigma)
  # zeta[i,] <- rtmvnorm(1, zeta[i-1,], sigma = Sigma,
  #                      lower = rep(0,M))
}

omega <- matrix(NA, nrow = N, ncol = M)
for (i in 1:N) {
  omega[i,] <- exp(zeta[i,])/sum(exp(zeta[i,]))
  # omega[i,] <- zeta[i,]/sum(zeta[i,])
}

plot(omega[,1] ~ t, type = "l")
lines(omega[,2] ~ t, col = "red")
lines(omega[,3] ~ t, col = "dodgerblue")
lines(omega[,4] ~ t, col = "purple")

a <- c(1.1,1.2,1.3,1.4)
X <- MASS::mvrnorm(N, rep(5,M), diag(rep(.5,M)))

# Xl <- apply(X, MARGIN = 2, FUN = function(x) {x*a})
Xl <- X

yt <- apply(Xl*omega, MARGIN = 1, FUN = sum) + rnorm(N)
plot(yt ~ t, type = "l")
# plot(yt ~ X[,4])
# points(yt ~ X[,2], col = "red")
# points(yt ~ X[,3], col = "dodgerblue")
# points(yt ~ X[,4], col = "purple")
p1 <- c()
p2 <- c()
p3 <- c()
p4 <- c()
for (i in 1:N) {
  p1[i] <- pnorm(yt[i], X[i,1])
  p2[i] <- pnorm(yt[i], X[i,2])
  p3[i] <- pnorm(yt[i], X[i,3])
  p4[i] <- pnorm(yt[i], X[i,4])
}
hist(c(p1, p2, p3, p4))


all_alphas <- array(NA, dim = c(M,M,N))
all_betas <- matrix(NA, nrow = N, ncol = M)
sigma <- rep(1, M)

for (i in 1:N) {
  
  mean_diff <- outer(Xl[i,], Xl[i,], "-")
  var_sum <- outer(sigma^2, sigma^2, "+")
  param_func <- array(c(mean_diff, var_sum), dim = c(M, M, 2))
  
  all_alphas[,,i] <- apply(param_func, MARGIN = c(1,2), FUN = alphaik)
  all_betas[i,] <- sapply(yt[i], FUN = betai, mu = Xl[i,], sigma = sigma)
  
}


pmixnorm <- function(y, mus, sigmas, wts) {
  return(sum(wts*pnorm(y,mus,sigmas)))
}




et <- 25
ps <- c()
pm <- c()
pb <- c()
for (d in 2:(nrow(all_betas) - 1)) {
  
  stan_dat <- list(
    T = d,
    y = yt[1:d],
    num_comp = nrow(all_alphas),
    eta = et,
    alphas = aperm(all_alphas[,,1:d]),
    betas = all_betas[1:d,]
  )
  
  fit <- gibbmod$variational(data = stan_dat)
  
  draws <- fit$draws(variables = "omegaT1", format = "df") %>% 
    select(contains("omega"))
  
  wt <- apply(draws, MARGIN = 2, FUN = mean)
  wt <- wt/sum(wt)
  
  bma <- c(mixtools::dmvnorm(yt[2:d], Xl[2:d,1]),
           mixtools::dmvnorm(yt[2:d], Xl[2:d,2]),
           mixtools::dmvnorm(yt[2:d], Xl[2:d,3]),
           mixtools::dmvnorm(yt[2:d], Xl[2:d,4]))
  bwt <- bma/sum(bma)
  
  
  
  ps[d] <- pmixnorm(yt[d + 1], Xl[d,], sigma, wt)
  pm[d] <- pmixnorm(yt[d + 1], Xl[d,], sigma, rep(1/M, M))
  pb[d] <- pmixnorm(yt[d + 1], Xl[d,], sigma, bwt)
  
  
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

















