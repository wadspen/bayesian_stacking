
vMv <- function(mat, wt) {
		return(wt%*%mat%*%wt)
}

mix_mat_crps <- function(wt, mse_mat, absdiff_arr) {
			#wt <- exp(wt)/sum(exp(wt))
			if (is.vector(mse_mat)) {
				one <- sum(wt*mse_mat)
				two <- wt%*%absdiff_arr%*%wt
			} else {
				one <- apply(mse_mat, MARGIN = 2, FUN = function(x) {sum(wt*x)})
				two <- apply(absdiff_arr, MARGIN = 3, FUN = vMv, wt = wt)
			}
			
			crpss <- one - (1/2)*two
			return(crpss)
}

sum_crps <- function(wt, mse_mat, absdiff_arr) {
	wt <- exp(wt)/sum(exp(wt))
	crpss <- mix_mat_crps(wt, mse_mat, absdiff_arr)
	return(sum(crpss))
}

sum_weight_crps <- function(wt, mse_mat, absdiff_arr, alpha = .98, T) {
  wt <- exp(wt)/sum(exp(wt))
  t <- 1:T
  crpss <- (alpha^(T - t))*mix_mat_crps(wt, mse_mat, absdiff_arr)
  return(sum(crpss))
}


alphaik <- function(par_fun) {
  mu_diff <- par_fun[1]
  sig2_sum <- par_fun[2]
  # wt_prod <- par_fun[3]
  alpha <- -(1/sqrt(2*pi))*sqrt(sig2_sum)*exp(-(mu_diff^2/(2*sig2_sum))) -
    ((mu_diff/2)*(2*pnorm(mu_diff/sqrt(sig2_sum)) - 1))
  return(alpha)
}

# betai <- function(y, mu, sigma, wt) {
#   beta <- sqrt(2/pi)*sigma*exp(-((mu - y)^2/(2*sigma^2))) +
#     (mu - y)*(2*pnorm(y, mu, sigma) - 1)
#   return(beta*wt)
# }

betai <- function(y, mu, sigma) {
  beta <- sqrt(2/pi)*sigma*exp(-(mu - y)^2/(2*sigma^2)) +
    ((mu - y))*(2*pnorm((mu - y)/sigma) - 1)
  return(beta)
}


# alphaikw <- function(par_fun, wt) {
#   mu_diff <- par_fun[1]
#   sig2_sum <- par_fun[2]
#   # wt_prod <- par_fun[3]
#   alpha <- -(1/sqrt(2*pi))*sqrt(sig2_sum)*exp(-(mu_diff^2/(2*sig2_sum))) -
#     ((mu_diff/2)*(2*pnorm(mu_diff/sqrt(sig2_sum)) - 1))
#   return(alpha*wt)
# }
# 
# betaiw <- function(y, mu, sigma, wt) {
#   beta <- sqrt(2/pi)*sigma*exp(-(mu - y)^2/(2*sigma^2)) +
#     ((mu - y))*(2*pnorm((mu - y)/sigma) - 1)
#   return(beta*wt)
# }


CRPS <- function(y, mu, sigma, wt) {
  mean_diff <- outer(mu, mu, "-")
  var_sum <- outer(sigma^2, sigma^2, "+")
  param_func <- array(c(mean_diff, var_sum), dim = c(M, M, 2))
  alpha <- apply(param_func, MARGIN = c(1,2), FUN = alphaik)
  beta <- betai(y, mu, sigma)
  
  return(sum(beta*wt) + sum(wt%*%alpha%*%wt))
}

abs_punit <- function(x, ppitd_est, d = 1) {
  abs(ppitd_est(x) - x)^d
}

unit_wass_dist <- function(ppitd_est, d = 1) {
  (d + 1)*integrate(abs_punit, lower = 0, upper = 1, ppitd_est, d = d,
                    subdivisions = 3000)$value
}
