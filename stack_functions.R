
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
