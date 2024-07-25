
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
