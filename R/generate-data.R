#' Generate data from mixed-type multivariate response regression model
#'
#' @param X An nr x p matrix of predictors. The first r elements are predictors
#'   for the r responses in the first independent response vector, and so on.
#' @param Beta A p-vector of latent regression coefficients
#' @param R An r x r square root of the covariance matrix of the latent vectors.
#'   Should be such that t(R) \%*\% R = solve(Omega).
#' @param type An r-vector indicating response types: 1 means normal, 2 means
#'   Bernoulli, and 3 means Poisson.
#' @param psi Vector of variance parameters, ignored for Bernoulli
#'   responses. For Poisson, generates Y| W as psi * Poi(exp(W) / psi)
#' @return An n x r matrix of responses.
#' @export
generate_mmrr <- function(X, Beta, R, type, psi)
{
  # R should be such that t(R) %*% R = Sigma = solve(Omega)
  # Beta and X should be supplied as matrices
  p <- ncol(X)
  r <- nrow(R)
  n <- nrow(X) / r
  W <- matrix(stats::rnorm(r * n), n, r) %*% R
  W <- W + matrix(X %*% Beta, nrow = n, ncol = r, byrow = T)

  # Replace W by Y, using same storage
  for(ii in 1:r){
    if(type[ii] == 1){
      W[, ii] <- stats::rnorm(n, mean = W[, ii], sd = sqrt(psi[ii]))
    } else if(type[ii] == 2){
      W[, ii] <- stats::rbinom(n, 1, 1 / (1 + exp(-W[, ii])))
    } else if(type[ii] == 3){
      W[, ii] <- psi[ii] * stats::rpois(n, exp(W[, ii]) / psi[ii])
    } else{
      stop("type must be 1 (normal) , 2 (Bernoulli), or 3 (Poisson).")
    }
  }
  # Return responses
  return(W)
}
