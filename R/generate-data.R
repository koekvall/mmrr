#' Generate Latent Variables Mixed-type Multivariate Regression
#'
#' @param X An n x p matrix of predictors (p = p_1 + ... + p_r).
#' @param type An r-vector indicating response types:
#'             1 means normal, 2 means Bernoulli,
#'             and 3 means Poisson.
#' @param n_preds An r-vector which ith element is the number of predictors for
#' the ith response.
#' @param Beta A p-vector of latent regression coefficients
#' @param R An r x r square root of the covariance matrix of the latent vectors.
#'          Should be such that t(R) \%*\% R = solve(Omega).
#' @return An n x r matrix of responses.
#' @export
generate_lvmmr <- function(X, type, n_preds, Beta, R)
{
  # R should be such that t(R) %*% R = Sigma = solve(Omega)
  # Beta and X should be supplied as matrices
  n_obs <- nrow(X)
  n_pred <- ncol(X)
  n_resp <- ncol(R)
  W <- matrix(stats::rnorm(n_resp * n_obs), n_obs, n_resp) %*% R
  W <- W + get_Xb_rcpp(t(X), Beta, n_preds)

  # Replace W by Y, using same storage
  for(ii in 1:n_resp){
    if(type[ii] == 1){
      W[, ii] <- stats::rnorm(n_obs, mean = W[, ii], sd = 1)
    } else if(type[ii] == 2){
      W[, ii] <- stats::rbinom(n_obs, 1, 1 / (1 + exp(-W[, ii])))
    } else if(type[ii] == 3){
      W[, ii] <- stats::rpois(n_obs, exp(W[, ii]))
    } else{
      stop("type must be 1 (normal) , 2 (Bernoulli), or 3 (Poisson).")
    }
  }
  # Return responses
  return(W)
}
