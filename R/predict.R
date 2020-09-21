#' Predict Latent Variables Mixed-type Multivariate Regression
#'
#' @param X An n x p matrix of predictors (p = p_1 + ... + p_r).
#' @param type An r-vector indicating response types:
#'             1 means normal, 2 means Bernoulli,
#'             and 3 means Poisson.
#' @param n_preds An r-vector whose ith element is the number of predictors for
#' the ith response.
#' @param Beta A p-vector of latent regression coefficients
#' @param sigma An r-vector of standard deviations for the responses,
#'              i.e. diag(solve(Omega)).
#' @param num_nodes Number of nodes to use in Gaussian quadrature used to
#'                  calculate predictions
#' @return An n x r matrix of predicted (fitted) values
#' @export
predict_lvmmr <- function(X, type, n_preds, Beta, sigma, num_nodes = 15)
{
  # Define constants
  n_obs <- nrow(X)
  n_resp <- length(type)
  XBeta <- get_Xb_rcpp(t(X), Beta, n_preds)
  # Standard normal quadrature
  quadrature <- statmod::gauss.quad.prob(n = num_nodes, dist = "normal")
  W0 <- matrix(rep(quadrature$nodes, each = n_obs), nrow = n_obs, byrow = FALSE)
  for(jj in 1:n_resp){
      W <- W0 * sigma[jj]
      W <- sweep(W, 1, STATS = XBeta[, jj], FUN = "+")
      W <- t(get_cumulant_diffs(W_T = t(W), type = rep(type[jj], num_nodes), order = 1))
      W <- sweep(W, 2, STATS = quadrature$weights, FUN = "*")
      XBeta[, jj] <- rowSums(W)
  }
  return(XBeta)
}
