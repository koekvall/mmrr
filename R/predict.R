#' Predict Latent Variables Mixed-type Multivariate Regression
#'
#' @param X An nr x p matrix of predictors. The first r elements are predictors
#'   for the r responses in the first independent response vector, and so on.
#' @param Beta A p-vector of latent regression coefficients
#' @param sigma An r-vector of standard deviations for the responses,
#'              i.e. sqrt(diag(Sigma)).
#' @param type An r-vector indicating response types:
#'             1 means normal, 2 means Bernoulli,
#'             and 3 means Poisson.
#' @param num_nodes Number of nodes to use in Gaussian quadrature used to
#'                  calculate predictions
#' @return An n x r matrix of predicted (fitted) values
#' @export
predict_lvmmr <- function(X, Beta, sigma, type, num_nodes = 15)
{
  # Define constants
  p <- ncol(X)
  r <- length(type)
  n <- nrow(X) / r

  Xb <- matrix(X %*% Beta, nrow = n, ncol = r, byrow = T)
  # Standard normal quadrature
  quadrature <- statmod::gauss.quad.prob(n = num_nodes, dist = "normal")
  W0 <- matrix(rep(quadrature$nodes, each = n), nrow = n, ncol = num_nodes,
               byrow = FALSE)
  for(jj in 1:r){
      W <- W0 * sigma[jj]
      W <- sweep(W, 1, STATS = Xb[, jj], FUN = "+")
      W <- t(get_cumulant_diffs(W_T = t(W), type = rep(type[jj], num_nodes), order = 1))
      W <- sweep(W, 2, STATS = quadrature$weights, FUN = "*")
      Xb[, jj] <- rowSums(W)
  }
  return(Xb)
}
