#' Calculate predictions (marginal expectations)
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
predict_mmrr <- function(X, Beta, sigma, type, num_nodes = 15)
{
  # Argument checking
  stopifnot(is.matrix(X),
            is.numeric(Beta), is.atomic(Beta),
            is.numeric(sigma), is.atomic(sigma), all(sigma >= 0),
            is.numeric(type), is.atomic(type), all(type %in% 1:3),
            is.numeric(num_nodes), is.atomic(num_nodes),
            length(num_nodes) == 1, floor(num_nodes) == num_nodes,
            num_nodes >= 1)

  # Define constants
  p <- ncol(X)
  r <- length(type)
  n <- nrow(X) / r
  stopifnot(floor(n) == n, length(Beta) == p, length(sigma) == r)

  Xb <- matrix(X %*% Beta, nrow = n, ncol = r, byrow = T)
  # Standard normal quadrature
  grid_gauss <- mvQuad::createNIGrid(dim = 1, type = "GHN", level = num_nodes,
                                     level.trans = FALSE)
  nodes <- as.vector(mvQuad::getNodes(grid_gauss))
  weights <- as.vector(mvQuad::getWeights(grid_gauss))
  W0 <- matrix(rep(nodes, each = n), nrow = n, ncol = num_nodes,
               byrow = FALSE)
  for(jj in 1:r){
      if(type[jj] == 1){
        # Prediction is linear predictor
      } else if (type[jj] == 3){
        Xb[, jj] <- exp(Xb[, jj] + 0.5 * sigma[jj]^2)
      } else{
        W <- W0 * sigma[jj]
        W <- sweep(W, 1, STATS = Xb[, jj], FUN = "+")
        W <- t(get_cumulant_diffs(W_T = t(W), type = rep(type[jj], num_nodes),
                                  order = 1))
        W <- sweep(W, 2, STATS = weights, FUN = "*")
        Xb[, jj] <- rowSums(W)
      }
  }
  return(Xb)
}

#' Calculate marginal covariance matrix of responses given predictors
#'
#' @param X An r x p matrix of predictors.
#' @param Beta A p-vector of latent regression coefficients
#' @param Sigma An r x r covariance matrix for latent vector
#' @param psi An r-vector of conditional variance parameters
#' @param type An r-vector indicating response types:
#'             1 means normal, 2 means Bernoulli,
#'             and 3 means Poisson.
#' @param num_nodes Number of nodes for Gaussian quadrature used to
#'                  calculate predictions
#' @return An r x r covariance matrix for responses given predictors (Y|X)
#' @export
cov_mmrr <- function(X, Beta, Sigma, psi, type, num_nodes = 10)
{
  # Argument checking
  stopifnot(is.matrix(X),
            is.numeric(Beta), is.atomic(Beta),
            is.matrix(Sigma), is.numeric(Sigma),
            is.numeric(psi), is.atomic(psi), all(psi > 0),
            is.numeric(type), is.atomic(type), all(type %in% 1:3),
            is.numeric(num_nodes), is.atomic(num_nodes),
            length(num_nodes) == 1, floor(num_nodes) == num_nodes,
            num_nodes >= 1)
  # Define constants
  p <- ncol(X)
  r <- length(type)
  stopifnot(length(Beta) == p, all(dim(Sigma) == r), length(psi) == r)


  Xb <- X %*% Beta
  mu <- predict_mmrr(X = X, Beta = Beta, sigma = sqrt(diag(Sigma)), type = type,
                    num_nodes = num_nodes)
  grid_gauss <- mvQuad::createNIGrid(dim = 2, type = "GHN", level = num_nodes,
                                     level.trans = FALSE)
  cov_mat <- matrix(0, r, r)
  for(jj in 1:r){
    for(kk in 1:jj){
      if(jj != kk){
        if((type[jj] == 1) & type[kk] == 1){
         cov_mat[jj, kk] <- Sigma[jj, kk]
       } else if((type[jj] == 1) & (type[kk] == 3)){
         cov_mat[jj, kk] <- Sigma[jj, kk] * exp(Xb[kk] + 0.5 * Sigma[kk, kk])
       } else if((type[jj] == 3) & (type[kk] == 1)){
         cov_mat[jj, kk] <- Sigma[kk, jj] * exp(Xb[jj] + 0.5 * Sigma[jj, jj])
       } else if ((type[jj] == 3) & (type[kk] == 3)){
         cov_mat[jj, kk] <- exp(Xb[jj] + Xb[kk] + 0.5 * Sigma[jj, jj] +
                                  0.5 *Sigma[kk, kk] + Sigma[jj, kk]) -
           exp(Xb[jj] + Xb[kk] + 0.5 * Sigma[jj, jj] + 0.5 * Sigma[kk, kk])
       } else{
         R <- chol(Sigma[c(jj, kk), c(jj, kk)])
         integrand <- function(w){
           w <- w %*% R + matrix(Xb[c(jj, kk)], nrow = nrow(w), ncol = 2, byrow = TRUE)
           w <- t(get_cumulant_diffs(W_T = t(w), type = type[c(jj, kk)], order = 1))
           w <- w - matrix(mu[c(jj, kk)], nrow = nrow(w), ncol = 2, byrow = TRUE)
           w[, 1] * w[, 2]
         }
         cov_mat[jj, kk] <- mvQuad::quadrature(integrand, grid = grid_gauss)
       }
       cov_mat[kk, jj] <- cov_mat[jj, kk]
      } else if(type[jj] == 1){ # Normal variance
       cov_mat[jj, jj] <- Sigma[jj, jj] + psi[jj]
     } else if(type[jj] == 2){ # Bernoulli variance
       cov_mat[jj, jj] <- mu[jj] - mu[jj]^2
     } else{ # Poisson variance
       cov_mat[jj, jj] <- psi[jj] * exp(Xb[jj] + 0.5 * Sigma[jj, jj]) +
         exp(2 * Xb[jj] + 2 * Sigma[jj, jj]) - mu[jj]^2
     }
    }
  }
  return(cov_mat)
}
