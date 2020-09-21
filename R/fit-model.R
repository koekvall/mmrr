#' Fit Latent Variables Mixed-type Multivariate Regression by PQL
#'
#' @param Y An n x r matrix of responses.
#' @param X An nr x p matrix of predictors.
#' @param type An r-vector indicating response types: 1 means normal, 2 means
#'   Bernoulli, and 3 means Poisson.
#' @param M An r x r matrix with restrictions, with NA for unrestricted
#' @param tol A 4-vector with tolerances for termination of: [1] overall
#'    algorithm, [2] update of Beta and Sigma with W fixed, [3] update
#'    of Beta, and [4] update of Sigma.
#'@param maxit A 4-vector with maximum number of iterations for the same steps
#'   as the tol vector.
#' @param quiet  A 4-vector indicating whether to print information for the
#'   same steps as the tol vector.
#' @param relative If TRUE, use relative decrease of parameters to determine
#'   convergence, otherwise use absolute.
#' @param Beta A p-vector with starting values for the regression coefficients.
#' @param W An n x r initial iterate for the expansion points.
#' @param psi An r-vector of variance parameters.
#' @return A list of final iterates
#' @useDynLib lvmmr, .registration = TRUE
#' @export
lvmmr_PQL <- function(Y, X, type, M, tol = rep(1e-8, 4), maxit = rep(1e2, 4),
                      quiet = rep(TRUE, 4), relative = TRUE, Beta, W,
                      psi = rep(1, ncol(Y)))
{
  # Define constants
  n <- nrow(Y)
  r <- ncol(Y)
  p <- ncol(X)

  # Check that model matrix has full rank
  if(qr(X)$rank != p) warning("X does not have full column rank.")

  # Get starting values from GLMs if missing
  if(missing(Beta)){
    # Fit separately for every response type and use
    # average coefficient as starting value
    unique_types <- unique(type)
    uni_coefs <- matrix(0, ncol = length(unique_types), nrow = p)
    for(ii in 1:length(unique_types)){
      fam <- c("gaussian", "binomial", "poisson")[unique_types[ii]]
      y_uni <- c(Y[, type == unique_types[ii]])
      X_uni <- X[rep(type == unique_types[ii], n), ]
      glm_fit <- glm(y_uni ~ 0 + X_uni, family = fam)
      uni_coef[, ii] <- coef(glm_fit)
    }
    Beta <- apply(uni_coefs, 1, mean, na.rm = T) # Can also weight by SE
  }

  if(missing(W)) W <- matrix(X %*% Beta, nrow = n, ncol = r, byrow = TRUE)

  out_iter <- 0
  iterate_outer <- out_iter < maxit[1] # Iterate updating (Beta, Sigma) and W
  while(iterate_outer){
    in_iter <- 0
    iterate_inner <- in_iter < maxit[2] # Iterate updating Beta and Sigma
    while(iterate_inner){
      D1 <- t(get_cumulant_diffs(W_T = t(W), type = type, order = 1))
      D2 <- t(get_cumulant_diffs(W_T = t(W), type = type, order = 2))
      start_obj <- -working_ll(Y, X, Beta, Sigma, W, psi, D1, D2)
    }


    new_Sigma <- update_Sigma_PQL(H = H, A = Psi * B, B = B, Sigma.init = Sigma,
                                  M = M, epsilon = 1e-4, tol.dykstra = tol_upd[1],
                                  tol.ipiano = tol_upd[1],
                                  max.iter.dykstra = maxit_upd[1],
                                  max.iter.ipiano = maxit_upd[1],
                                  quiet = quiet_upd[1])

    new_Beta <- update_beta_pql_rcpp(Sigma = new_Sigma, Psi = Psi, Y_T = t(Y),
                                     X_T = t(X), W_T = t(W), type = type,
                                     n_preds = n_preds, mu_1 = 0)

    beta_diff <- sum(abs(Beta - new_Beta))
    sigma_diff <- sum(abs(Sigma - new_Sigma))
    if(relative){
      beta_diff <- beta_diff / sum(abs(Beta))
      sigma_diff <- sigma_diff / sum(abs(Sigma))
    }

    Beta <- new_Beta
    Sigma <- new_Sigma

    W <- update_W_pql(Beta = Beta, Omega = qr.solve(Sigma + diag(tol, n_resp)), Y = Y, X = X, W = W,
                      type = type, n_preds = n_preds, pen = 1, quiet = quiet_upd[2])

    difference <- max(beta_diff, sigma_diff)
    iter <- iter + 1

    if(!quiet){
      cat("Difference: ", difference, " iter: ", iter, "\n")
    }

    if(iter >= maxit){
      warning("Reached maximum number of iterations")
      break
    }
  }
  return(list(Beta = Beta, Sigma = Sigma, W = W, iter  = iter, diff = difference))
}
