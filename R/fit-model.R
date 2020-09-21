#' Fit Latent Variables Mixed-type Multivariate Regression by PQL
#'
#' @param Y An n x r matrix of responses.
#' @param X An nr x p matrix of predictors.
#' @param type An r-vector indicating response types: 1 means normal, 2 means
#'   Bernoulli, and 3 means Poisson.
#' @param M An r x r matrix with restrictions, with NA for unrestricted
#' @param tol A 4-vector with tolerances for termination of: [1] overall
#'    algorithm, [2] update of Beta and Sigma with W fixed, [3] update
#'    of Sigma, and [4] update of W.
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

  if(missing(Sigma)) Sigma <- diag(1e-6, r) # For small var, model approx. GLM

  out_iter <- 0
  iterate_outer <- out_iter < maxit[1] # Iterate updating (Beta, Sigma) and W
  while(iterate_outer){

    # Joint update of Beta and Sigma
    in_iter <- 0
    iterate_inner <- in_iter < maxit[2]
    while(iterate_inner){
      D1 <- t(get_cumulant_diffs(W_T = t(W), type = type, order = 1))
      D2 <- t(get_cumulant_diffs(W_T = t(W), type = type, order = 2))

      # Keep track of progress of joint Beta and Sigma update step
      start_obj <- -working_ll(Y, X, Beta, Sigma, W, psi, D1, D2)

      # Update Beta
      new_Beta <- update_beta(Y = Y, X = X, W = W, Sigma = Sigma, psi = psi,
                              type = type)

      # Update Sigma
      Xb <- matrix(X %*% Beta, nrow = n, ncol = r, byrow = T)
      new_Sigma <- update_Sigma_PQL(H = Xb, A = psi * D2, B = D2, Sigma.init = Sigma,
                                    M = M, epsilon = 0, tol.dykstra = tol[3],
                                    tol.ipiano = tol[3],
                                    max.iter.dykstra = maxit[3],
                                    max.iter.ipiano = maxit[3],
                                    quiet = quiet[3])

      # Keep track of progress of joint Beta and Sigma update step
      end_obj <- -working_ll(Y, X, new_Beta, new_Sigma, W, psi, D1, D2)
      if(!quiet[2]){
        cat("Change from Beta and Sigma update: ", end_obj - start_obj, "\n")
        if(end_obj - start_obj > 0){
          warning("Joint Beta and Sigma update increased objective function.")
        }
      }

      # Check whether to terminate inner loop
      in_iter <- in_iter + 1
      change <- abs(end_obj - start_obj)
      if(relative) change <- change / abs(start_obj)
      iterate_inner <- ((change > tol[2]) & (in_iter < maxit[2]))
    }

    # Update W
    W <- update_W(Y = Y, X = X, W = W, Beta = new_Beta, Sigma = new_Sigma,
                 psi = psi, type = type, pen = 1e-6, tol = tol[4],
                 maxit = maxit[4], quiet = quiet[4])

    # Check whether to terminate inner loop
    # Seems elementwise relative change could be unstable here?
    if(relative){
      change <- max(abs(c((Sigma - new_Sigma) / max(abs(Sigma)),
                          (Beta - new_Beta) / max(abs(Beta)))))
    } else{
      change <- max(abs(c(Sigma - new_Sigma, Beta - new_Beta)))
    }
    out_iter <- out_iter + 1
    iterate_outer <-  ((change > tol[1]) & (out_iter < maxit[1]))

    # Prepare next iteration
    Beta <- new_Beta
    Sigma <- new_Sigma
  }
  return(list(Beta = Beta, Sigma = Sigma, W = W, iter = out_iter,
              change = change))
}


