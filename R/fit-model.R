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
#' @param same If TRUE, the Taylor expansion point for the ith observation is
#'   W[ii, ] = u + Xb[ii, ], for an u that does not depend on i.
#' @param Beta A p-vector with starting values for the regression coefficients.
#' @param Sigma An r x r initial iterate for the latent covariance matrix.
#' @param W An n x r initial iterate for the expansion points.
#' @param psi An r-vector of variance parameters.
#' @return A list of final iterates
#' @useDynLib lvmmrPQL, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom Rcpp evalCpp
#' @export
#' @export
lvmmr_PQL <- function(Y, X, type, M, tol = rep(1e-8, 4), maxit = rep(1e2, 4),
                      quiet = rep(TRUE, 4), relative = TRUE, same = FALSE,
                      Beta, Sigma, W, psi = rep(1, ncol(Y)))
{
  # Define constants
  n <- nrow(Y)
  r <- ncol(Y)
  p <- ncol(X)

  # Check that model matrix has full rank
  if(qr(X)$rank != p) warning("X does not have full column rank.")

  # Check that restrictions are symmetric
  if(!isSymmetric(M)){
    warning("M is not symmetric, using lower half only")
    M[upper.tri(M)] <- t(M)[upper.tri(M)]
  }

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
      glm_fit <- stats::glm(y_uni ~ 0 + X_uni, family = fam)
      uni_coefs[, ii] <- stats::coef(glm_fit)
    }
    Beta <- apply(uni_coefs, 1, mean, na.rm = T) # Can also weight by SE
  }

  if(missing(W)) W <- matrix(X %*% Beta, nrow = n, ncol = r, byrow = TRUE)

  if(missing(Sigma)) Sigma <- diag(1e-3, r) # For small var, model approx. GLM

  out_iter <- 0
  iterate_outer <- out_iter < maxit[1] # Iterate updating (Beta, Sigma) and W
  while(iterate_outer){

    # Inner loop with joint update of Beta and Sigma
    in_iter <- 0
    iterate_inner <- in_iter < maxit[2]

    # Avoid using Beta and Sigma in inner loop
    new_Beta <- Beta
    new_Sigma <- Sigma

    # Pre-compute
    D1 <- t(get_cumulant_diffs(W_T = t(W), type = type, order = 1))
    D2 <- t(get_cumulant_diffs(W_T = t(W), type = type, order = 2))
    A <- sweep(D2, 2, psi, FUN = "*")
    while(iterate_inner){
      # Track progress of Beta and Sigma update
      start_obj <- -working_ll(Y = Y, X = X, Beta = new_Beta,
                                 Sigma = new_Sigma, W = W,
                                 psi = psi, D1 = D1, D2 = D2)

      # Update Sigma
      H <- matrix(X %*% new_Beta, nrow = n, ncol = r, byrow = T) # called Xb elsewhere
      H <- D1 + D2 * (H - W)
      H <- Y - H
      new_Sigma <- update_Sigma_trust(Sigma = new_Sigma, R = H, D2 = D2,
                                      psi = psi, M = M, use_idx = 1:n,
                                      pen = 1e-8)
      # new_Sigma <- update_Sigma_PQL(H = H, A =  A, B = D2,
      #                               Sigma.init = new_Sigma,
      #                               M = M, epsilon = 1e-8, tol.dykstra = tol[3],
      #                               tol.ipiano = tol[3],
      #                               max.iter.dykstra = maxit[3],
      #                               max.iter.ipiano = maxit[3],
      #                               quiet = quiet[3])

      if(!quiet[2]){
        mid_obj <-  -working_ll(Y = Y, X = X, Beta = new_Beta,
                                Sigma = new_Sigma, W = W,
                                psi = psi, D1 = D1, D2 = D2)
        cat("Change from Sigma update: ", mid_obj - start_obj, "\n")
        if(mid_obj - start_obj > tol[2]){
          warning("Sigma update increased objective function more than tol[2].
                  \n")
        }
      }



      # Update Beta
      new_Beta <- update_beta(Y = Y, X = X, W = W, Sigma = new_Sigma, psi = psi,
                              type = type)

      # Track progress of Beta and Sigma update
      end_obj <- -working_ll(Y = Y, X = X, Beta = new_Beta,
                               Sigma = new_Sigma, W = W,
                               psi = psi, D1 = D1, D2 = D2)
      if(!quiet[2]){
        cat("Change from Beta update: ", end_obj - mid_obj, "\n")
        if(end_obj - mid_obj > tol[2]){
          warning("Beta update increased objective function more than tol[2].
                  \n")
        }
      }

      # Check whether to terminate inner loop
      in_iter <- in_iter + 1
      change <- abs(end_obj - start_obj)
      if(relative) change <- change / abs(start_obj)
      iterate_inner <- ((change > tol[2]) & (in_iter < maxit[2]))
    } # End inner loop

    # Update W
    if(same){
      W <- update_W_same(Y = Y, X = X, W = W, Beta = new_Beta, Sigma = new_Sigma,
                    psi = psi, type = type, tol = tol[4],
                    maxit = maxit[4], quiet = quiet[4])
    } else{
      W <- update_W(Y = Y, X = X, W = W, Beta = new_Beta, Sigma = new_Sigma,
                    psi = psi, type = type, tol = tol[4],
                    maxit = maxit[4], quiet = quiet[4])
    }

    # Check whether to terminate outer loop
    # (Seems elementwise relative change could be unstable here?)
    if(relative){
      change <- max(abs(c((Sigma - new_Sigma) / max(abs(Sigma)),
                          (Beta - new_Beta) / max(abs(Beta)))))
    } else{
      change <- max(abs(c(Sigma - new_Sigma, Beta - new_Beta)))
    }
    out_iter <- out_iter + 1
    iterate_outer <-  ((change > tol[1]) & (out_iter < maxit[1]))

    # Print progress of outer loop
    if(!quiet[1]){
      cat("Change in parameters: ", change, "\n")
    }

    # Prepare next iteration
    Beta <- new_Beta
    Sigma <- new_Sigma
  } # End outer loop
  return(list(Beta = Beta, Sigma = Sigma, W = W, iter = out_iter,
              change = change))
}


