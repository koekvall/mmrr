#' Fit Latent Variables Mixed-type Multivariate Regression by PQL
#'
#' @param Y An n x r matrix of responses.
#' @param X An nr x p matrix of predictors or a list of length r whose ith
#'   element is an n x p_i design matrix for the ith response.
#' @param type An r-vector indicating response types: 1 means Normal, 2 means
#'    Bernoulli, and 3 means (quasi-)Poisson.
#' @param psi An r-vector of conditional variance parameters.
#' @param M An r x r matrix with restrictions for Sigma, with NA for unrestricted.
#' @param tol A 4-vector with tolerances for termination of: [1] overall
#'    algorithm, [2] update of Beta and Sigma with W fixed, [3] update
#'    of Sigma, and [4] update of W.
#'@param maxit A 4-vector with maximum number of iterations for the same steps
#'   as the tol vector.
#' @param quiet  A 4-vector indicating whether to print information for the
#'   same steps as the tol vector.
#' @param relative If TRUE, use relative decrease of parameters to determine
#'   convergence, otherwise use absolute.
#' @param pgd If TRUE, use projected gradient descent; ensures SPSD Sigma.
#' @param eps Lower bound for the smallest eigenvalue of Sigma, only used if
#'   pgd = TRUE.
#' @param uni_fit If TRUE, fit r separate models. This requires (i) X is
#'    a list or (ii) X is a matrix and r is a divisor of p. If (ii), it is
#'    assumed that the first p / r columns of X correspond to the first
#'    response, and so on.
#' @param Beta Initial iterate of regression coefficient vector. Either
#'    a p-vector or a list of length r, where each element is the
#'    coefficient vector for the ith response. Is obtained by fitting separate
#'    GLMs if not supplied.
#' @param Sigma An r x r initial iterate for the latent covariance matrix.
#'    Is set to diag(1e-3, ncol(Y)) if not supplied.
#' @param W An n x r initial iterate for the expansion points.
#'    Is set to matrix(X %*% Beta, nrow = n, ncol = r, byrow = TRUE) if not
#'    supplied.
#' @param w_pen Ridge penalty in W update; often useful to avoid overflows.
#'   Defaults to largest eigenvalue of current Sigma iterate if not supplied.
#' @return A list of final iterates and other information about the fit.
#' @useDynLib lvmmr, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom Rcpp evalCpp
#' @export
lvmmr <- function(Y,
                  X,
                  type,
                  psi = rep(1, ncol(Y)),
                  M = matrix(NA, nrow = ncol(Y), ncol = ncol(Y)),
                  tol = rep(1e-8, 4),
                  maxit = rep(1e2, 4),
                  quiet = rep(TRUE, 4),
                  relative = TRUE,
                  pgd = FALSE,
                  eps = 0,
                  uni_fit = FALSE,
                  Beta,
                  Sigma,
                  W,
                  w_pen)
{
  #############################################################################
  # Argument validation and preparations
  #############################################################################
  stopifnot(is.matrix(Y))
  r <- ncol(Y)
  n <- nrow(Y)

  stopifnot(is.list(X) | is.matrix(X))
  if(is.list(X)){
    stopifnot(length(X) == ncol(Y),
              all(sapply(X, is.matrix)),
              all(sapply(X, nrow) == nrow(Y)))
    X_list <- X
    n_pred <- unlist(sapply(X_list, ncol))
    X <- as.matrix(Matrix::bdiag(X_list))
  }
  p <- ncol(X)

  # Create list for univariate fitting if it does not exist
  if(uni_fit & !exists("X_list") & r > 1){
    pii <- p / r
    if(pii != floor(pii)){
      stop("Cannot create list of design matrices for univariate fitting
            because r does not divide p.")
    }
    X_list <- list()
    for(ii in 1:r){
      X_list[[ii]] <- X[seq(ii, n * r, r), seq((ii - 1) * pii + 1,
                                               length.out = pii), drop = FALSE]
    }
    n_pred <- unlist(sapply(X_list, ncol))
  }

  stopifnot(nrow(X) == (n * r),
            is.numeric(type), length(type) == r,
            is.matrix(M), ncol(M) == r, nrow(M) == r,
            is.numeric(tol), length(tol) == 4,
            is.numeric(maxit), length(maxit) == 4,
            is.logical(quiet), length(quiet) == 4,
            is.logical(relative), length(relative) == 1,
            is.logical(pgd), length(pgd) == 1,
            is.numeric(eps), length(eps) == 1
            )

  if(uni_fit){
    M[upper.tri(M)] <-0
    M[lower.tri(M)] <- 0
  }

  if(!missing(Sigma)) stopifnot(is.matrix(Sigma),
                                ncol(Sigma) == r,
                                nrow(Sigma) == r)
  if(!missing(W)) stopifnot(is.matrix(W),
                            ncol(W) == r,
                            nrow(W) == n)
  if(!missing(w_pen)) stopifnot(is.numeric(w_pen),
                                length(w_pen) ==1)

  if(qr(X)$rank != p) warning("X does not have full column rank.")

  if(!isSymmetric(M)){
    warning("M is not symmetric, using lower half only")
    M[upper.tri(M)] <- t(M)[upper.tri(M)]
  }

  if(all(!is.na(M))){
    message("Skipping Sigma update because all elements constrained.")
    maxit[3] <- 0
    Sigma <- M
  }

  test_M <- M
  test_M[is.na(test_M)] <- 1
  diag(test_M) <- 0
  if(any(colSums(abs(test_M)) == 0) & (r != 1)){
    message("Restrictions imply one or more responses are independent of all
            others; consider separate models.")
  }
  rm(test_M)

  # Get starting values from GLMs if missing
  if(missing(Beta)){
    # Fit separately for every response type and use
    # average coefficient as starting value
    unique_types <- unique(type)
    uni_coefs <- matrix(0, ncol = length(unique_types), nrow = p)
    for(ii in 1:length(unique_types)){
      fam <- c("gaussian", "binomial", "quasipoisson")[unique_types[ii]]
      y_uni <- c(Y[, type == unique_types[ii]])
      X_uni <- X[rep(type == unique_types[ii], n), ]
      glm_fit <- stats::glm(y_uni ~ 0 + X_uni, family = fam)
      uni_coefs[, ii] <- stats::coef(glm_fit)
    }
    Beta <- apply(uni_coefs, 1, mean, na.rm = T) # Could also weight by SE
  }

  stopifnot(is.list(Beta) | is.numeric(Beta))
  if(is.list(Beta)){
    stopifnot(length(Beta) == r,
              all(sapply(Beta, is.numeric)),
              sum(sapply(Beta, length)) == p)
    Beta <- unlist(Beta)
  }

  if(missing(W)) W <- matrix(X %*% Beta, nrow = n, ncol = r, byrow = TRUE)

  if(missing(Sigma)) Sigma <- diag(1e-3, r) # For small var, model approx. GLM

  #############################################################################
  # Univariate fitting (recursive calling and only if uni_fit = TRUE)
  #############################################################################
  if(uni_fit & r > 1){
    Sigma <- diag(diag(Sigma), r)
    for(ii in 1:r){
      Xi <- X_list[[ii]]
      Yi <- Y[, ii, drop = F]
      beta_idx <- seq(cumsum(n_pred)[ii] - n_pred[ii] + 1,
                      length.out = n_pred[ii])
      fit_uni <- lvmmr(Y = Yi,
                           X = Xi,
                           type = type[ii],
                           M = M[ii, ii, drop = F],
                           relative = relative,
                           quiet = quiet,
                           maxit = maxit,
                           tol = tol,
                           psi = psi[ii],
                           pgd = pgd,
                           Beta = Beta[beta_idx],
                           Sigma = Sigma[ii, ii, drop = F],
                           W = W[, ii, drop = F])
      Beta[beta_idx] <- fit_uni$Beta
      Sigma[ii, ii] <- fit_uni$Sigma
      W[, ii] <- fit_uni$W
    }
    D1 <- t(get_cumulant_diffs(W_T = t(W), type = type, order = 1))
    D2 <- t(get_cumulant_diffs(W_T = t(W), type = type, order = 2))
    end_obj <- -working_ll_rcpp(Y_T = t(Y), X_T = t(X), beta = Beta,
                                Sigma = Sigma, W_T = t(W), psi = psi,
                                D1_T = t(D1), D2_T = t(D2))
    return(list(Beta = unname(Beta), Sigma = Sigma, W = W, iter = NULL,
                change = NULL, obj = end_obj,
                X = X, Y = Y, M = M, type = type, psi = psi))
  }

  #############################################################################
  # Main algorithm
  #############################################################################
  out_iter <- 0
  iterate_outer <- out_iter < maxit[1] # Iterate updating (Beta, Sigma) and W
  while(iterate_outer){
    # For inner loop with joint update of Beta and Sigma
    in_iter <- 0
    iterate_inner <- in_iter < maxit[2]

    # Avoid using Beta and Sigma storage in inner loop.
    # Starting value of Sigma is set to be PD and satisfy constraints
    # NB: This is only for the starting value, and is done for stability
    # since the W update may have led to indefinite (for some i)
    #     Ci = diag(D2[ii,]) Sigma diag(D2[ii,]) + diag(D2[ii,]) diag(psi)
    # at the previous iterate of Sigma. Not an issue if pgd = TRUE.
    new_Beta <- Beta
    new_Sigma <- project_rcpp(X = Sigma,
                              restr_idx = which(c(!is.na(M))),
                              restr = as.numeric(M[!is.na(M)]),
                              eps = sqrt(.Machine$double.eps),
                              tol = sqrt(.Machine$double.eps),
                              maxit = 1e4)
    # Pre-compute
    D1 <- t(get_cumulant_diffs(W_T = t(W), type = type, order = 1))
    D2 <- t(get_cumulant_diffs(W_T = t(W), type = type, order = 2))

    while(iterate_inner){
      # Track progress of Beta and Sigma update
      start_obj <- -working_ll_rcpp(Y_T = t(Y), X_T = t(X), beta = new_Beta,
                                    Sigma = new_Sigma, W_T = t(W), psi = psi,
                                    D1_T = t(D1), D2_T = t(D2))
      # Update Sigma
      R <- matrix(X %*% new_Beta, nrow = n, ncol = r, byrow = T)
      R <- D1 + D2 * (R - W)
      R <- Y - R
      if(pgd){
        new_Sigma <- update_Sigma_pgd(R = R, D2 = D2, psi = psi,
                                       Sigma.init = new_Sigma,
                                       M = M, epsilon = eps,
                                       tol.dykstra = tol[3],
                                       tol.ipiano = tol[3],
                                       max.iter.dykstra = maxit[3],
                                       max.iter.ipiano = maxit[3],
                                       quiet = quiet[3])
      } else{
        new_Sigma <- update_Sigma_trust(Sigma_start = new_Sigma, R = R,
                                        D2 = D2, psi = psi, M = M,
                                        use_idx = 1:n, maxit = maxit[3],
                                        tol = tol[3])
      }

      if(!quiet[2]){
        mid_obj <- -working_ll_rcpp(Y_T = t(Y), X_T = t(X), beta = new_Beta,
                                   Sigma = new_Sigma, W_T = t(W), psi = psi,
                                   D1_T = t(D1), D2_T = t(D2))
        cat("Change from Sigma update: ", mid_obj - start_obj, "\n")
        if(mid_obj - start_obj > sqrt(tol[2])){
          warning("Sigma update increased objective function more than
          sqrt(tol[2]). \n")
        }
      }

      # Update Beta
      new_Beta <- update_beta(Y = Y, X = X, W = W, Sigma = new_Sigma,
                              psi = psi, type = type)

      # Track progress of Beta and Sigma update
      end_obj <- -working_ll_rcpp(Y_T = t(Y), X_T = t(X), beta = new_Beta,
                                  Sigma = new_Sigma, W_T = t(W), psi = psi,
                                  D1_T = t(D1), D2_T = t(D2))
      if(!quiet[2]){
        cat("Change from Beta update: ", end_obj - mid_obj, "\n")
        if(end_obj - mid_obj > sqrt(tol[2])){
          warning("Beta update increased objective function more than
          sqrt(tol[2]). \n")
        }
      }

      # Check whether to terminate inner loop
      in_iter <- in_iter + 1
      change <- abs(end_obj - start_obj)
      if(relative) change <- change / abs(start_obj)
      iterate_inner <- ((change > tol[2]) & (in_iter < maxit[2]))
    } # End inner loop

    # Update W
    pen <- ifelse(missing(w_pen),
                  max(abs(eigen(new_Sigma, only.values = TRUE)$val)),
                  w_pen)
    W <- update_W(Y = Y, X = X, W = W, Beta = new_Beta, Sigma = new_Sigma,
                    psi = psi, type = type, tol = tol[4],
                    maxit = maxit[4], quiet = quiet[4],
                    pen = pen)

    # Check whether to terminate outer loop
    if(relative){
      change <- max(abs(c((Sigma - new_Sigma) / max(abs(Sigma)),
                          (Beta - new_Beta) / max(abs(Beta)))))
    } else{
      change <- max(abs(c(Sigma - new_Sigma, Beta - new_Beta)))
    }
    out_iter <- out_iter + 1
    iterate_outer <-  ((change > tol[1]) & (out_iter < maxit[1]))

    if(!quiet[1]){
      cat("Change in parameters: ", change, "\n")
    }

    # Prepare next iteration
    Beta <- new_Beta
    Sigma <- new_Sigma
  } # End outer loop
  return(list(Beta = unname(Beta), Sigma = Sigma, W = W, iter = out_iter,
              change = change, obj = end_obj,
              X = X, Y = Y, M = M, type = type, psi = psi))
}

