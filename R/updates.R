update_beta <- function(Y, X, W, Sigma, psi, type)
{
  n <- nrow(Y)
  p <- ncol(X)
  r <- ncol(Y)

  D1 <- t(get_cumulant_diffs(t(W), type, 1)) # t(); C++ code has obs by column
  D2 <- t(get_cumulant_diffs(t(W), type, 2))

  # Allocate
  H <- matrix(0, nrow = p, ncol = p + 1)
  for(ii in 1:n){
    # Working covariance
    C <- Sigma * tcrossprod(D2[ii, ]) # Works because drop = TRUE in D2 index
    diag(C) <- diag(C) + psi * D2[ii, ]

    # Working predictor
    start_idx <- (ii - 1) * r + 1
    Zi <- X[seq(start_idx, start_idx + r - 1), , drop = F]

    # Scale jth row by \nabla c_j(w_j), j = 1, ... r, by recycling.
    Zi <- Zi * D2[ii, ]

    # Working response
    yi <- Y[ii, ] - D1[ii, ] + D2[ii, ] * W[ii, ]

    # Sum up individual contributions
    H <- H + crossprod(Zi, qr.solve(C, cbind(Zi, yi)))
  }
  return(qr.solve(H[, 1:p], H[, p + 1]))
}

update_W <- function(Y, X, W, Beta, Sigma, psi, type, pen = 1e-4, tol = 1e-8,
                     maxit = 100, quiet = T)
{
  if(maxit <= 0) return(W) # Allows maxit[4] = 0 to skip W update
  n <- nrow(Y)
  r <- ncol(Y)

  # Set lower limit on eigenvalues for stability
  e_S <- eigen(Sigma)
  e_S$values <- pmax(e_S$values, sqrt(.Machine$double.eps))

  # Replace by "inverse"
  Sigma <- e_S$vectors %*% (e_S$values * t(e_S$vectors))

  # Precompute
  Xb <- matrix(X %*% Beta, nrow = n, ncol = r, byrow = T)

  # This update uses ui = W[ii, ] - Xb[ii, ]
  # The ith objective is:
  #     -Y[ii, ] * u + c(u + Xb[ii, ]) + 0.5 * t(u) %*% solve(Sigma) %*% u +
  #         0.5 * pen ||u||^2

  for(ii in 1:n){
    trust_obj <- function(u){
      u <- as.matrix(u, ncol = 1)
      d0 <- as.numeric(get_cumulant_diffs(u + Xb[ii, ], type, 0))
      d1 <- as.numeric(get_cumulant_diffs(u + Xb[ii, ], type, 1))
      d2 <- as.numeric(get_cumulant_diffs(u + Xb[ii, ], type, 2))

      # Objective
      val <-  sum(d0 - Y[ii, ] * u)
      val <- as.numeric(val + 0.5 * crossprod(u, Sigma %*% u))
      val <- val + 0.5 * pen * sum(u^2)

      # Gradient
      g <- d1 - Y[ii, ]

      g <- as.numeric(g + Sigma %*% u + pen * u)

      # Hessian (modified)
      H <- Sigma

      diag(H) <- diag(H) + d2 + pen

      return(list(value = val, gradient = g, hessian = H))
    }
    opt <- trust::trust(trust_obj, W[ii, ] - Xb[ii, ], rinit = 1, rmax = 100,
                        iterlim = maxit, fterm = tol, mterm = tol)

    W[ii, ] <- opt$argument + Xb[ii, ]
    if(!opt$converged){
      warning(paste0("w_", ii, " update did not converge \n"))
    }
  }
  return(W)
}

update_Sigma_proj <- function(R, D2, psi, Sigma.init, M, epsilon = 0,
                              tol.dykstra = 1e-12, tol.ipiano = 1e-10,
                              max.iter.dykstra = 1e3, max.iter.ipiano = 1e3,
                              quiet = TRUE)
{
  Sigmakm1 <- Sigma.init
  Sigma <- Sigma.init
  L0 <- 2
  delta <- 1e-4
  c2 <- 1e-6
  # obj.prev <- evalObj(H, A, B, Sigma)
  obj.prev <- obj_sigma_rcpp(Sigma= Sigma,
                             R_T = t(R),
                             D2_T = t(D2),
                             psi = psi,
                             use_idx = 1:nrow(R),
                             order = 0)$value
  obj.orig <- obj.prev

  for(kk in 1:max.iter.ipiano){

    #tempGrad <- getGrad(H, A, B, Sigma)
    tempGrad <- obj_sigma_rcpp(Sigma = Sigma,
                              R_T = t(R),
                              D2_T = t(D2),
                              psi = psi,
                              use_idx = 1:nrow(R),
                              order = 1)$gradient
    Ln <- L0
    linesearch <- TRUE

    while(linesearch){

      # -- interial step projected grad ---
      # b <- (delta + Ln/2)/(c2 + Ln/2)
      # Bn <- (b-1)/(b-.5)
      # alpha <- 2*(1 - Bn)/(2*c2 + Ln)
      Bn <- .95
      alpha <- 1.9*(1 - Bn)/Ln
      temp <- Sigma - alpha*tempGrad + Bn*(Sigma - Sigmakm1)
      #Sigma.temp <- CorrelationProjection(temp, M = M, epsilon= epsilon)
      Sigma.temp <- project_rcpp(X = temp,
                                 restr_idx = which(c(!is.na(M))),
                                 restr = as.numeric(M[!is.na(M)]),
                                 eps = epsilon,
                                 tol = tol.dykstra,
                                 maxit = max.iter.dykstra)
      # obj.temp <- evalObj(H, A, B, Sigma.temp)
      obj.temp <- obj_sigma_rcpp(Sigma = Sigma.temp,
                                 R_T = t(R),
                                 D2_T = t(D2),
                                 psi = psi,
                                 use_idx = 1:nrow(R),
                                 order = 0)$value
      if(obj.temp < obj.prev + sum(tempGrad*t(Sigma.temp - Sigma)) + (Ln/2)*sum((Sigma.temp - Sigma)^2)){
        Sigmakm1 <- Sigma
        Sigma <- Sigma.temp
        linesearch <- FALSE
      } else {
        Ln <- Ln*5
      }
    }

    if(abs(obj.temp - obj.prev) < tol.ipiano*abs(obj.orig)){
      break
    }
    obj.prev <- obj.temp
    if(!quiet){
      cat(obj.prev, "\n")
    }

  }

  return(Sigma)
}

update_Sigma_trust <- function(Sigma_start, R, D2, psi, M, use_idx)
{
  n <- nrow(R)
  r <- ncol(R)

  m <- matrixcalc::vech(M)
  opt_idx <- lower.tri(M, diag = T) & is.na(M)

  theta0 <- as.numeric(Sigma_start[opt_idx])

  D <- matrixcalc::D.matrix(r) # This is very inefficient

  obj <- function(theta){
    Sigma <- M
    Sigma[opt_idx] <- theta
    Sigma[upper.tri(Sigma)] <- t(Sigma)[upper.tri(Sigma)]

    all_outs <- obj_sigma_rcpp(Sigma = Sigma, R_T = t(R), D2_T = t(D2),
                               psi = psi, use_idx = 1:n, order = 2)
    g <- as.numeric(crossprod(D, as.numeric(all_outs$gradient)))[is.na(m)]

    H <- crossprod(D, all_outs$hessian %*% D)[is.na(m), is.na(m)]
    return(list(value = all_outs$value, gradient = g, hessian = H))
  }

    opt <- trust::trust(objfun = obj, parinit = theta0, rinit = 1, rmax = 100)

    if(!opt$converged){
      warning("Sigma update did not converge \n")
    }
    Sigma <- Sigma_start
    Sigma[opt_idx] <- opt$argument
    Sigma[!is.na(M)] <- M[!is.na(M)]
    Sigma[upper.tri(Sigma)] <- t(Sigma)[upper.tri(Sigma)]

    return(Sigma)
}
