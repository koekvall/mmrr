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


# -------------------------------------------------
# ipiano alg for covariance update with M
# -------------------------------------------------
# H: n x r matrix with ith row being the mean for ith subject
# A: n x r matrix with ith row being Psi-nabla^2 c(w)i
# B: n x r matrix with ith row being nabla^2 c(w)_i
# Sigma.init: initial value for Sigma (e.g., previous iteration)
# M: r x r matrix with constraints -- NA means unconstrained, non-NA gives constrained value
# epsilon: lower bound on Sigma eigenvalues (zero by default)
# tol.dykstra: tolerance for projection algorithm
# tol.ipiano: tolerance for projected gradient descent algorithm
# -
update_Sigma_PQL <- function(H, A, B, Sigma.init, M, epsilon = 0, tol.dykstra = 1e-12, tol.ipiano = 1e-10,
                             max.iter.dykstra = 1e3, max.iter.ipiano = 1e3, quiet=TRUE){

  CorrelationProjection <- function(input, M = M, epsilon = epsilon, tol = tol.dykstra, max.iter = max.iter.dykstra){

    eval.objective <- function(X, input){
      return(0.5*sum((X - input)^2))
    }
    # -----------------------------------------------
    # Dykstra algorithm for alternating projections
    # ----------------------------------------------
    Xkm1 <- input
    Pkm1 <- matrix(0, nrow=nrow(input), ncol=ncol(input))
    Qkm1 <- Pkm1
    obj.old <- 0

    for(kk in 1:max.iter){

      # -----------------------
      # Project onto M
      # -----------------------
      Ykm1 <- Xkm1 + Pkm1
      Ykm1[!is.na(M)] <- M[!is.na(M)]
      Pk <- Xkm1 + Pkm1 - Ykm1
      eo <- eigen(Ykm1 + Qkm1)
      Xk <- tcrossprod(eo$vec*(rep(1, nrow(input))%*%t(pmax(eo$val, epsilon))),
                       eo$vec)
      Qk <- Ykm1 + Qkm1 - Xk
      obj.new <- eval.objective(Xk, input)
      #cat(obj.new, "\n")
      if(abs(obj.new - obj.old) < tol){
        break
      }
      Xkm1 <- Xk
      Pkm1 <- Pk
      Qkm1 <- Qk

      obj.old <- obj.new
    }
    #Xk[!is.na(M)] <- Xk[!is.na(M)]
    return(Xk)
  }

  # -----------------------------------------
  # Projected gradient descent
  # -----------------------------------------
  getGrad <- function(H, A, B, Sigma){
    n <- dim(H)[1]
    out <- matrix(0, nrow=nrow(Sigma), ncol=ncol(Sigma))
    r <- nrow(Sigma)
    for(k in 1:n){
      Omega <- chol2inv(chol(diag(A[k,]) + tcrossprod(B[k,])*Sigma))
      out <- out + tcrossprod(B[k,])*crossprod(Omega, diag(1, r) - tcrossprod(H[k,])%*%Omega)
    }
    return(out)
  }

  evalObj <- function(H, A, B, Sigma){
    out <- 0
    n <- dim(H)[1]
    for(k in 1:n){
      Omega <- chol2inv(chol(diag(A[k,]) + tcrossprod(B[k,])*Sigma))
      out <- out + tcrossprod(crossprod(H[k,],Omega), H[k,]) - determinant(Omega, logarithm=TRUE)$modulus[1]
    }
    return(out)
  }

  Sigmakm1 <- Sigma.init
  Sigma <- Sigma.init
  L0 <- 2
  delta <- 1e-4
  c2 <- 1e-6
  obj.prev <- evalObj(H, A, B, Sigma)
  obj.orig <- obj.prev

  for(kk in 1:max.iter.ipiano){

    tempGrad <- getGrad(H, A, B, Sigma)
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
      Sigma.temp <- CorrelationProjection(temp, M = M, epsilon= epsilon)
      obj.temp <- evalObj(H, A, B, Sigma.temp)
      if(obj.temp < obj.prev + sum(tempGrad*t(Sigma.temp - Sigma)) + (Ln/2)*sum((Sigma.temp - Sigma)^2)){
        Sigmakm1 <- Sigma
        Sigma <- Sigma.temp
        linesearch <- FALSE
      } else {
        Ln <- Ln*5
        # cat("decreasing Ln", "\n")
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
