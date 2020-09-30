working_ll <- function(Y, X, Beta, Sigma, W, psi, D1, D2)
{
  n <- nrow(Y)
  r <- ncol(Y)

  Xb <- matrix(X %*% Beta, nrow = n, ncol = r, byrow = T)

  ll <- 0
  for(ii in 1:n){
    m <- D1[ii, ] + D2[ii, ] * (Xb[ii, ] - W[ii, ])
    C <- Sigma * tcrossprod(D2[ii, ]) # Works because drop = TRUE in D2 index
    diag(C) <- diag(C) + psi * D2[ii, ]
    ll <- ll + mvtnorm::dmvnorm(x = Y[ii, ],
                                mean = m,
                                sigma = C,
                                log = TRUE,
                                checkSymmetry = FALSE)
  }
  return(ll)
}

obj_Sigma <- function(Sigma, R, D2, psi, use_idx)
{
  n <- nrow(D2)
  r <- ncol(D2)
  out <- 0
  if(missing(use_idx)) use_idx <- 1:n
  for(ii in use_idx){
    C <- Sigma * tcrossprod(D2[ii, ])
    diag(C) <- diag(C) + psi * D2[ii, ]
    e_C <- eigen(C)
    if(e_C$values[r] < 1e-12) return(Inf)
    out <- out + sum(log(e_C$values))
    out <- out + sum((sqrt(1 / e_C$values) *
                        crossprod(e_C$vectors, R[ii, ]))^2)
  }
  return(out)
}

gradient_Sigma <- function(Sigma, R, D2, psi, use_idx)
{
  n <- nrow(D2)
  r <- ncol(D2)
  G <- matrix(0, r, r)
  if(missing(use_idx)) use_idx <- 1:n
  for(ii in use_idx){
    C <- Sigma * tcrossprod(D2[ii, ])
    diag(C) <- diag(C) + psi * D2[ii, ]
    e_C <- eigen(C)
    I <-  e_C$vectors %*% (1 / e_C$values * t(e_C$vectors))
    C <- C - tcrossprod(R[ii, ])
    C <- I %*% C %*% I
    C <- D2[ii, ] * C
    C <- t(D2[ii, ] * t(C))
    G <- G + C
  }
  return(G)
}

fisher_Sigma <- function(Sigma, D2, psi, use_idx)
{
  n <- nrow(D2)
  r <- ncol(D2)

  if(missing(use_idx)) use_idx <- 1:n

  H <- matrix(0, r^2, r^2)
  for(ii in use_idx){
    C <- Sigma
    diag(C) <- diag(C) + psi / D2[ii, ]
    e_C <- eigen(C)
    C <- e_C$vectors %*% (1 / e_C$values * t(e_C$vectors))
    H <- H + kronecker(C, C)
  }
  return(H)
}

hessian_Sigma <- function(Sigma, R, D2, psi, use_idx){
  n <- nrow(D2)
  r <- ncol(D2)
  H <- matrix(0, r^2, r^2)
  if(missing(use_idx)) use_idx <- 1:n
  for(ii in use_idx){
    C <- Sigma
    diag(C) <- diag(C) + psi / D2[ii, ]
    e_C <- eigen(C)
    C <- e_C$vectors %*% (1 / e_C$values * t(e_C$vectors))
    v <- matrix(C %*% (R[ii, ] / D2[ii, ]), ncol = 1)
    H <- H + kronecker(v, kronecker(t(v), C))
    H <- H + kronecker(kronecker(C, v), t(v))
    H <- H - kronecker(C, C)
  }
  return(H)
}

# vech <- function(X)
# {
#   return(c(X[lower.tri(X, diag = T)]))
# }
#
# vech_inv <- function(x)
# {
#   d <- 0.5 * (-1  + sqrt(1 + 8 * length(x)))
#   X <- matrix(0, d, d)
#   X[lower.tri(X, diag = T)] <- x
#   X[upper.tri(X)] <- t(X)[upper.tri(X)]
#   return(X)
# }
#
# vec_to_vech <- function(x){
#   d <- sqrt(length(x))
#   return(x[c(lower.tri(matrix(NA, d, d), diag = T))])
# }
#
# vech_to_vec <- function(x){
#   return(c(vech_inv(x)))
# }
