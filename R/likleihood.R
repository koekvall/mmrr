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
