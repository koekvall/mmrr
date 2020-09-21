working_ll <- function(Y, X, Beta, Sigma, W, psi, D1, D2)
{
  n <- nrow(Y)
  p <- ncol(X)
  r <- nrol(Y)
  
  Xb <- matrix(X %*% Beta, nrow = n, ncol = r, byrow = T)
  
  ll <- 0
  for(ii in 1:n){
    m <- D1[ii, ] + D2[ii, ] * (Xb[ii, ] - W[ii, ])
    C <- sweep(Sigma, 2, D2[ii, ], FUN = "*")
    C <- sweep(C, 1, D2[ii, ], FUN = "*")
    diag(C) <- diag(C) + psi * D2[ii, ]
    ll <- ll + mvtnorm::dmvnorm(x = m, 
                                sigma = C,
                                log = TRUE,
                                checkSymmetry = FALSE)
  }
  return(ll)
}