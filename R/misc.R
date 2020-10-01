# -----------------------------------------------
# Dykstra algorithm for alternating projections
# ----------------------------------------------
project_pd_M <- function(input, M, epsilon, tol, max.iter){

  eval.objective <- function(X, input){
    return(0.5*sum((X - input)^2))
  }
  Xkm1 <- input
  Pkm1 <- matrix(0, nrow=nrow(input), ncol=ncol(input))
  Qkm1 <- Pkm1
  obj.old <- 0

  for(kk in 1:max.iter){
    if(kk >= max.iter){
      warning("Projection algorithm reached max.iter \n")
    }
    # -----------------------
    # Project onto M
    # -----------------------
    Ykm1 <- Xkm1 + Pkm1
    Ykm1[!is.na(M)] <- M[!is.na(M)]
    Pk <- Xkm1 + Pkm1 - Ykm1
    eo <- eigen(Ykm1 + Qkm1)
    Xk <- eo$vec %*% (pmax(eo$val, epsilon) * t(eo$vec))
    Qk <- Ykm1 + Qkm1 - Xk
    obj.new <- eval.objective(Xk, input)
    if(abs(obj.new - obj.old) < tol){
      break
    }
    Xkm1 <- Xk
    Pkm1 <- Pk
    Qkm1 <- Qk

    obj.old <- obj.new
  }

  return(Xk)
}
