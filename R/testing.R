#' Approximate likelihood ratio testing for Latent Variables Multivariate
#'    Mixed-type Response Regression
#'
#' @param fit_null The output from lvmmr() fit under the null hypothesis
#' @param fit_full The output from lvmmr() fit under the alternative hypothesis
#' @param df The degrees of freedom; if not supplied, set to the sum of the
#'    difference in number of predictors and number of restrictions on Sigma.
#' @param ... Other arguments passed to lvmmr() when re-fitting the full model
#'   holding W fixed at fit_null$W.
#' @return A list with the test statistic, p-value, and degrees of freedom.
#' @export
lrt_approx <- function(fit_null, fit_full, df, ...)
{
  # Unconstrained fit with null W
  fit_new <- lvmmr(Y = fit_full$Y, X = fit_full$X, type = fit_full$type,
                   psi = fit_full$psi, M = fit_full$M, maxit = c(10, 500, 500, 0),
                   Beta = fit_null$Beta, Sigma = fit_null$Sigma,
                   W = fit_null$Sigma, ...)
  if(missing(df)){
    df <- ncol(fit_full$X) - ncol(fit_null$X)
    df <- df + sum(is.na(fit_full$M[upper.tri(fit_full$M, diag = TRUE)]))
    df <- df - sum(is.na(fit_null$M[upper.tri(fit_null$M, diag = TRUE)]))
  }
  test_stat <- 2 * (fit_null$obj - fit_new$obj)
  return(list(stat = test_stat,
              p = stats::pchisq(test_stat, df = df, lower.tail = FALSE),
              df = df))
}
