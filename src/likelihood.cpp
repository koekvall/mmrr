#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double working_ll_rcpp(const arma::mat Y_T,
                       const arma::mat X_T,
                       const arma::vec beta,
                       const arma::mat Sigma,
                       const arma::mat W_T,
                       const arma::vec psi,
                       const arma::mat D1_T,
                       const arma::mat D2_T)
{
  uint r = Y_T.n_rows;
  uint n = Y_T.n_cols;
  double val = 0.0;

  arma::vec Xb = X_T.t() * beta;

  for (size_t ii = 0; ii < n; ii++){
    // ith working covariance matrix
    arma::mat C = Sigma;
    C.each_row() %= D2_T.col(ii).t();
    C.each_col() %= D2_T.col(ii);
    C.diag() += psi % D2_T.col(ii);
    C = arma::symmatu(C);

    arma::mat U(r, r);
    arma::vec d(r);
    arma::eig_sym(d, U, C);

    // Return infinite objective if C is numerically ~singular
    if(d(0) < 1e-12){
      val = std::numeric_limits<double>::infinity();
      break;
    }

    // ith working mean and residual
    arma::vec m = Xb.subvec(ii * r, (ii + 1) * r - 1);
    m = D1_T.col(ii) + D2_T.col(ii) % (m - W_T.col(ii));
    m = Y_T.col(ii) - m;
    // ith contribution do objective function
    val += 0.5 * arma::sum(arma::log(d));
    val += 0.5 * arma::sum(arma::square((1.0 / arma::sqrt(d)) %
           (U.t() * m)));
  }
  val = -val - 0.5 * n * r * std::log(2.0 * arma::datum::pi);
  return val;
}
