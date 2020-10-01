#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List obj_sigma_rcpp(arma::mat Sigma,
                          const arma::mat R_T,
                          const arma::mat D2_T,
                          const arma::vec psi,
                          const arma::ivec use_idx,
                          const uint order)
{
  uint r = R_T.n_rows;
  uint n_eff = use_idx.n_elem;
  double val = 0.0;
  arma::mat G(r, r, arma::fill::zeros);
  arma::mat H(r * r, r * r, arma::fill::zeros);
  arma::mat I(r, r); // Storage used in calculations
  for (size_t ii = 0; ii < n_eff; ii++){
    size_t idx = use_idx(ii) - 1; // C++ starts indexing at zero

    // Create ith covariance matrix
    arma::mat C = Sigma;
    C.each_row() %= D2_T.col(idx).t();
    C.each_col() %= D2_T.col(idx);
    C.diag() += psi % D2_T.col(idx);

    arma::mat U(r, r);
    arma::vec d(r);
    arma::eig_sym(d, U, C);

    // Return infinite objective if C is numerically ~singular
    if(d(0) < 1e-12){
      val = std::numeric_limits<double>::infinity();
      break;
    }

    // ith contribution do objective function
    val += arma::sum(arma::log(d));
    val += arma::sum(arma::square((1.0 / arma::sqrt(d)) %
           (U.t() * R_T.col(idx))));

    // ith contribution to gradient
    if(order >= 1){
      I = U * arma::diagmat(1.0 / d) * U.t();
      I.each_col() %= D2_T.col(idx);
      C -= R_T.col(idx) * R_T.col(idx).t();
      C = I * C * I.t();
      G += C;
    }

    // ith contribution to Hessian
    if(order >= 2){
      I.each_row() %= D2_T.col(idx).t();
      arma::mat v = I * (R_T.col(idx) / D2_T.col(idx));
      H += arma::kron(v, arma::kron(v.t(), I));
      H += arma::kron(arma::kron(I, v), v.t());
      H -= arma::kron(I, I);
    }
  }
  return Rcpp::List::create(Rcpp::Named("value") = val,
                          Rcpp::Named("gradient") = G,
                          Rcpp::Named("hessian") = H);
}
