#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat project_rcpp(arma::mat X, const arma::uvec restr_idx,
  const arma::vec restr, const double eps, const double tol, uint maxit)
{
  const uint d = X.n_cols;
  const double Inf = std::numeric_limits<double>::infinity();

  arma::mat Xk(d, d);
  if(restr_idx.n_elem == 0){
    arma::vec e(d);
    X = arma::symmatu(X);
    arma::eig_sym(e, X, X);
    Xk = X * arma::diagmat(arma::clamp(e, eps, Inf)) * X.t();
  } else{
    arma::mat Xkm1 = X;
    arma::mat Pkm1(d, d, arma::fill::zeros);
    arma::mat Qkm1(d, d, arma::fill::zeros);
    double obj_old = 0.0;
    for (size_t kk = 0; kk < maxit; kk++) {
      if(kk >= maxit - 1){
        Rcpp::warning("Projection algorithm reached max. iter. \n");
      }
      arma::mat Ykm1 = Xkm1 + Pkm1;
      Ykm1.elem(restr_idx - 1) = restr; // C++ index start at zero
      arma::mat Pk = Xkm1 + Pkm1 - Ykm1;
      arma::mat U = Ykm1 + Qkm1;
      arma::vec e(d);
      U = arma::symmatu(U);
      arma::eig_sym(e, U, U);
      Xk = U * arma::diagmat(arma::clamp(e, eps, Inf)) * U.t();
      arma::mat Qk = Ykm1 + Qkm1 - Xk;
      double obj_new = 0.5 * arma::accu(arma::square((Xk - X)));
      if(std::abs(obj_new - obj_old) < tol){
        break;
      }
      // prepare next iteration
      Xkm1 = Xk;
      Pkm1 = Pk;
      Qkm1 = Qk;
      obj_old = obj_new;
    }
  }
  return Xk;
}
