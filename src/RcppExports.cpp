// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// get_cumulant_diffs
arma::mat get_cumulant_diffs(arma::mat W_T, arma::ivec type, int order);
RcppExport SEXP _lvmmr_get_cumulant_diffs(SEXP W_TSEXP, SEXP typeSEXP, SEXP orderSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type W_T(W_TSEXP);
    Rcpp::traits::input_parameter< arma::ivec >::type type(typeSEXP);
    Rcpp::traits::input_parameter< int >::type order(orderSEXP);
    rcpp_result_gen = Rcpp::wrap(get_cumulant_diffs(W_T, type, order));
    return rcpp_result_gen;
END_RCPP
}
// working_ll_rcpp
double working_ll_rcpp(const arma::mat Y_T, const arma::mat X_T, const arma::vec beta, const arma::mat Sigma, const arma::mat W_T, const arma::vec psi, const arma::mat D1_T, const arma::mat D2_T);
RcppExport SEXP _lvmmr_working_ll_rcpp(SEXP Y_TSEXP, SEXP X_TSEXP, SEXP betaSEXP, SEXP SigmaSEXP, SEXP W_TSEXP, SEXP psiSEXP, SEXP D1_TSEXP, SEXP D2_TSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type Y_T(Y_TSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type X_T(X_TSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type Sigma(SigmaSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type W_T(W_TSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type psi(psiSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type D1_T(D1_TSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type D2_T(D2_TSEXP);
    rcpp_result_gen = Rcpp::wrap(working_ll_rcpp(Y_T, X_T, beta, Sigma, W_T, psi, D1_T, D2_T));
    return rcpp_result_gen;
END_RCPP
}
// project_rcpp
arma::mat project_rcpp(arma::mat X, const arma::uvec restr_idx, const arma::vec restr, const double eps, const double tol, uint maxit);
RcppExport SEXP _lvmmr_project_rcpp(SEXP XSEXP, SEXP restr_idxSEXP, SEXP restrSEXP, SEXP epsSEXP, SEXP tolSEXP, SEXP maxitSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::uvec >::type restr_idx(restr_idxSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type restr(restrSEXP);
    Rcpp::traits::input_parameter< const double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< const double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< uint >::type maxit(maxitSEXP);
    rcpp_result_gen = Rcpp::wrap(project_rcpp(X, restr_idx, restr, eps, tol, maxit));
    return rcpp_result_gen;
END_RCPP
}
// obj_sigma_rcpp
Rcpp::List obj_sigma_rcpp(arma::mat Sigma, const arma::mat R_T, const arma::mat D2_T, const arma::vec psi, const arma::ivec use_idx, const uint order);
RcppExport SEXP _lvmmr_obj_sigma_rcpp(SEXP SigmaSEXP, SEXP R_TSEXP, SEXP D2_TSEXP, SEXP psiSEXP, SEXP use_idxSEXP, SEXP orderSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Sigma(SigmaSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type R_T(R_TSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type D2_T(D2_TSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type psi(psiSEXP);
    Rcpp::traits::input_parameter< const arma::ivec >::type use_idx(use_idxSEXP);
    Rcpp::traits::input_parameter< const uint >::type order(orderSEXP);
    rcpp_result_gen = Rcpp::wrap(obj_sigma_rcpp(Sigma, R_T, D2_T, psi, use_idx, order));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_lvmmr_get_cumulant_diffs", (DL_FUNC) &_lvmmr_get_cumulant_diffs, 3},
    {"_lvmmr_working_ll_rcpp", (DL_FUNC) &_lvmmr_working_ll_rcpp, 8},
    {"_lvmmr_project_rcpp", (DL_FUNC) &_lvmmr_project_rcpp, 6},
    {"_lvmmr_obj_sigma_rcpp", (DL_FUNC) &_lvmmr_obj_sigma_rcpp, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_lvmmr(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
