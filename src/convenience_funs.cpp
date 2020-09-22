#include <RcppArmadillo.h>
// type = 1 means Normal, type = 2 Bernoulli, and type = 3 Poisson

// [[Rcpp::depends(RcppArmadillo)]]

arma::vec log1pexp(arma::vec x)
{
  for(int ii = 0; ii < x.n_elem; ++ii){
    if(x(ii) < 18){
      x(ii) = std::log1p(std::exp(x(ii)));
    } else{
      x(ii) = x(ii) + std::exp(-x(ii));
    }
  }
  return x;
}

double cumulant(double w, int type)
{
  if(type == 1){
    w = 0.5 * std::pow(w, 2);
  } else if(type == 2){
    w = log1pexp(w);
  } else if (type == 3){
    w = std::exp(w);
  } else{
    Rcpp::stop("type must be 1 (normal), 2 (Bernoulli), or 3 (Poisson).");
  }
  return w;
}

double cumulant_d(double w, int type)
{
  if(type == 1){
    // Identity map for normal
  } else if(type == 2){
    w = 1.0 / (1 + std::exp(-w));
  } else if (type == 3){
    w = std::exp(w);
  } else{
    Rcpp::stop("type must be 1 (normal), 2 (Bernoulli), or 3 (Poisson).");
  }
  return w;
}

double cumulant_dd(double w, int type)
{
  if(type == 1){
    w = 1.0;
  } else if(type == 2){
    w = 1.0 / (1 + std::exp(-w));
    w *= (1.0 - w);
  } else if (type == 3){
    w = std::exp(w);
  } else{
    Rcpp::stop("type must be 1 (normal), 2 (Bernoulli), or 3 (Poisson).");
  }
  return w;
}

double cumulant_d3(double w, int type)
{
  if(type == 1){
    w = 0.0;
  } else if(type == 2){
    w = 1.0 / (1.0 + std::exp(-w));
    w = w * ( 1.0 - w) * (1 - 2.0 * w);
  } else if (type == 3){
    w = std::exp(w);
  } else{
    Rcpp::stop("type must be 1 (normal), 2 (Bernoulli), or 3 (Poisson).");
  }
  return w;
}

double cumulant_d4(double w, int type)
{
  if(type == 1){
    w = 0.0;
  } else if(type == 2){
    w = 1.0 / (1.0 + std::exp(-w));
    w = w * ( 1.0 - w) * (6.0 * std::pow(w, 2) - 6.0 * w + 1.0);
  } else if (type == 3){
    w = std::exp(w);
  } else{
    Rcpp::stop("type must be 1 (normal), 2 (Bernoulli), or 3 (Poisson).");
  }
  return w;
}

// Supply transposes of W (better memory accessing pattern)

// [[Rcpp::export]]
arma::mat get_cumulant_diffs(arma::mat W_T, arma::ivec type, int order)
{
  if(order > 4){
    Rcpp::warning("Only order <= 4 supported, returning fourth derivatives.");
  }

  for(int ii = 0; ii < W_T.n_cols; ++ii)
  {
    for(int jj = 0; jj < W_T.n_rows; ++jj){
      if(order == 0){
        W_T(jj, ii) = cumulant(W_T(jj, ii), type(jj));
      } else if(order == 1){
        W_T(jj, ii) = cumulant_d(W_T(jj, ii), type(jj));
      } else if(order == 2){
        W_T(jj, ii) = cumulant_dd(W_T(jj, ii), type(jj));
      } else if(order == 3){
        W_T(jj, ii) = cumulant_d3(W_T(jj, ii), type(jj));
      } else{
        W_T(jj, ii) = cumulant_d4(W_T(jj, ii), type(jj));
      }
    }
  }
  return W_T;
}
