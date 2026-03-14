#include <RcppArmadillo.h>
#include <vector>
#include <cmath>
#include <algorithm>

#ifdef _OPENMP
#include <omp.h>
#endif


// [[Rcpp::depends(RcppArmadillo)]]


using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List gw_reg(const arma::mat& x, const arma::vec& y, const arma::vec& w, const arma::mat& d, bool hatmatrix, int focus) {
  arma::vec wd = w % d.diag();
  arma::mat xtw = trans(x) * diagmat(wd);
  arma::mat xtwx = xtw * x;
  arma::mat xtwy = trans(x) * diagmat(wd) * y;
  arma::mat xtwx_inv = inv(xtwx);
  arma::vec beta = xtwx_inv * xtwy;
  
  if (hatmatrix) {
    arma::mat ci = xtwx_inv * xtw;
    arma::rowvec s_ri = x.row(focus - 1) * ci; // hatmatrix
    return List::create(
      Named("beta") = beta,
      Named("S_ri") = s_ri,
      Named("Ci") = ci
    );
  }
  else {
    return List::create(
      Named("beta") = beta
    );
  }
}



// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List update_all_rows_cpp(const arma::mat& Si, const arma::mat& S_arrayi, const arma::mat& Shat, int num_threads = 1) {
  int n = Si.n_rows;
  arma::mat temp1(n, n, arma::fill::zeros);
  arma::mat temp2(n, n, arma::fill::zeros);
  
  omp_set_num_threads(num_threads);
  
#pragma omp parallel for
  for (int row = 0; row < n; ++row) {
    for (int col = 0; col < n; ++col) {
      temp1(row, col) = arma::dot(Si.row(row), S_arrayi.col(col));
      temp2(row, col) = arma::dot(Si.row(row), Shat.col(col)); 
    }
  }
  
  return Rcpp::List::create(Rcpp::Named("temp1") = temp1,
                            Rcpp::Named("temp2") = temp2);
}




// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat calculateBetaSE_Cpp(const arma::mat& x1, const arma::vec& dd, const arma::cube& S_arrays) {
  int n = x1.n_rows;
  int var_n = x1.n_cols;
  arma::mat Beta_SE(n, var_n, arma::fill::zeros);
  
  for (int i = 0; i < var_n; i++) {
    arma::vec x1_col = x1.col(i);
    arma::mat Ci = arma::diagmat(1 / x1_col) * S_arrays.slice(i);
    arma::mat Ci_dd_CiT = Ci * arma::diagmat(1 / dd) * Ci.t();
    Beta_SE.col(i) = arma::sqrt(Ci_dd_CiT.diag());
  }
  
  return Beta_SE;
}
