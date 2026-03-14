#include <Rcpp.h>
#include <omp.h>
#include <unordered_set>
#include <unordered_map>
using namespace Rcpp;

// [[Rcpp::export]]
List cox_derivatives_cpp(NumericVector eta, NumericVector time, NumericVector status, int num_threads = 1) {
  omp_set_num_threads(num_threads);
  
  int n = eta.size();
  NumericVector first_derivative(n);
  NumericVector second_derivative(n);
  
  // R_r 계산: 각 시간에서 risk set에 있는 i들의 index
  std::vector<std::vector<int>> R_list(n);
#pragma omp parallel for schedule(dynamic)
  for (int r = 0; r < n; ++r) {
    for (int i = 0; i < n; ++i) {
      if (time[i] >= time[r]) {
#pragma omp critical
        R_list[r].push_back(i);
      }
    }
  }
  
  // C_i 계산: 각 i가 포함된 모든 risk set의 index
  std::vector<std::vector<int>> C_list(n);
#pragma omp parallel for schedule(dynamic)
  for (int i = 0; i < n; ++i) {
    for (int r = 0; r < n; ++r) {
      if (std::find(R_list[r].begin(), R_list[r].end(), i) != R_list[r].end()) {
#pragma omp critical
        C_list[i].push_back(r);
      }
    }
  }
  
#pragma omp parallel for schedule(dynamic)
  for (int i = 0; i < n; ++i) {
    double sum_term_first = 0;
    double sum_term_second1 = 0;
    double sum_term_second2 = 0;
    
    for (int r : C_list[i]) {
      const std::vector<int>& R_r = R_list[r];
      double sum_exp_eta_R_r = 0;
      double d_r = 0;
      
      // risk set의 exp(eta) 합 계산
#pragma omp parallel for reduction(+:sum_exp_eta_R_r)
      for (int j = 0; j < R_r.size(); ++j) {
        sum_exp_eta_R_r += exp(eta[R_r[j]]);
      }
      
      // d_r 계산 (time[r]에서 event가 발생한 individual 수)
#pragma omp parallel for reduction(+:d_r)
      for (int j = 0; j < n; ++j) {
        if (time[j] == time[r] && status[j] == 1) {
          d_r += 1;
        }
      }
      
      double inv_sum_exp_eta_R_r = 1.0 / sum_exp_eta_R_r;
      double term = d_r * inv_sum_exp_eta_R_r;
      
      // first_derivative의 second term 계산
      sum_term_first += term;
      
      // second_derivative의 summation 부분 계산
      sum_term_second1 += term;
      sum_term_second2 += term * inv_sum_exp_eta_R_r;
    }
    
    double exp_eta_i = exp(eta[i]);
    double exp_2eta_i = exp_eta_i * exp_eta_i;
    
    first_derivative[i] = status[i] - exp_eta_i * sum_term_first;
    second_derivative[i] = -exp_eta_i * sum_term_second1 + exp_2eta_i * sum_term_second2;
  }
  
  return List::create(Named("first_derivative") = first_derivative,
                      Named("second_derivative") = second_derivative);
} 