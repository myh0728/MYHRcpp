#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::vec testfunction(arma::vec a,
                       Rcpp::Function my_r_func){
  arma::vec s = Rcpp::as<arma::vec>(my_r_func(a));
  return(s);
}
