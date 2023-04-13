#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat outer_times_rcpp(arma::vec a,
                           arma::vec b){

  arma::uword l_a = a.n_elem;
  arma::uword l_b = b.n_elem;
  arma::mat outer(l_a, l_b);

  for (size_t i = 0; i < l_a; ++i){

    for (size_t j = 0; j < l_b; ++j){

      outer(i, j) = a(i) * b(j);
    }
  }

  return outer;
}

// [[Rcpp::export]]
arma::umat outer_leq_rcpp_v1(arma::vec a,
                             arma::vec b){

  arma::uword l_a = a.n_elem;
  arma::uword l_b = b.n_elem;
  arma::umat outer(l_a, l_b);

  for (size_t i = 0; i < l_a; ++i){

    for (size_t j = 0; j < l_b; ++j){

      outer(i, j) = (a(i) <= b(j));
    }
  }

  return outer;
}

