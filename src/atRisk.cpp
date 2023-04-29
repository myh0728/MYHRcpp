#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat atRisk_integral_rcpp(arma::mat integrand,
                               arma::vec t_start,
                               arma::vec t_stop,
                               arma::vec t_event){

  arma::uword n_obs = t_start.n_elem;
  arma::uword n_t_event = t_event.n_elem;
  arma::uword n_p = integrand.n_cols;
  arma::mat result(n_obs, n_p);

  for (size_t i = 0; i < n_t_event; ++i){

    for (size_t j = 0; j < n_obs; ++j){

      if (t_start[j] < t_event[i] and t_stop[j] >= t_event[i]){

        result.row(j) = result.row(j) + integrand.row(i);
      }
    }
  }
  return result;
}

// [[Rcpp::export]]
arma::mat sum_atRisk_rcpp(arma::mat summand,
                          arma::vec t_start,
                          arma::vec t_stop,
                          arma::vec t_event){

  arma::uword n_obs = t_start.n_elem;
  arma::uword n_t_event = t_event.n_elem;
  arma::uword n_p = summand.n_cols;
  arma::mat result(n_t_event, n_p);

  for (size_t i = 0; i < n_t_event; ++i){

    for (size_t j = 0; j < n_obs; ++j){

      if (t_start[j] < t_event[i] and t_stop[j] >= t_event[i]){

        result.row(i) = result.row(i) + summand.row(j);
      }
    }
  }
  return result;
}






