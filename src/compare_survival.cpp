#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// t_stop must be in ascending order

// [[Rcpp::export]]
arma::vec KME_rcpp_n1(arma::vec t_stop,
                      arma::uvec is_event,
                      arma::vec t_event){

  arma::uword n_n = t_stop.n_elem;
  arma::uword n_t_event = t_event.n_elem;
  arma::vec Nhat(n_t_event);
  arma::vec Dhat(n_t_event);

  for (size_t k = 0; k < n_t_event; ++k){

    for (size_t i = n_n; i > 0; --i){

      if (t_stop(i - 1) >= t_event(k)){

        Dhat(k) += 1;

        if (t_stop(i - 1) == t_event(k) and is_event(i - 1) == 1){

          Nhat(k) += 1;
        }
      }else{

        break;
      }
    }
  }

  arma::vec dLhat = Nhat / Dhat;
  return(dLhat);
}

// [[Rcpp::export]]
arma::vec KME_w_rcpp_n1(arma::vec t_stop,
                        arma::uvec is_event,
                        arma::vec t_event,
                        arma::vec w){

  arma::uword n_n = t_stop.n_elem;
  arma::uword n_t_event = t_event.n_elem;
  arma::vec Nhat(n_t_event);
  arma::vec Dhat(n_t_event);

  for (size_t k = 0; k < n_t_event; ++k){

    for (size_t i = n_n; i > 0; --i){

      if (t_stop(i - 1) >= t_event(k)){

        Dhat(k) += w(i - 1);

        if (t_stop(i - 1) == t_event(k) and is_event(i - 1) == 1){

          Nhat(k) += w(i - 1);
        }
      }else{

        break;
      }
    }
  }

  arma::vec dLhat = Nhat/Dhat;
  return(dLhat);
}
