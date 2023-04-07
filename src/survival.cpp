#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::vec KME_rcpp(arma::vec t_stop,
                   arma::uvec is_event,
                   arma::vec t_event){
  arma::uword n_n = t_stop.n_elem;
  arma::uword n_t_event = t_event.n_elem;
  arma::vec Nhat(n_t_event);
  arma::vec Dhat(n_t_event);
  for (arma::uword k = 0; k < n_t_event; ++k){
    for (arma::uword i = n_n; i > 0; --i){
      if (t_stop(i-1) >= t_event(k)){
        Dhat(k) += 1;
        if (t_stop(i-1) == t_event(k) and is_event(i-1) == 1){
          Nhat(k) += 1;
        }
      }else{
        break;
      }
    }
  }
  arma::vec dLhat = Nhat/Dhat;
  return(dLhat);
}

// [[Rcpp::export]]
arma::vec KME_w_rcpp(arma::vec t_stop,
                   arma::uvec is_event,
                   arma::vec t_event,
                   arma::vec w){
  arma::uword n_n = t_stop.n_elem;
  arma::uword n_t_event = t_event.n_elem;
  arma::vec Nhat(n_t_event);
  arma::vec Dhat(n_t_event);
  for (arma::uword k = 0; k < n_t_event; ++k){
    for (arma::uword i = n_n; i > 0; --i){
      if (t_stop(i-1) >= t_event(k)){
        Dhat(k) += w(i-1);
        if (t_stop(i-1) == t_event(k) and is_event(i-1) == 1){
          Nhat(k) += w(i-1);
        }
      }else{
        break;
      }
    }
  }
  arma::vec dLhat = Nhat/Dhat;
  return(dLhat);
}

// [[Rcpp::export]]
arma::mat SKME_K2B_rcpp(arma::vec t_stop,
                        arma::uvec is_event,
                        arma::vec t_event,
                        arma::mat X,
                        arma::mat x,
                        arma::vec h){
  arma::uword n_n = t_stop.n_elem;
  arma::uword n_t_event = t_event.n_elem;
  arma::uword n_k = x.n_rows;
  arma::mat Nhat(n_k, n_t_event);
  arma::mat Dhat(n_k, n_t_event);
  arma::mat dLhat(n_k, n_t_event);
  for (arma::uword k = 0; k < n_k; ++k){
    for (arma::uword i = 0; i < n_n; ++i){
      arma::vec Dik_h = vectorise(X.row(i)-x.row(k))/h;
      double Kik_h = prod((abs(Dik_h) < 1)%pow(1-pow(Dik_h, 2), 2)*15/16);
      for (arma::uword m = 0; m < n_t_event; ++m){
        if (t_event(m) <= t_stop(i)){
          Dhat(k, m) += Kik_h;
          if (t_event(m) == t_stop(i) and is_event(i) == 1){
            Nhat(k, m) += Kik_h;
          }
        }else{
          break;
        }
      }
    }
    dLhat.row(k) = Nhat.row(k)/Dhat.row(k);
  }
  return(dLhat);
}

// [[Rcpp::export]]
arma::mat SKME_K2B_w_rcpp(arma::vec t_stop,
                          arma::uvec is_event,
                          arma::vec t_event,
                          arma::mat X,
                          arma::mat x,
                          arma::vec h,
                          arma::vec w){
  arma::uword n_n = t_stop.n_elem;
  arma::uword n_t_event = t_event.n_elem;
  arma::uword n_k = x.n_rows;
  arma::mat Nhat(n_k, n_t_event);
  arma::mat Dhat(n_k, n_t_event);
  arma::mat dLhat(n_k, n_t_event);
  for (arma::uword k = 0; k < n_k; ++k){
    for (arma::uword i = 0; i < n_n; ++i){
      arma::vec Dik_h = vectorise(X.row(i)-x.row(k))/h;
      double Kik_h = prod((abs(Dik_h) < 1)%pow(1-pow(Dik_h, 2), 2)*15/16)*w(i);
      for (arma::uword m = 0; m < n_t_event; ++m){
        if (t_event(m) <= t_stop(i)){
          Dhat(k, m) += Kik_h;
          if (t_event(m) == t_stop(i) and is_event(i) == 1){
            Nhat(k, m) += Kik_h;
          }
        }else{
          break;
        }
      }
    }
    dLhat.row(k) = Nhat.row(k)/Dhat.row(k);
  }
  return(dLhat);
}

// [[Rcpp::export]]
arma::mat SKME_K4B_rcpp(arma::vec t_stop,
                        arma::uvec is_event,
                        arma::vec t_event,
                        arma::mat X,
                        arma::mat x,
                        arma::vec h){
  arma::uword n_n = t_stop.n_elem;
  arma::uword n_t_event = t_event.n_elem;
  arma::uword n_k = x.n_rows;
  arma::mat Nhat(n_k, n_t_event);
  arma::mat Dhat(n_k, n_t_event);
  arma::mat dLhat(n_k, n_t_event);
  for (arma::uword k = 0; k < n_k; ++k){
    for (arma::uword i = 0; i < n_n; ++i){
      arma::vec Dik_h = vectorise(X.row(i)-x.row(k))/h;
      double Kik_h = prod((abs(Dik_h) < 1)%(1-3*pow(Dik_h, 2))%
                          pow(1-pow(Dik_h, 2), 2)*105/64);
      for (arma::uword m = 0; m < n_t_event; ++m){
        if (t_event(m) <= t_stop(i)){
          Dhat(k, m) += Kik_h;
          if (t_event(m) == t_stop(i) and is_event(i) == 1){
            Nhat(k, m) += Kik_h;
          }
        }else{
          break;
        }
      }
    }
    dLhat.row(k) = Nhat.row(k)/Dhat.row(k);
  }
  return(dLhat);
}

// [[Rcpp::export]]
arma::mat SKME_K4B_w_rcpp(arma::vec t_stop,
                          arma::uvec is_event,
                          arma::vec t_event,
                          arma::mat X,
                          arma::mat x,
                          arma::vec h,
                          arma::vec w){
  arma::uword n_n = t_stop.n_elem;
  arma::uword n_t_event = t_event.n_elem;
  arma::uword n_k = x.n_rows;
  arma::mat Nhat(n_k, n_t_event);
  arma::mat Dhat(n_k, n_t_event);
  arma::mat dLhat(n_k, n_t_event);
  for (arma::uword k = 0; k < n_k; ++k){
    for (arma::uword i = 0; i < n_n; ++i){
      arma::vec Dik_h = vectorise(X.row(i)-x.row(k))/h;
      double Kik_h = prod((abs(Dik_h) < 1)%(1-3*pow(Dik_h, 2))%
                          pow(1-pow(Dik_h, 2), 2)*105/64)*w(i);
      for (arma::uword m = 0; m < n_t_event; ++m){
        if (t_event(m) <= t_stop(i)){
          Dhat(k, m) += Kik_h;
          if (t_event(m) == t_stop(i) and is_event(i) == 1){
            Nhat(k, m) += Kik_h;
          }
        }else{
          break;
        }
      }
    }
    dLhat.row(k) = Nhat.row(k)/Dhat.row(k);
  }
  return(dLhat);
}
