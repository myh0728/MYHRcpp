#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// t_stop must be in ascending order

// [[Rcpp::export]]
arma::vec KME_rcpp(arma::vec t_stop,
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
arma::vec KME_w_rcpp(arma::vec t_stop,
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

  for (size_t k = 0; k < n_k; ++k){

    const auto xrow_k = x.row(k);

    for (size_t i = 0; i < n_n; ++i){

      arma::vec Dik_h = arma::vectorise(X.row(i) - xrow_k) / h;
      double Kik_h = arma::prod(
        (arma::abs(Dik_h) < 1) %
          arma::square(1.0 - arma::square(Dik_h)) * 15.0 / 16.0
      );

      for (size_t m = 0; m < n_t_event; ++m){

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

    for (size_t m = 0; m < n_t_event; ++m){

      if (Dhat(k, m) != 0){

        dLhat(k, m) = Nhat(k, m) / Dhat(k, m);
      }
    }
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

  for (size_t k = 0; k < n_k; ++k){

    const auto xrow_k = x.row(k);

    for (size_t i = 0; i < n_n; ++i){

      arma::vec Dik_h = arma::vectorise(X.row(i) - xrow_k) / h;
      double Kik_h = arma::prod(
        (arma::abs(Dik_h) < 1) %
          arma::square(1.0 - arma::square(Dik_h)) * 15.0 / 16.0
      ) * w(i);

      for (size_t m = 0; m < n_t_event; ++m){
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

    for (size_t m = 0; m < n_t_event; ++m){

      if (Dhat(k, m) != 0){

        dLhat(k, m) = Nhat(k, m) / Dhat(k, m);
      }
    }
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

  for (size_t k = 0; k < n_k; ++k){

    const auto xrow_k = x.row(k);

    for (size_t i = 0; i < n_n; ++i){

      arma::vec Dik_h_sq = arma::square(
        arma::vectorise(X.row(i) - xrow_k) / h
      );
      double Kik_h = arma::prod(
        (Dik_h_sq < 1) % (1.0 - Dik_h_sq * 3.0) %
          arma::square(1.0 - Dik_h_sq) * 105.0 / 64.0
      );

      for (size_t m = 0; m < n_t_event; ++m){

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

    for (size_t m = 0; m < n_t_event; ++m){

      if (Dhat(k, m) != 0){

        dLhat(k, m) = Nhat(k, m) / Dhat(k, m);
      }
    }
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

  for (size_t k = 0; k < n_k; ++k){

    const auto xrow_k = x.row(k);

    for (size_t i = 0; i < n_n; ++i){

      arma::vec Dik_h_sq = arma::square(
        arma::vectorise(X.row(i) - xrow_k) / h
      );
      double Kik_h = arma::prod(
        (Dik_h_sq < 1) % (1.0 - Dik_h_sq * 3.0) %
          arma::square(1.0 - Dik_h_sq) * 105.0 / 64.0
      ) * w(i);

      for (size_t m = 0; m < n_t_event; ++m){

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

    for (size_t m = 0; m < n_t_event; ++m){

      if (Dhat(k, m) != 0){

        dLhat(k, m) = Nhat(k, m) / Dhat(k, m);
      }
    }
  }

  return(dLhat);
}

// [[Rcpp::export]]
arma::mat SKME_KG_rcpp(arma::vec t_stop,
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

  for (size_t k = 0; k < n_k; ++k){

    const auto xrow_k = x.row(k);

    for (size_t i = 0; i < n_n; ++i){

      arma::vec Dik_h = arma::vectorise(X.row(i) - xrow_k) / h;
      double Kik_h = arma::prod(
        arma::exp(-arma::square(Dik_h) / 2.0) * 0.39894228
      );

      for (size_t m = 0; m < n_t_event; ++m){

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

    for (size_t m = 0; m < n_t_event; ++m){

      if (Dhat(k, m) != 0){

        dLhat(k, m) = Nhat(k, m) / Dhat(k, m);
      }
    }
  }

  return(dLhat);
}

// [[Rcpp::export]]
arma::mat SKME_KG_w_rcpp(arma::vec t_stop,
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

  for (size_t k = 0; k < n_k; ++k){

    const auto xrow_k = x.row(k);

    for (size_t i = 0; i < n_n; ++i){

      arma::vec Dik_h = arma::vectorise(X.row(i) - xrow_k) / h;
      double Kik_h = arma::prod(
        arma::exp(-arma::square(Dik_h) / 2.0) * 0.39894228
      ) * w(i);

      for (size_t m = 0; m < n_t_event; ++m){

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

    for (size_t m = 0; m < n_t_event; ++m){

      if (Dhat(k, m) != 0){

        dLhat(k, m) = Nhat(k, m) / Dhat(k, m);
      }
    }
  }

  return(dLhat);
}




