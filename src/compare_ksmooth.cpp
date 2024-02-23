#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::vec KDE_rcpp_kernel(arma::mat X,
                          arma::mat x,
                          Rcpp::Function K,
                          arma::vec h){

  arma::uword n_n = X.n_rows;
  arma::uword n_k = x.n_rows;
  arma::vec Dhat(n_k);

  for (size_t k = 0; k < n_k; ++k){

    for (size_t i = 0; i < n_n; ++i){

      arma::vec Dik_h = arma::vectorise(X.row(i) - x.row(k)) / h;
      double Kik_h = arma::prod(Rcpp::as<arma::vec>(K(Dik_h)) / h);
      Dhat(k) += Kik_h;
    }
  }

  Dhat /= n_n;
  return(Dhat);
}

// (i, j) runs over n_n times n_n iterations.
// [[Rcpp::export]]
arma::vec KDEcv_K2B_rcpp_o1(arma::mat X,
                            arma::vec h){

  arma::uword n_n = X.n_rows;
  arma::vec Dhat(n_n);

  for (size_t i = 0; i < n_n; ++i){

    const auto xrow_i = X.row(i);

    for (size_t j = 0; j < n_n; ++j){

      if (j != i){

        arma::vec Dji_h = arma::vectorise(X.row(j) - xrow_i) / h;
        double Kji_h = arma::prod(
          (arma::abs(Dji_h) < 1) %
            arma::square(1.0 - arma::square(Dji_h)) * 15.0 / 16.0 / h
        );
        Dhat(i) += Kji_h;
      }
    }
  }

  Dhat /= n_n - 1;
  return(Dhat);
}

// X.row(i) is not set as const.
// [[Rcpp::export]]
arma::vec KDEcv_K2B_rcpp_o2(arma::mat X,
                            arma::vec h){

  arma::uword n_n = X.n_rows;
  arma::vec Dhat(n_n);

  for (size_t i = 0; i < n_n; ++i){

    for (size_t j = 0; j < i; ++j){

      arma::vec Dji_h = arma::vectorise(X.row(j) - X.row(i)) / h;
      double Kji_h = arma::prod(
        (arma::abs(Dji_h) < 1) %
          arma::square(1.0 - arma::square(Dji_h)) * 15.0 / 16.0 / h
      );

      Dhat(i) += Kji_h;
      Dhat(j) += Kji_h;
    }
  }

  Dhat /= n_n - 1;
  return(Dhat);
}

// yhat.row is used.
// [[Rcpp::export]]
arma::mat NW_K2B_rcpp_o1(arma::mat X,
                         arma::mat Y,
                         arma::mat x,
                         arma::vec h){

  arma::uword n_n = X.n_rows;
  arma::uword n_k = x.n_rows;
  arma::uword n_m = Y.n_cols;
  arma::mat Nhat(n_k, n_m);
  arma::vec Dhat(n_k);
  arma::mat yhat(n_k, n_m);

  for (size_t k = 0; k < n_k; ++k){

    const auto xrow_k = x.row(k);

    for (size_t i = 0; i < n_n; ++i){

      arma::vec Dik_h = arma::vectorise(X.row(i) - xrow_k) / h;
      double Kik_h = arma::prod(
        (arma::abs(Dik_h) < 1) %
          arma::square(1.0 - arma::square(Dik_h)) * 15.0 / 16.0
      );
      Dhat(k) += Kik_h;

      for (size_t m = 0; m < n_m; ++m){

        Nhat(k, m) += Y(i, m) * Kik_h;
      }
    }

    if (Dhat(k) != 0){

      yhat.row(k) = Nhat.row(k) / Dhat(k);
    }
  }
  return(yhat);
}

// Nhat.row is used.
// [[Rcpp::export]]
arma::mat NW_K2B_rcpp_o2(arma::mat X,
                         arma::mat Y,
                         arma::mat x,
                         arma::vec h){

  arma::uword n_n = X.n_rows;
  arma::uword n_k = x.n_rows;
  arma::uword n_m = Y.n_cols;
  arma::mat Nhat(n_k, n_m);
  arma::vec Dhat(n_k);
  arma::mat yhat(n_k, n_m);

  for (size_t k = 0; k < n_k; ++k){

    const auto xrow_k = x.row(k);

    for (size_t i = 0; i < n_n; ++i){

      arma::vec Dik_h = arma::vectorise(X.row(i) - xrow_k) / h;
      double Kik_h = arma::prod(
        (arma::abs(Dik_h) < 1) %
          arma::square(1.0 - arma::square(Dik_h)) * 15.0 / 16.0
      );
      Dhat(k) += Kik_h;
      Nhat.row(k) += Y.row(i) * Kik_h;
    }

    if (Dhat(k) != 0){

      for (size_t m = 0; m < n_m; ++m){

        yhat(k, m) = Nhat(k, m) / Dhat(k);
      }
    }
  }

  return(yhat);
}

// Using const reference
// [[Rcpp::export]]
arma::mat NW_K2B_rcpp_n1(const arma::mat & X,
                         const arma::mat & Y,
                         const arma::mat & x,
                         const arma::vec & h){

  arma::uword n_n = X.n_rows;
  arma::uword n_k = x.n_rows;
  arma::uword n_m = Y.n_cols;
  arma::mat Nhat(n_k, n_m);
  arma::vec Dhat(n_k);
  arma::mat yhat(n_k, n_m);

  for (size_t k = 0; k < n_k; ++k){

    const auto xrow_k = x.row(k);

    for (size_t i = 0; i < n_n; ++i){

      arma::vec Dik_h = arma::vectorise(X.row(i) - xrow_k) / h;
      double Kik_h = arma::prod(
        (arma::abs(Dik_h) < 1) %
          arma::square(1.0 - arma::square(Dik_h)) * 15.0 / 16.0
      );
      Dhat(k) += Kik_h;

      for (size_t m = 0; m < n_m; ++m){

        Nhat(k, m) += Y(i, m) * Kik_h;
      }
    }

    if (Dhat(k) != 0){

      for (size_t m = 0; m < n_m; ++m){

        yhat(k, m) = Nhat(k, m) / Dhat(k);
      }
    }
  }

  return(yhat);
}

// Using arma::all
// [[Rcpp::export]]
arma::mat NWD_K2B_rcpp(const arma::mat& X,
                       const arma::mat& Y,
                       const arma::mat& x,
                       const arma::mat& y,
                       const arma::vec& h){

  arma::uword n_n = X.n_rows;
  arma::uword n_k = x.n_rows;
  arma::uword n_m = y.n_rows;
  arma::mat Nhat(n_k, n_m);
  arma::vec Dhat(n_k);
  arma::mat Fhat(n_k, n_m);

  for (size_t i = 0; i < n_n; ++i){

    for (size_t k = 0; k < n_k; ++k){

      arma::vec Dik_h = arma::vectorise(X.row(i) - x.row(k)) / h;
      double Kik_h = arma::prod(
        (arma::abs(Dik_h) < 1) %
          arma::square(1.0 - arma::square(Dik_h)) * 15.0 / 16.0
      );
      Dhat(k) += Kik_h;

      for (size_t m = 0; m < n_m; ++m){

        if (arma::all(Y.row(i) <= y.row(m)) == 1){

          Nhat(k, m) += Kik_h;
        }
      }
    }
  }

  for (size_t k = 0; k < n_k; ++k){

    if (Dhat(k) != 0){

      for (size_t m = 0; m < n_m; ++m){

        Fhat(k, m) = Nhat(k, m) / Dhat(k);
      }
    }
  }

  return(Fhat);
}

// Using arma::prod
// [[Rcpp::export]]
arma::mat NWD_K2B_rcpp_v1(arma::mat X,
                          arma::mat Y,
                          arma::mat x,
                          arma::mat y,
                          arma::vec h){

  arma::uword n_n = X.n_rows;
  arma::uword n_k = x.n_rows;
  arma::uword n_m = y.n_rows;
  arma::mat Nhat(n_k, n_m);
  arma::vec Dhat(n_k);
  arma::mat Fhat(n_k, n_m);

  for (size_t k = 0; k < n_k; ++k){

    const auto xrow_k = x.row(k);

    for (size_t i = 0; i < n_n; ++i){

      arma::vec Dik_h = arma::vectorise(X.row(i) - xrow_k) / h;
      double Kik_h = arma::prod(
        (arma::abs(Dik_h) < 1) %
          arma::square(1.0 - arma::square(Dik_h)) * 15.0 / 16.0
      );
      Dhat(k) += Kik_h;

      for (size_t m = 0; m < n_m; ++m){

        Nhat(k, m) += arma::prod(Y.row(i) <= y.row(m)) * Kik_h;
      }
    }

    if (Dhat(k) != 0){

      for (size_t m = 0; m < n_m; ++m){

        Fhat(k, m) = Nhat(k, m) / Dhat(k);
      }
    }
  }

  return(Fhat);
}

// [[Rcpp::export]]
double CVMNW_K2B_rcpp_n1(arma::mat X,
                         arma::mat Y,
                         arma::vec h,
                         arma::vec p_Y){

  arma::uword n_n = X.n_rows;
  arma::uword n_m = Y.n_cols;
  double cv = 0.0;

  for (size_t i = 0; i < n_n; ++i){

    arma::vec Nhat_i(n_m);
    double Dhat_i = 0;
    const auto xrow_i = X.row(i);

    for (size_t j = 0; j < n_n; ++j){

      if (j != i){

        arma::vec Dji_h = arma::vectorise(X.row(j) - xrow_i) / h;
        double Kji_h = arma::prod(
          (arma::abs(Dji_h) < 1.0) %
            arma::square(1.0 - arma::square(Dji_h)) * 15.0 / 16.0
        );
        Dhat_i += Kji_h;
        Nhat_i += Y.row(j) * Kji_h;
      }
    }

    if (Dhat_i != 0){

      cv += arma::sum(arma::pow(Y.row(i) - Nhat_i / Dhat_i, 2) % p_Y);

    }else{

      cv += arma::sum(arma::pow(Y.row(i), 2) % p_Y);
    }
  }

  cv /= n_n;

  return cv;
}

// [[Rcpp::export]]
double CVDNWuniY_K2B_rcpp_n1(arma::mat X,
                             arma::vec Y,
                             arma::vec h,
                             arma::vec y,
                             arma::vec p_y){

  arma::uword n_n = X.n_rows;
  arma::uword n_m = y.n_elem;
  double cv = 0;

  for (size_t i = 0; i < n_n; ++i){

    arma::vec Nhat_i(n_m);
    double Dhat_i = 0;
    const auto xrow_i = X.row(i);

    for (size_t j = 0; j < n_n; ++j){

      if (j != i){

        arma::vec Dji_h = arma::vectorise(X.row(j) - xrow_i) / h;
        double Kji_h = arma::prod(
          (arma::abs(Dji_h) < 1.0) %
            arma::square(1.0 - arma::square(Dji_h)) * 15.0 / 16.0
        );
        Dhat_i += Kji_h;
        Nhat_i.elem(arma::find(y >= Y(j))) += Kji_h;
      }
    }

    if (Dhat_i != 0){

      // arma::vec Fhat_i = arma::clamp(Nhat_i / Dhat_i, 0.0, 1.0);
      cv += arma::sum(arma::pow((y >= Y(i)) - Nhat_i / Dhat_i, 2) % p_y);

    }else{

      cv += arma::sum((y >= Y(i)) % p_y);
    }
  }

  cv /= n_n;

  return(cv);
}


















