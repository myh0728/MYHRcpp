#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// kernel density estimation

// Using K2_Biweight kernel

// [[Rcpp::export]]
arma::vec KDE_K2B_rcpp(arma::mat X,
                       arma::mat x,
                       arma::vec h){

  arma::uword n_n = X.n_rows;
  arma::uword n_k = x.n_rows;
  arma::vec Dhat(n_k);

  for (size_t k = 0; k < n_k; ++k){

    const auto xrow_k = x.row(k);

    for (size_t i = 0; i < n_n; ++i){

      arma::vec Dik_h = arma::vectorise(X.row(i) - xrow_k) / h;
      double Kik_h = arma::prod(
        (arma::abs(Dik_h) < 1) %
          arma::square(1.0 - arma::square(Dik_h)) * 15.0 / 16.0 / h
      );
      Dhat(k) += Kik_h;
    }
  }

  Dhat /= n_n;
  return(Dhat);
}

// [[Rcpp::export]]
arma::vec KDE_K2B_w_rcpp(arma::mat X,
                         arma::mat x,
                         arma::vec h,
                         arma::vec w){

  arma::uword n_n = X.n_rows;
  arma::uword n_k = x.n_rows;
  arma::vec Dhat(n_k);

  for (size_t k = 0; k < n_k; ++k){

    const auto xrow_k = x.row(k);

    for (size_t i = 0; i < n_n; ++i){

      arma::vec Dik_h = arma::vectorise(X.row(i) - xrow_k) / h;
      double Kik_h = arma::prod(
        (arma::abs(Dik_h) < 1) %
          arma::square(1.0 - arma::square(Dik_h)) * 15.0 / 16.0 / h
      );
      Dhat(k) += Kik_h * w(i);
    }
  }

  Dhat /= n_n;
  return(Dhat);
}

// [[Rcpp::export]]
arma::vec KDEcv_K2B_rcpp(arma::mat X,
                         arma::vec h){

  arma::uword n_n = X.n_rows;
  arma::vec Dhat(n_n);

  for (size_t i = 0; i < n_n; ++i){

    const auto xrow_i = X.row(i);

    for (size_t j = 0; j < i; ++j){

      arma::vec Dji_h = arma::vectorise(X.row(j) - xrow_i) / h;
      double Kji_h = arma::prod(
        (arma::abs(Dji_h) < 1) %
          arma::square(1.0 - arma::square(Dji_h)) * 15.0 / 16.0 / h
      );

      Dhat(i) += Kji_h;
      Dhat(j) += Kji_h;
    }
  }

  Dhat /= n_n-1;
  return(Dhat);
}

// [[Rcpp::export]]
arma::vec KDEcv_K2B_w_rcpp(arma::mat X,
                           arma::vec h,
                           arma::vec w){

  arma::uword n_n = X.n_rows;
  arma::vec Dhat(n_n);

  for (size_t i = 0; i < n_n; ++i){

    const auto xrow_i = X.row(i);

    for (size_t j = 0; j < i; ++j){

      arma::vec Dji_h = arma::vectorise(X.row(j) - xrow_i) / h;
      double Kji_h = arma::prod(
        (arma::abs(Dji_h) < 1) %
          arma::square(1.0 - arma::square(Dji_h)) * 15.0 / 16.0 / h
      );

      Dhat(i) += Kji_h * w(j);
      Dhat(j) += Kji_h * w(i);
    }
  }

  Dhat /= n_n-1;
  return(Dhat);
}

// Using K4_Biweight kernel

// [[Rcpp::export]]
arma::vec KDE_K4B_rcpp(arma::mat X,
                       arma::mat x,
                       arma::vec h){

  arma::uword n_n = X.n_rows;
  arma::uword n_k = x.n_rows;
  arma::vec Dhat(n_k);

  for (size_t k = 0; k < n_k; ++k){

    const auto xrow_k = x.row(k);

    for (size_t i = 0; i < n_n; ++i){

      arma::vec Dik_h_sq = arma::square(
        arma::vectorise(X.row(i) - xrow_k) / h
      );
      double Kik_h = arma::prod(
        (Dik_h_sq < 1) % (1.0 - Dik_h_sq * 3.0) %
          arma::square(1.0 - Dik_h_sq) * 105.0 / 64.0 / h
      );
      Dhat(k) += Kik_h;
    }
  }

  Dhat /= n_n;
  Dhat(find(Dhat < 0)).fill(0);
  return(Dhat);
}

// [[Rcpp::export]]
arma::vec KDE_K4B_w_rcpp(arma::mat X,
                         arma::mat x,
                         arma::vec h,
                         arma::vec w){

  arma::uword n_n = X.n_rows;
  arma::uword n_k = x.n_rows;
  arma::vec Dhat(n_k);

  for (size_t k = 0; k < n_k; ++k){

    const auto xrow_k = x.row(k);

    for (size_t i = 0; i < n_n; ++i){

      arma::vec Dik_h_sq = arma::square(
        arma::vectorise(X.row(i) - xrow_k) / h
      );
      double Kik_h = arma::prod(
        (Dik_h_sq < 1) % (1.0 - Dik_h_sq * 3.0) %
          arma::square(1.0 - Dik_h_sq) * 105.0 / 64.0 / h
      );
      Dhat(k) += Kik_h * w(i);
    }
  }

  Dhat /= n_n;
  Dhat(find(Dhat < 0)).fill(0);
  return(Dhat);
}

// [[Rcpp::export]]
arma::vec KDEcv_K4B_rcpp(arma::mat X,
                         arma::vec h){

  arma::uword n_n = X.n_rows;
  arma::vec Dhat(n_n);

  for (size_t i = 0; i < n_n; ++i){

    const auto xrow_i = X.row(i);

    for (size_t j = 0; j < i; ++j){

      arma::vec Dji_h_sq = arma::square(
        arma::vectorise(X.row(j) - xrow_i) / h
      );
      double Kji_h = arma::prod(
        (Dji_h_sq < 1) % (1 - Dji_h_sq * 3.0) %
          arma::square(1.0 - Dji_h_sq) * 105.0 / 64.0 / h
      );

      Dhat(i) += Kji_h;
      Dhat(j) += Kji_h;
    }
  }

  Dhat /= n_n-1;
  Dhat(find(Dhat < 0)).fill(0);
  return(Dhat);
}

// [[Rcpp::export]]
arma::vec KDEcv_K4B_w_rcpp(arma::mat X,
                           arma::vec h,
                           arma::vec w){

  arma::uword n_n = X.n_rows;
  arma::vec Dhat(n_n);

  for (size_t i = 0; i < n_n; ++i){

    const auto xrow_i = X.row(i);

    for (size_t j = 0; j < i; ++j){

      arma::vec Dji_h_sq = arma::square(
        arma::vectorise(X.row(j) - xrow_i) / h
      );
      double Kji_h = arma::prod(
        (Dji_h_sq < 1) % (1 - Dji_h_sq * 3.0) %
          arma::square(1.0 - Dji_h_sq) * 105.0 / 64.0 / h
      );

      Dhat(i) += Kji_h * w(j);
      Dhat(j) += Kji_h * w(i);
    }
  }

  Dhat /= n_n-1;
  Dhat(find(Dhat < 0)).fill(0);
  return(Dhat);
}

// Using Gaussian kernel

// [[Rcpp::export]]
arma::vec KDE_KG_rcpp(arma::mat X,
                      arma::mat x,
                      arma::vec h){

  arma::uword n_n = X.n_rows;
  arma::uword n_k = x.n_rows;
  arma::vec Dhat(n_k);

  for (size_t k = 0; k < n_k; ++k){

    const auto xrow_k = x.row(k);

    for (size_t i = 0; i < n_n; ++i){

      arma::vec Dik_h = arma::vectorise(X.row(i) - xrow_k) / h;
      double Kik_h = arma::prod(
        arma::exp(-arma::square(Dik_h) / 2.0) * 0.39894228 / h
      );
      Dhat(k) += Kik_h;
    }
  }

  Dhat /= n_n;
  return(Dhat);
}

// [[Rcpp::export]]
arma::vec KDE_KG_w_rcpp(arma::mat X,
                        arma::mat x,
                        arma::vec h,
                        arma::vec w){

  arma::uword n_n = X.n_rows;
  arma::uword n_k = x.n_rows;
  arma::vec Dhat(n_k);

  for (size_t k = 0; k < n_k; ++k){

    const auto xrow_k = x.row(k);

    for (size_t i = 0; i < n_n; ++i){

      arma::vec Dik_h = arma::vectorise(X.row(i) - xrow_k) / h;
      double Kik_h = arma::prod(
        arma::exp(-arma::square(Dik_h) / 2.0) * 0.39894228 / h
      );
      Dhat(k) += Kik_h * w(i);
    }
  }

  Dhat /= n_n;
  return(Dhat);
}

// [[Rcpp::export]]
arma::vec KDEcv_KG_rcpp(arma::mat X,
                        arma::vec h){

  arma::uword n_n = X.n_rows;
  arma::vec Dhat(n_n);

  for (size_t i = 0; i < n_n; ++i){

    const auto xrow_i = X.row(i);

    for (size_t j = 0; j < i; ++j){

      arma::vec Dji_h = arma::vectorise(X.row(j) - xrow_i) / h;
      double Kji_h = arma::prod(
        arma::exp(-arma::square(Dji_h) / 2.0) * 0.39894228 / h
      );

      Dhat(i) += Kji_h;
      Dhat(j) += Kji_h;
    }
  }

  Dhat /= n_n-1;
  return(Dhat);
}

// [[Rcpp::export]]
arma::vec KDEcv_KG_w_rcpp(arma::mat X,
                          arma::vec h,
                          arma::vec w){

  arma::uword n_n = X.n_rows;
  arma::vec Dhat(n_n);

  for (size_t i = 0; i < n_n; ++i){

    const auto xrow_i = X.row(i);

    for (size_t j = 0; j < i; ++j){

      arma::vec Dji_h = arma::vectorise(X.row(j) - xrow_i) / h;
      double Kji_h = arma::prod(
        arma::exp(-arma::square(Dji_h) / 2.0) * 0.39894228 / h
      );

      Dhat(i) += Kji_h * w(j);
      Dhat(j) += Kji_h * w(i);
    }
  }

  Dhat /= n_n-1;
  return(Dhat);
}

// Nadaraya-Watson estimator

// Using K2_Biweight kernel

// [[Rcpp::export]]
arma::mat NW_K2B_rcpp(arma::mat X,
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

      for (size_t m = 0; m < n_m; ++m){

        yhat(k, m) = Nhat(k, m) / Dhat(k);
      }
    }
  }

  return(yhat);
}

// [[Rcpp::export]]
arma::mat NW_K2B_w_rcpp(arma::mat X,
                        arma::mat Y,
                        arma::mat x,
                        arma::vec h,
                        arma::vec w){

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
      ) * w(i);
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

// [[Rcpp::export]]
arma::mat NWcv_K2B_rcpp(arma::mat X,
                        arma::mat Y,
                        arma::vec h){

  arma::uword n_n = X.n_rows;
  arma::uword n_m = Y.n_cols;
  arma::mat Nhat(n_n, n_m);
  arma::vec Dhat(n_n);
  arma::mat Yhat(n_n, n_m);

  for (size_t i = 0; i < n_n; ++i){

    const auto xrow_i = X.row(i);

    for (size_t j = 0; j < i; ++j)
    {

      arma::vec Dji_h = arma::vectorise(X.row(j) - xrow_i) / h;
      double Kji_h = arma::prod(
        (arma::abs(Dji_h) < 1.0) %
          arma::square(1.0 - arma::square(Dji_h)) * 15.0 / 16.0
      );

      Dhat(i) += Kji_h;
      Dhat(j) += Kji_h;

      for (size_t m = 0; m < n_m; ++m){

        Nhat(i, m) += Y(j, m) * Kji_h;
        Nhat(j, m) += Y(i, m) * Kji_h;
      }
    }
  }

  for (size_t i = 0; i < n_n; ++i){

    if (Dhat(i) != 0){

      for (size_t m = 0; m < n_m; ++m){

        Yhat(i, m) = Nhat(i, m) / Dhat(i);
      }
    }
  }

  return(Yhat);
}

// [[Rcpp::export]]
arma::mat NWcv_K2B_w_rcpp(arma::mat X,
                          arma::mat Y,
                          arma::vec h,
                          arma::vec w){

  arma::uword n_n = X.n_rows;
  arma::uword n_m = Y.n_cols;
  arma::mat Nhat(n_n, n_m);
  arma::vec Dhat(n_n);
  arma::mat Yhat(n_n, n_m);

  for (size_t i = 0; i < n_n; ++i){

    const auto xrow_i = X.row(i);

    for (size_t j = 0; j < i; ++j){

      arma::vec Dji_h = arma::vectorise(X.row(j) - xrow_i) / h;
      double Kji_h = arma::prod(
        (arma::abs(Dji_h) < 1.0) %
          arma::square(1.0 - arma::square(Dji_h)) * 15.0 / 16.0
      );

      Dhat(i) += Kji_h * w(j);
      Dhat(j) += Kji_h * w(i);

      for (size_t m = 0; m < n_m; ++m){

        Nhat(i, m) += Y(j, m) * Kji_h * w(j);
        Nhat(j, m) += Y(i, m) * Kji_h * w(i);
      }
    }
  }

  for (size_t i = 0; i < n_n; ++i){

    if (Dhat(i) != 0){

      for (size_t m = 0; m < n_m; ++m){

        Yhat(i, m) = Nhat(i, m) / Dhat(i);
      }
    }
  }

  return(Yhat);
}

// Using K4_Biweight kernel

// [[Rcpp::export]]
arma::mat NW_K4B_rcpp(arma::mat X,
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

      arma::vec Dik_h_sq = arma::square(
        arma::vectorise(X.row(i) - xrow_k) / h
      );
      double Kik_h = arma::prod(
        (Dik_h_sq < 1) % (1.0 - Dik_h_sq * 3.0) %
          arma::square(1.0 - Dik_h_sq) * 105.0 / 64.0
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

// [[Rcpp::export]]
arma::mat NW_K4B_w_rcpp(arma::mat X,
                        arma::mat Y,
                        arma::mat x,
                        arma::vec h,
                        arma::vec w){

  arma::uword n_n = X.n_rows;
  arma::uword n_k = x.n_rows;
  arma::uword n_m = Y.n_cols;
  arma::mat Nhat(n_k, n_m);
  arma::vec Dhat(n_k);
  arma::mat yhat(n_k, n_m);

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

// [[Rcpp::export]]
arma::mat NWcv_K4B_rcpp(arma::mat X,
                        arma::mat Y,
                        arma::vec h){

  arma::uword n_n = X.n_rows;
  arma::uword n_m = Y.n_cols;
  arma::mat Nhat(n_n, n_m);
  arma::vec Dhat(n_n);
  arma::mat Yhat(n_n, n_m);

  for (size_t i = 0; i < n_n; ++i){

    const auto xrow_i = X.row(i);

    for (size_t j = 0; j < i; ++j){

      arma::vec Dji_h_sq = arma::square(
        arma::vectorise(X.row(j) - xrow_i) / h
      );
      double Kji_h = arma::prod(
        (Dji_h_sq < 1) % (1.0 - Dji_h_sq * 3.0) %
          arma::square(1.0 - Dji_h_sq) * 105.0 / 64.0
      );

      Dhat(i) += Kji_h;
      Dhat(j) += Kji_h;

      for (size_t m = 0; m < n_m; ++m){

        Nhat(i, m) += Y(j, m) * Kji_h;
        Nhat(j, m) += Y(i, m) * Kji_h;
      }
    }
  }

  for (size_t i = 0; i < n_n; ++i){

    if (Dhat(i) != 0){

      for (size_t m = 0; m < n_m; ++m){

        Yhat(i, m) = Nhat(i, m) / Dhat(i);
      }
    }
  }

  return(Yhat);
}

// [[Rcpp::export]]
arma::mat NWcv_K4B_w_rcpp(arma::mat X,
                          arma::mat Y,
                          arma::vec h,
                          arma::vec w){

  arma::uword n_n = X.n_rows;
  arma::uword n_m = Y.n_cols;
  arma::mat Nhat(n_n, n_m);
  arma::vec Dhat(n_n);
  arma::mat Yhat(n_n, n_m);

  for (size_t i = 0; i < n_n; ++i){

    const auto xrow_i = X.row(i);

    for (size_t j = 0; j < i; ++j){

      arma::vec Dji_h_sq = arma::square(
        arma::vectorise(X.row(j) - xrow_i) / h
      );
      double Kji_h = arma::prod(
        (Dji_h_sq < 1) % (1.0 - Dji_h_sq * 3.0) %
          arma::square(1.0 - Dji_h_sq) * 105.0 / 64.0
      );

      Dhat(i) += Kji_h * w(j);
      Dhat(j) += Kji_h * w(i);

      for (size_t m = 0; m < n_m; ++m){

        Nhat(i, m) += Y(j, m) * Kji_h * w(j);
        Nhat(j, m) += Y(i, m) * Kji_h * w(i);
      }
    }
  }

  for (size_t i = 0; i < n_n; ++i){

    if (Dhat(i) != 0){

      for (size_t m = 0; m < n_m; ++m){

        Yhat(i, m) = Nhat(i, m) / Dhat(i);
      }
    }
  }

  return(Yhat);
}

// Using Gaussian kernel

// [[Rcpp::export]]
arma::mat NW_KG_rcpp(arma::mat X,
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
        arma::exp(-arma::square(Dik_h) / 2.0) * 0.39894228
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

// [[Rcpp::export]]
arma::mat NW_KG_w_rcpp(arma::mat X,
                       arma::mat Y,
                       arma::mat x,
                       arma::vec h,
                       arma::vec w){

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
        arma::exp(-arma::square(Dik_h) / 2.0) * 0.39894228
      ) * w(i);
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

// [[Rcpp::export]]
arma::mat NWcv_KG_rcpp(arma::mat X,
                       arma::mat Y,
                       arma::vec h){

  arma::uword n_n = X.n_rows;
  arma::uword n_m = Y.n_cols;
  arma::mat Nhat(n_n, n_m);
  arma::vec Dhat(n_n);
  arma::mat Yhat(n_n, n_m);

  for (size_t i = 0; i < n_n; ++i){

    const auto xrow_i = X.row(i);

    for (size_t j = 0; j < i; ++j){

      arma::vec Dji_h = arma::vectorise(X.row(j) - xrow_i) / h;
      double Kji_h = arma::prod(
        arma::exp(-arma::square(Dji_h) / 2.0) * 0.39894228
      );

      Dhat(i) += Kji_h;
      Dhat(j) += Kji_h;

      for (size_t m = 0; m < n_m; ++m){

        Nhat(i, m) += Y(j, m) * Kji_h;
        Nhat(j, m) += Y(i, m) * Kji_h;
      }
    }
  }

  for (size_t i = 0; i < n_n; ++i){

    if (Dhat(i) != 0){

      for (size_t m = 0; m < n_m; ++m){

        Yhat(i, m) = Nhat(i, m) / Dhat(i);
      }
    }
  }

  return(Yhat);
}

// [[Rcpp::export]]
arma::mat NWcv_KG_w_rcpp(arma::mat X,
                         arma::mat Y,
                         arma::vec h,
                         arma::vec w){

  arma::uword n_n = X.n_rows;
  arma::uword n_m = Y.n_cols;
  arma::mat Nhat(n_n, n_m);
  arma::vec Dhat(n_n);
  arma::mat Yhat(n_n, n_m);

  for (size_t i = 0; i < n_n; ++i){

    const auto xrow_i = X.row(i);

    for (size_t j = 0; j < i; ++j){

      arma::vec Dji_h = arma::vectorise(X.row(j) - xrow_i) / h;
      double Kji_h = arma::prod(
        arma::exp(-arma::square(Dji_h) / 2.0) * 0.39894228
      );

      Dhat(i) += Kji_h * w(j);
      Dhat(j) += Kji_h * w(i);

      for (size_t m = 0; m < n_m; ++m){

        Nhat(i, m) += Y(j, m) * Kji_h * w(j);
        Nhat(j, m) += Y(i, m) * Kji_h * w(i);
      }
    }
  }

  for (size_t i = 0; i < n_n; ++i){

    if (Dhat(i) != 0){

      for (size_t m = 0; m < n_m; ++m){

        Yhat(i, m) = Nhat(i, m) / Dhat(i);
      }
    }
  }

  return(Yhat);
}

// cross-validation criterion for Nadaraya-Watson estimator

// conditional mean

// Using K2_Biweight kernel

// [[Rcpp::export]]
double CVMNW_K2B_rcpp(arma::mat X,
                      arma::mat Y,
                      arma::vec h,
                      arma::vec p_Y){

  arma::uword n_n = X.n_rows;
  arma::uword n_m = Y.n_cols;
  arma::mat Nhat(n_n, n_m);
  arma::vec Dhat(n_n);
  double cv = 0.0;

  for (size_t i = 0; i < n_n; ++i){

    const auto xrow_i = X.row(i);

    for (size_t j = 0; j < i; ++j)
    {

      arma::vec Dji_h = arma::vectorise(X.row(j) - xrow_i) / h;
      double Kji_h = arma::prod(
        (arma::abs(Dji_h) < 1.0) %
          arma::square(1.0 - arma::square(Dji_h)) * 15.0 / 16.0
      );

      Dhat(i) += Kji_h;
      Dhat(j) += Kji_h;

      for (size_t m = 0; m < n_m; ++m){

        Nhat(i, m) += Y(j, m) * Kji_h;
        Nhat(j, m) += Y(i, m) * Kji_h;
      }
    }
  }

  for (size_t m = 0; m < n_m; ++m){

    for (size_t i = 0; i < n_n; ++i){

      if (Dhat(i) != 0){

        cv += pow(Y(i, m) - Nhat(i, m) / Dhat(i), 2) * p_Y(m);

      }else{

        cv += pow(Y(i, m), 2) * p_Y(m);
      }
    }
  }

  cv /= n_n;

  return cv;
}

// [[Rcpp::export]]
double CVMNW_K2B_w_rcpp(arma::mat X,
                        arma::mat Y,
                        arma::vec h,
                        arma::vec p_Y,
                        arma::vec w){

  arma::uword n_n = X.n_rows;
  arma::uword n_m = Y.n_cols;
  arma::mat Nhat(n_n, n_m);
  arma::vec Dhat(n_n);
  double cv = 0.0;

  for (size_t i = 0; i < n_n; ++i){

    const auto xrow_i = X.row(i);

    for (size_t j = 0; j < i; ++j)
    {

      arma::vec Dji_h = arma::vectorise(X.row(j) - xrow_i) / h;
      double Kji_h = arma::prod(
        (arma::abs(Dji_h) < 1.0) %
          arma::square(1.0 - arma::square(Dji_h)) * 15.0 / 16.0
      );

      Dhat(i) += Kji_h * w(j);
      Dhat(j) += Kji_h * w(i);

      for (size_t m = 0; m < n_m; ++m){

        Nhat(i, m) += Y(j, m) * Kji_h * w(j);
        Nhat(j, m) += Y(i, m) * Kji_h * w(i);
      }
    }
  }

  for (size_t m = 0; m < n_m; ++m){

    for (size_t i = 0; i < n_n; ++i){

      if (Dhat(i) != 0){

        cv += pow(Y(i, m) - Nhat(i, m) / Dhat(i), 2) * w(i) * p_Y(m);

      }else{

        cv += pow(Y(i, m), 2) * w(i) * p_Y(m);
      }
    }
  }

  cv /= n_n;

  return cv;
}

// Using K4_Biweight kernel

// [[Rcpp::export]]
double CVMNW_K4B_rcpp(arma::mat X,
                      arma::mat Y,
                      arma::vec h,
                      arma::vec p_Y){

  arma::uword n_n = X.n_rows;
  arma::uword n_m = Y.n_cols;
  arma::mat Nhat(n_n, n_m);
  arma::vec Dhat(n_n);
  double cv = 0.0;

  for (size_t i = 0; i < n_n; ++i){

    const auto xrow_i = X.row(i);

    for (size_t j = 0; j < i; ++j)
    {

      arma::vec Dji_h_sq = arma::square(
        arma::vectorise(X.row(j) - xrow_i) / h
      );
      double Kji_h = arma::prod(
        (Dji_h_sq < 1) % (1.0 - Dji_h_sq * 3.0) %
          arma::square(1.0 - Dji_h_sq) * 105.0 / 64.0
      );

      Dhat(i) += Kji_h;
      Dhat(j) += Kji_h;

      for (size_t m = 0; m < n_m; ++m){

        Nhat(i, m) += Y(j, m) * Kji_h;
        Nhat(j, m) += Y(i, m) * Kji_h;
      }
    }
  }

  for (size_t m = 0; m < n_m; ++m){

    for (size_t i = 0; i < n_n; ++i){

      if (Dhat(i) != 0){

        cv += pow(Y(i, m) - Nhat(i, m) / Dhat(i), 2) * p_Y(m);

      }else{

        cv += pow(Y(i, m), 2) * p_Y(m);
      }
    }
  }

  cv /= n_n;

  return cv;
}

// [[Rcpp::export]]
double CVMNW_K4B_w_rcpp(arma::mat X,
                        arma::mat Y,
                        arma::vec h,
                        arma::vec p_Y,
                        arma::vec w){

  arma::uword n_n = X.n_rows;
  arma::uword n_m = Y.n_cols;
  arma::mat Nhat(n_n, n_m);
  arma::vec Dhat(n_n);
  double cv = 0.0;

  for (size_t i = 0; i < n_n; ++i){

    const auto xrow_i = X.row(i);

    for (size_t j = 0; j < i; ++j)
    {

      arma::vec Dji_h_sq = arma::square(
        arma::vectorise(X.row(j) - xrow_i) / h
      );
      double Kji_h = arma::prod(
        (Dji_h_sq < 1) % (1.0 - Dji_h_sq * 3.0) %
          arma::square(1.0 - Dji_h_sq) * 105.0 / 64.0
      );

      Dhat(i) += Kji_h * w(j);
      Dhat(j) += Kji_h * w(i);

      for (size_t m = 0; m < n_m; ++m){

        Nhat(i, m) += Y(j, m) * Kji_h * w(j);
        Nhat(j, m) += Y(i, m) * Kji_h * w(i);
      }
    }
  }

  for (size_t m = 0; m < n_m; ++m){

    for (size_t i = 0; i < n_n; ++i){

      if (Dhat(i) != 0){

        cv += pow(Y(i, m) - Nhat(i, m) / Dhat(i), 2) * w(i) * p_Y(m);

      }else{

        cv += pow(Y(i, m), 2) * w(i) * p_Y(m);
      }
    }
  }

  cv /= n_n;

  return cv;
}

// Using Gaussian kernel

// [[Rcpp::export]]
double CVMNW_KG_rcpp(arma::mat X,
                     arma::mat Y,
                     arma::vec h,
                     arma::vec p_Y){

  arma::uword n_n = X.n_rows;
  arma::uword n_m = Y.n_cols;
  arma::mat Nhat(n_n, n_m);
  arma::vec Dhat(n_n);
  double cv = 0.0;

  for (size_t i = 0; i < n_n; ++i){

    const auto xrow_i = X.row(i);

    for (size_t j = 0; j < i; ++j)
    {

      arma::vec Dji_h = arma::vectorise(X.row(j) - xrow_i) / h;
      double Kji_h = arma::prod(
        arma::exp(-arma::square(Dji_h) / 2.0) * 0.39894228
      );

      Dhat(i) += Kji_h;
      Dhat(j) += Kji_h;

      for (size_t m = 0; m < n_m; ++m){

        Nhat(i, m) += Y(j, m) * Kji_h;
        Nhat(j, m) += Y(i, m) * Kji_h;
      }
    }
  }

  for (size_t m = 0; m < n_m; ++m){

    for (size_t i = 0; i < n_n; ++i){

      if (Dhat(i) != 0){

        cv += pow(Y(i, m) - Nhat(i, m) / Dhat(i), 2) * p_Y(m);

      }else{

        cv += pow(Y(i, m), 2) * p_Y(m);
      }
    }
  }

  cv /= n_n;

  return cv;
}

// [[Rcpp::export]]
double CVMNW_KG_w_rcpp(arma::mat X,
                       arma::mat Y,
                       arma::vec h,
                       arma::vec p_Y,
                       arma::vec w){

  arma::uword n_n = X.n_rows;
  arma::uword n_m = Y.n_cols;
  arma::mat Nhat(n_n, n_m);
  arma::vec Dhat(n_n);
  double cv = 0.0;

  for (size_t i = 0; i < n_n; ++i){

    const auto xrow_i = X.row(i);

    for (size_t j = 0; j < i; ++j)
    {

      arma::vec Dji_h = arma::vectorise(X.row(j) - xrow_i) / h;
      double Kji_h = arma::prod(
        arma::exp(-arma::square(Dji_h) / 2.0) * 0.39894228
      );

      Dhat(i) += Kji_h * w(j);
      Dhat(j) += Kji_h * w(i);

      for (size_t m = 0; m < n_m; ++m){

        Nhat(i, m) += Y(j, m) * Kji_h * w(j);
        Nhat(j, m) += Y(i, m) * Kji_h * w(i);
      }
    }
  }

  for (size_t m = 0; m < n_m; ++m){

    for (size_t i = 0; i < n_n; ++i){

      if (Dhat(i) != 0){

        cv += pow(Y(i, m) - Nhat(i, m) / Dhat(i), 2) * w(i) * p_Y(m);

      }else{

        cv += pow(Y(i, m), 2) * w(i) * p_Y(m);
      }
    }
  }

  cv /= n_n;

  return cv;
}

// Trimming to [0, 1]

// Using K4_Biweight kernel

// [[Rcpp::export]]
double CVMNWdist_K4B_rcpp(arma::mat X,
                          arma::mat Y_CP,
                          arma::vec h,
                          arma::vec p_Y){

  arma::uword n_n = X.n_rows;
  arma::uword n_m = Y_CP.n_cols;
  arma::mat Nhat(n_n, n_m);
  arma::vec Dhat(n_n);
  double cv = 0.0;

  for (size_t i = 0; i < n_n; ++i){

    const auto xrow_i = X.row(i);

    for (size_t j = 0; j < i; ++j)
    {

      arma::vec Dji_h_sq = arma::square(
        arma::vectorise(X.row(j) - xrow_i) / h
      );
      double Kji_h = arma::prod(
        (Dji_h_sq < 1) % (1.0 - Dji_h_sq * 3.0) %
          arma::square(1.0 - Dji_h_sq) * 105.0 / 64.0
      );

      Dhat(i) += Kji_h;
      Dhat(j) += Kji_h;

      for (size_t m = 0; m < n_m; ++m){

        Nhat(i, m) += Y_CP(j, m) * Kji_h;
        Nhat(j, m) += Y_CP(i, m) * Kji_h;
      }
    }
  }

  for (size_t m = 0; m < n_m; ++m){

    for (size_t i = 0; i < n_n; ++i){

      if (Dhat(i) != 0){

        double Yhat_i_m = Nhat(i, m) / Dhat(i);

        if (Yhat_i_m > 1){

          cv += pow(Y_CP(i, m) - 1, 2) * p_Y(m);

        }else if (Yhat_i_m < 0){

          cv += pow(Y_CP(i, m), 2) * p_Y(m);

        }else{

          cv += pow(Y_CP(i, m) - Yhat_i_m, 2) * p_Y(m);
        }
      }else{

        cv += pow(Y_CP(i, m), 2) * p_Y(m);
      }
    }
  }

  cv /= n_n;

  return cv;
}

// [[Rcpp::export]]
double CVMNWdist_K4B_w_rcpp(arma::mat X,
                            arma::mat Y_CP,
                            arma::vec h,
                            arma::vec p_Y,
                            arma::vec w){

  arma::uword n_n = X.n_rows;
  arma::uword n_m = Y_CP.n_cols;
  arma::mat Nhat(n_n, n_m);
  arma::vec Dhat(n_n);
  double cv = 0.0;

  for (size_t i = 0; i < n_n; ++i){

    const auto xrow_i = X.row(i);

    for (size_t j = 0; j < i; ++j)
    {

      arma::vec Dji_h_sq = arma::square(
        arma::vectorise(X.row(j) - xrow_i) / h
      );
      double Kji_h = arma::prod(
        (Dji_h_sq < 1) % (1.0 - Dji_h_sq * 3.0) %
          arma::square(1.0 - Dji_h_sq) * 105.0 / 64.0
      );

      Dhat(i) += Kji_h * w(j);
      Dhat(j) += Kji_h * w(i);

      for (size_t m = 0; m < n_m; ++m){

        Nhat(i, m) += Y_CP(j, m) * Kji_h * w(j);
        Nhat(j, m) += Y_CP(i, m) * Kji_h * w(i);
      }
    }
  }

  for (size_t m = 0; m < n_m; ++m){

    for (size_t i = 0; i < n_n; ++i){

      if (Dhat(i) != 0){

        double Yhat_i_m = Nhat(i, m) / Dhat(i);

        if (Yhat_i_m > 1){

          cv += pow(Y_CP(i, m) - 1, 2) * w(i) * p_Y(m);

        }else if (Yhat_i_m < 0){

          cv += pow(Y_CP(i, m), 2) * w(i) * p_Y(m);

        }else{

          cv += pow(Y_CP(i, m) - Yhat_i_m, 2) * w(i) * p_Y(m);
        }
      }else{

        cv += pow(Y_CP(i, m), 2) * w(i) * p_Y(m);
      }
    }
  }

  cv /= n_n;

  return cv;
}

// conditional distribution with univariate response

// Y must be in ascending order

// Using K2_Biweight kernel

// [[Rcpp::export]]
double CVDNWuniY_K2B_rcpp(arma::mat X,
                          arma::vec Y,
                          arma::vec h,
                          arma::uvec rank_y_in_Y,
                          arma::vec p_y){

  arma::uword n_n = X.n_rows;
  arma::uword n_m = rank_y_in_Y.n_elem;
  double cv = 0;

  for (size_t i = 0; i < n_n; ++i){

    arma::vec Ki_h(n_n);
    const auto xrow_i = X.row(i);

    for (size_t j = 0; j < n_n; ++j){

      if (j != i){
        arma::vec Dji_h = arma::vectorise(X.row(j) - xrow_i) / h;
        Ki_h(j) = arma::prod(
          (arma::abs(Dji_h) < 1.0) %
            arma::square(1.0 - arma::square(Dji_h)) * 15.0 / 16.0
        );
      }
    }

    double Dhat_i = arma::sum(Ki_h);
    arma::vec cKi_h = arma::cumsum(Ki_h);
    arma::vec Nhat_i = cKi_h.elem(rank_y_in_Y - 1);

    if (Dhat_i != 0){

      // arma::vec Fhat_i = arma::clamp(Nhat_i / Dhat_i, 0.0, 1.0);
      cv += arma::sum(arma::pow((rank_y_in_Y >= (i + 1)) - Nhat_i / Dhat_i, 2) % p_y);

    }else{

      cv += arma::sum((rank_y_in_Y >= (i + 1)) % p_y);
    }
  }

  cv /= n_n;

  return(cv);
}

// [[Rcpp::export]]
double CVDNWuniY_K2B_w_rcpp(arma::mat X,
                            arma::vec Y,
                            arma::vec h,
                            arma::uvec rank_y_in_Y,
                            arma::vec p_y,
                            arma::vec w){

  arma::uword n_n = X.n_rows;
  arma::uword n_m = rank_y_in_Y.n_elem;
  double cv = 0;

  for (size_t i = 0; i < n_n; ++i){

    arma::vec Ki_h(n_n);
    const auto xrow_i = X.row(i);

    for (size_t j = 0; j < n_n; ++j){

      if (j != i){
        arma::vec Dji_h = arma::vectorise(X.row(j) - xrow_i) / h;
        Ki_h(j) = arma::prod(
          (arma::abs(Dji_h) < 1.0) %
            arma::square(1.0 - arma::square(Dji_h)) * 15.0 / 16.0
        ) * w(j);
      }
    }

    double Dhat_i = arma::sum(Ki_h);
    arma::vec cKi_h = arma::cumsum(Ki_h);
    arma::vec Nhat_i = cKi_h.elem(rank_y_in_Y - 1);

    if (Dhat_i != 0){

      // arma::vec Fhat_i = arma::clamp(Nhat_i / Dhat_i, 0.0, 1.0);
      cv += arma::sum(arma::pow((rank_y_in_Y >= (i + 1)) - Nhat_i / Dhat_i, 2) % p_y) * w(i);

    }else{

      cv += arma::sum((rank_y_in_Y >= (i + 1)) % p_y) * w(i);
    }
  }

  cv /= n_n;

  return(cv);
}

// Using K4_Biweight kernel

// [[Rcpp::export]]
double CVDNWuniY_K4B_rcpp(arma::mat X,
                          arma::vec Y,
                          arma::vec h,
                          arma::uvec rank_y_in_Y,
                          arma::vec p_y){

  arma::uword n_n = X.n_rows;
  arma::uword n_m = rank_y_in_Y.n_elem;
  double cv = 0;

  for (size_t i = 0; i < n_n; ++i){

    arma::vec Ki_h(n_n);
    const auto xrow_i = X.row(i);

    for (size_t j = 0; j < n_n; ++j){

      if (j != i){
        arma::vec Dji_h_sq = arma::square(
          arma::vectorise(X.row(j) - xrow_i) / h
        );
        Ki_h(j) = arma::prod(
          (Dji_h_sq < 1) % (1.0 - Dji_h_sq * 3.0) %
            arma::square(1.0 - Dji_h_sq) * 105.0 / 64.0
        );
      }
    }

    double Dhat_i = arma::sum(Ki_h);
    arma::vec cKi_h = arma::cumsum(Ki_h);
    arma::vec Nhat_i = cKi_h.elem(rank_y_in_Y - 1);

    if (Dhat_i != 0){

      arma::vec Fhat_i = arma::clamp(Nhat_i / Dhat_i, 0.0, 1.0);
      cv += arma::sum(arma::pow((rank_y_in_Y >= (i + 1)) - Fhat_i, 2) % p_y);

    }else{

      cv += arma::sum((rank_y_in_Y >= (i + 1)) % p_y);
    }
  }

  cv /= n_n;

  return(cv);
}

// [[Rcpp::export]]
double CVDNWuniY_K4B_w_rcpp(arma::mat X,
                            arma::vec Y,
                            arma::vec h,
                            arma::uvec rank_y_in_Y,
                            arma::vec p_y,
                            arma::vec w){

  arma::uword n_n = X.n_rows;
  arma::uword n_m = rank_y_in_Y.n_elem;
  double cv = 0;

  for (size_t i = 0; i < n_n; ++i){

    arma::vec Ki_h(n_n);
    const auto xrow_i = X.row(i);

    for (size_t j = 0; j < n_n; ++j){

      if (j != i){
        arma::vec Dji_h_sq = arma::square(
          arma::vectorise(X.row(j) - xrow_i) / h
        );
        Ki_h(j) = arma::prod(
          (Dji_h_sq < 1) % (1.0 - Dji_h_sq * 3.0) %
            arma::square(1.0 - Dji_h_sq) * 105.0 / 64.0
        ) * w(j);
      }
    }

    double Dhat_i = arma::sum(Ki_h);
    arma::vec cKi_h = arma::cumsum(Ki_h);
    arma::vec Nhat_i = cKi_h.elem(rank_y_in_Y - 1);

    if (Dhat_i != 0){

      arma::vec Fhat_i = arma::clamp(Nhat_i / Dhat_i, 0.0, 1.0);
      cv += arma::sum(arma::pow((rank_y_in_Y >= (i + 1)) - Fhat_i, 2) % p_y) * w(i);

    }else{

      cv += arma::sum((rank_y_in_Y >= (i + 1)) % p_y) * w(i);
    }
  }

  cv /= n_n;

  return(cv);
}

// Using Gaussian kernel

// [[Rcpp::export]]
double CVDNWuniY_KG_rcpp(arma::mat X,
                         arma::vec Y,
                         arma::vec h,
                         arma::uvec rank_y_in_Y,
                         arma::vec p_y){

  arma::uword n_n = X.n_rows;
  arma::uword n_m = rank_y_in_Y.n_elem;
  double cv = 0;

  for (size_t i = 0; i < n_n; ++i){

    arma::vec Ki_h(n_n);
    const auto xrow_i = X.row(i);

    for (size_t j = 0; j < n_n; ++j){

      if (j != i){
        arma::vec Dji_h = arma::vectorise(X.row(j) - xrow_i) / h;
        Ki_h(j) = arma::prod(
          arma::exp(-arma::square(Dji_h) / 2.0) * 0.39894228
        );
      }
    }

    double Dhat_i = arma::sum(Ki_h);
    arma::vec cKi_h = arma::cumsum(Ki_h);
    arma::vec Nhat_i = cKi_h.elem(rank_y_in_Y - 1);

    if (Dhat_i != 0){

      cv += arma::sum(arma::pow((rank_y_in_Y >= (i + 1)) - Nhat_i / Dhat_i, 2) % p_y);

    }else{

      cv += arma::sum((rank_y_in_Y >= (i + 1)) % p_y);
    }
  }

  cv /= n_n;

  return(cv);
}

// [[Rcpp::export]]
double CVDNWuniY_KG_w_rcpp(arma::mat X,
                           arma::vec Y,
                           arma::vec h,
                           arma::uvec rank_y_in_Y,
                           arma::vec p_y,
                           arma::vec w){

  arma::uword n_n = X.n_rows;
  arma::uword n_m = rank_y_in_Y.n_elem;
  double cv = 0;

  for (size_t i = 0; i < n_n; ++i){

    arma::vec Ki_h(n_n);
    const auto xrow_i = X.row(i);

    for (size_t j = 0; j < n_n; ++j){

      if (j != i){
        arma::vec Dji_h = arma::vectorise(X.row(j) - xrow_i) / h;
        Ki_h(j) = arma::prod(
          arma::exp(-arma::square(Dji_h) / 2.0) * 0.39894228
        ) * w(j);
      }
    }

    double Dhat_i = arma::sum(Ki_h);
    arma::vec cKi_h = arma::cumsum(Ki_h);
    arma::vec Nhat_i = cKi_h.elem(rank_y_in_Y - 1);

    if (Dhat_i != 0){

      cv += arma::sum(arma::pow((rank_y_in_Y >= (i + 1)) - Nhat_i / Dhat_i, 2) % p_y) * w(i);

    }else{

      cv += arma::sum((rank_y_in_Y >= (i + 1)) % p_y) * w(i);
    }
  }

  cv /= n_n;

  return(cv);
}

// semiparametric sum of squares criterion for Nadaraya-Watson estimator

// conditional mean

// Using K2_Biweight kernel

// [[Rcpp::export]]
double SSMNW_K2B_rcpp(arma::mat X,
                      arma::mat Y,
                      arma::vec h,
                      arma::vec p_Y){

  arma::uword n_n = X.n_rows;
  arma::uword n_m = Y.n_cols;
  arma::mat Nhat(n_n, n_m);
  arma::vec Dhat(n_n);
  double cv = 0.0;

  for (size_t i = 0; i < n_n; ++i){

    const auto xrow_i = X.row(i);

    for (size_t j = 0; j <= i; ++j)
    {

      arma::vec Dji_h = arma::vectorise(X.row(j) - xrow_i) / h;
      double Kji_h = arma::prod(
        (arma::abs(Dji_h) < 1.0) %
          arma::square(1.0 - arma::square(Dji_h)) * 15.0 / 16.0
      );

      Dhat(i) += Kji_h;
      Dhat(j) += Kji_h;

      for (size_t m = 0; m < n_m; ++m){

        Nhat(i, m) += Y(j, m) * Kji_h;
        Nhat(j, m) += Y(i, m) * Kji_h;
      }
    }
  }

  for (size_t m = 0; m < n_m; ++m){

    for (size_t i = 0; i < n_n; ++i){

      if (Dhat(i) != 0){

        cv += pow(Y(i, m) - Nhat(i, m) / Dhat(i), 2) * p_Y(m);

      }else{

        cv += pow(Y(i, m), 2) * p_Y(m);
      }
    }
  }

  cv /= n_n;

  return cv;
}

// [[Rcpp::export]]
double SSMNW_K2B_w_rcpp(arma::mat X,
                        arma::mat Y,
                        arma::vec h,
                        arma::vec p_Y,
                        arma::vec w){

  arma::uword n_n = X.n_rows;
  arma::uword n_m = Y.n_cols;
  arma::mat Nhat(n_n, n_m);
  arma::vec Dhat(n_n);
  double cv = 0.0;

  for (size_t i = 0; i < n_n; ++i){

    const auto xrow_i = X.row(i);

    for (size_t j = 0; j <= i; ++j)
    {

      arma::vec Dji_h = arma::vectorise(X.row(j) - xrow_i) / h;
      double Kji_h = arma::prod(
        (arma::abs(Dji_h) < 1.0) %
          arma::square(1.0 - arma::square(Dji_h)) * 15.0 / 16.0
      );

      Dhat(i) += Kji_h * w(j);
      Dhat(j) += Kji_h * w(i);

      for (size_t m = 0; m < n_m; ++m){

        Nhat(i, m) += Y(j, m) * Kji_h * w(j);
        Nhat(j, m) += Y(i, m) * Kji_h * w(i);
      }
    }
  }

  for (size_t m = 0; m < n_m; ++m){

    for (size_t i = 0; i < n_n; ++i){

      if (Dhat(i) != 0){

        cv += pow(Y(i, m) - Nhat(i, m) / Dhat(i), 2) * w(i) * p_Y(m);

      }else{

        cv += pow(Y(i, m), 2) * w(i) * p_Y(m);
      }
    }
  }

  cv /= n_n;

  return cv;
}

// Using K4_Biweight kernel

// [[Rcpp::export]]
double SSMNW_K4B_rcpp(arma::mat X,
                      arma::mat Y,
                      arma::vec h,
                      arma::vec p_Y){

  arma::uword n_n = X.n_rows;
  arma::uword n_m = Y.n_cols;
  arma::mat Nhat(n_n, n_m);
  arma::vec Dhat(n_n);
  double cv = 0.0;

  for (size_t i = 0; i < n_n; ++i){

    const auto xrow_i = X.row(i);

    for (size_t j = 0; j <= i; ++j)
    {

      arma::vec Dji_h_sq = arma::square(
        arma::vectorise(X.row(j) - xrow_i) / h
      );
      double Kji_h = arma::prod(
        (Dji_h_sq < 1) % (1.0 - Dji_h_sq * 3.0) %
          arma::square(1.0 - Dji_h_sq) * 105.0 / 64.0
      );

      Dhat(i) += Kji_h;
      Dhat(j) += Kji_h;

      for (size_t m = 0; m < n_m; ++m){

        Nhat(i, m) += Y(j, m) * Kji_h;
        Nhat(j, m) += Y(i, m) * Kji_h;
      }
    }
  }

  for (size_t m = 0; m < n_m; ++m){

    for (size_t i = 0; i < n_n; ++i){

      if (Dhat(i) != 0){

        cv += pow(Y(i, m) - Nhat(i, m) / Dhat(i), 2) * p_Y(m);

      }else{

        cv += pow(Y(i, m), 2) * p_Y(m);
      }
    }
  }

  cv /= n_n;

  return cv;
}

// [[Rcpp::export]]
double SSMNW_K4B_w_rcpp(arma::mat X,
                        arma::mat Y,
                        arma::vec h,
                        arma::vec p_Y,
                        arma::vec w){

  arma::uword n_n = X.n_rows;
  arma::uword n_m = Y.n_cols;
  arma::mat Nhat(n_n, n_m);
  arma::vec Dhat(n_n);
  double cv = 0.0;

  for (size_t i = 0; i < n_n; ++i){

    const auto xrow_i = X.row(i);

    for (size_t j = 0; j <= i; ++j)
    {

      arma::vec Dji_h_sq = arma::square(
        arma::vectorise(X.row(j) - xrow_i) / h
      );
      double Kji_h = arma::prod(
        (Dji_h_sq < 1) % (1.0 - Dji_h_sq * 3.0) %
          arma::square(1.0 - Dji_h_sq) * 105.0 / 64.0
      );

      Dhat(i) += Kji_h * w(j);
      Dhat(j) += Kji_h * w(i);

      for (size_t m = 0; m < n_m; ++m){

        Nhat(i, m) += Y(j, m) * Kji_h * w(j);
        Nhat(j, m) += Y(i, m) * Kji_h * w(i);
      }
    }
  }

  for (size_t m = 0; m < n_m; ++m){

    for (size_t i = 0; i < n_n; ++i){

      if (Dhat(i) != 0){

        cv += pow(Y(i, m) - Nhat(i, m) / Dhat(i), 2) * w(i) * p_Y(m);

      }else{

        cv += pow(Y(i, m), 2) * w(i) * p_Y(m);
      }
    }
  }

  cv /= n_n;

  return cv;
}

// Using Gaussian kernel

// [[Rcpp::export]]
double SSMNW_KG_rcpp(arma::mat X,
                     arma::mat Y,
                     arma::vec h,
                     arma::vec p_Y){

  arma::uword n_n = X.n_rows;
  arma::uword n_m = Y.n_cols;
  arma::mat Nhat(n_n, n_m);
  arma::vec Dhat(n_n);
  double cv = 0.0;

  for (size_t i = 0; i < n_n; ++i){

    const auto xrow_i = X.row(i);

    for (size_t j = 0; j <= i; ++j)
    {

      arma::vec Dji_h = arma::vectorise(X.row(j) - xrow_i) / h;
      double Kji_h = arma::prod(
        arma::exp(-arma::square(Dji_h) / 2.0) * 0.39894228
      );

      Dhat(i) += Kji_h;
      Dhat(j) += Kji_h;

      for (size_t m = 0; m < n_m; ++m){

        Nhat(i, m) += Y(j, m) * Kji_h;
        Nhat(j, m) += Y(i, m) * Kji_h;
      }
    }
  }

  for (size_t m = 0; m < n_m; ++m){

    for (size_t i = 0; i < n_n; ++i){

      if (Dhat(i) != 0){

        cv += pow(Y(i, m) - Nhat(i, m) / Dhat(i), 2) * p_Y(m);

      }else{

        cv += pow(Y(i, m), 2) * p_Y(m);
      }
    }
  }

  cv /= n_n;

  return cv;
}

// [[Rcpp::export]]
double SSMNW_KG_w_rcpp(arma::mat X,
                       arma::mat Y,
                       arma::vec h,
                       arma::vec p_Y,
                       arma::vec w){

  arma::uword n_n = X.n_rows;
  arma::uword n_m = Y.n_cols;
  arma::mat Nhat(n_n, n_m);
  arma::vec Dhat(n_n);
  double cv = 0.0;

  for (size_t i = 0; i < n_n; ++i){

    const auto xrow_i = X.row(i);

    for (size_t j = 0; j <= i; ++j)
    {

      arma::vec Dji_h = arma::vectorise(X.row(j) - xrow_i) / h;
      double Kji_h = arma::prod(
        arma::exp(-arma::square(Dji_h) / 2.0) * 0.39894228
      );

      Dhat(i) += Kji_h * w(j);
      Dhat(j) += Kji_h * w(i);

      for (size_t m = 0; m < n_m; ++m){

        Nhat(i, m) += Y(j, m) * Kji_h * w(j);
        Nhat(j, m) += Y(i, m) * Kji_h * w(i);
      }
    }
  }

  for (size_t m = 0; m < n_m; ++m){

    for (size_t i = 0; i < n_n; ++i){

      if (Dhat(i) != 0){

        cv += pow(Y(i, m) - Nhat(i, m) / Dhat(i), 2) * w(i) * p_Y(m);

      }else{

        cv += pow(Y(i, m), 2) * w(i) * p_Y(m);
      }
    }
  }

  cv /= n_n;

  return cv;
}

// Trimming to [0, 1]

// Using K4_Biweight kernel

// [[Rcpp::export]]
double SSMNWdist_K4B_rcpp(arma::mat X,
                          arma::mat Y_CP,
                          arma::vec h,
                          arma::vec p_Y){

  arma::uword n_n = X.n_rows;
  arma::uword n_m = Y_CP.n_cols;
  arma::mat Nhat(n_n, n_m);
  arma::vec Dhat(n_n);
  double cv = 0.0;

  for (size_t i = 0; i < n_n; ++i){

    const auto xrow_i = X.row(i);

    for (size_t j = 0; j <= i; ++j)
    {

      arma::vec Dji_h_sq = arma::square(
        arma::vectorise(X.row(j) - xrow_i) / h
      );
      double Kji_h = arma::prod(
        (Dji_h_sq < 1) % (1.0 - Dji_h_sq * 3.0) %
          arma::square(1.0 - Dji_h_sq) * 105.0 / 64.0
      );

      Dhat(i) += Kji_h;
      Dhat(j) += Kji_h;

      for (size_t m = 0; m < n_m; ++m){

        Nhat(i, m) += Y_CP(j, m) * Kji_h;
        Nhat(j, m) += Y_CP(i, m) * Kji_h;
      }
    }
  }

  for (size_t m = 0; m < n_m; ++m){

    for (size_t i = 0; i < n_n; ++i){

      if (Dhat(i) != 0){

        double Yhat_i_m = Nhat(i, m) / Dhat(i);

        if (Yhat_i_m > 1){

          cv += pow(Y_CP(i, m) - 1, 2) * p_Y(m);

        }else if (Yhat_i_m < 0){

          cv += pow(Y_CP(i, m), 2) * p_Y(m);

        }else{

          cv += pow(Y_CP(i, m) - Yhat_i_m, 2) * p_Y(m);
        }
      }else{

        cv += pow(Y_CP(i, m), 2) * p_Y(m);
      }
    }
  }

  cv /= n_n;

  return cv;
}

// [[Rcpp::export]]
double SSMNWdist_K4B_w_rcpp(arma::mat X,
                            arma::mat Y_CP,
                            arma::vec h,
                            arma::vec p_Y,
                            arma::vec w){

  arma::uword n_n = X.n_rows;
  arma::uword n_m = Y_CP.n_cols;
  arma::mat Nhat(n_n, n_m);
  arma::vec Dhat(n_n);
  double cv = 0.0;

  for (size_t i = 0; i < n_n; ++i){

    const auto xrow_i = X.row(i);

    for (size_t j = 0; j <= i; ++j)
    {

      arma::vec Dji_h_sq = arma::square(
        arma::vectorise(X.row(j) - xrow_i) / h
      );
      double Kji_h = arma::prod(
        (Dji_h_sq < 1) % (1.0 - Dji_h_sq * 3.0) %
          arma::square(1.0 - Dji_h_sq) * 105.0 / 64.0
      );

      Dhat(i) += Kji_h * w(j);
      Dhat(j) += Kji_h * w(i);

      for (size_t m = 0; m < n_m; ++m){

        Nhat(i, m) += Y_CP(j, m) * Kji_h * w(j);
        Nhat(j, m) += Y_CP(i, m) * Kji_h * w(i);
      }
    }
  }

  for (size_t m = 0; m < n_m; ++m){

    for (size_t i = 0; i < n_n; ++i){

      if (Dhat(i) != 0){

        double Yhat_i_m = Nhat(i, m) / Dhat(i);

        if (Yhat_i_m > 1){

          cv += pow(Y_CP(i, m) - 1, 2) * w(i) * p_Y(m);

        }else if (Yhat_i_m < 0){

          cv += pow(Y_CP(i, m), 2) * w(i) * p_Y(m);

        }else{

          cv += pow(Y_CP(i, m) - Yhat_i_m, 2) * w(i) * p_Y(m);
        }
      }else{

        cv += pow(Y_CP(i, m), 2) * w(i) * p_Y(m);
      }
    }
  }

  cv /= n_n;

  return cv;
}

// conditional distribution with univariate response

// Y must be in ascending order

// Using K2_Biweight kernel

// [[Rcpp::export]]
double SSDNWuniY_K2B_rcpp(arma::mat X,
                          arma::vec Y,
                          arma::vec h,
                          arma::uvec rank_y_in_Y,
                          arma::vec p_y){

  arma::uword n_n = X.n_rows;
  arma::uword n_m = rank_y_in_Y.n_elem;
  double cv = 0;

  for (size_t i = 0; i < n_n; ++i){

    arma::vec Ki_h(n_n);
    const auto xrow_i = X.row(i);

    for (size_t j = 0; j < n_n; ++j){

      arma::vec Dji_h = arma::vectorise(X.row(j) - xrow_i) / h;
      Ki_h(j) = arma::prod(
        (arma::abs(Dji_h) < 1.0) %
          arma::square(1.0 - arma::square(Dji_h)) * 15.0 / 16.0
      );
    }

    double Dhat_i = arma::sum(Ki_h);
    arma::vec cKi_h = arma::cumsum(Ki_h);
    arma::vec Nhat_i = cKi_h.elem(rank_y_in_Y - 1);

    if (Dhat_i != 0){

      // arma::vec Fhat_i = arma::clamp(Nhat_i / Dhat_i, 0.0, 1.0);
      cv += arma::sum(arma::pow((rank_y_in_Y >= (i + 1)) - Nhat_i / Dhat_i, 2) % p_y);

    }else{

      cv += arma::sum((rank_y_in_Y >= (i + 1)) % p_y);
    }
  }

  cv /= n_n;

  return(cv);
}

// [[Rcpp::export]]
double SSDNWuniY_K2B_w_rcpp(arma::mat X,
                            arma::vec Y,
                            arma::vec h,
                            arma::uvec rank_y_in_Y,
                            arma::vec p_y,
                            arma::vec w){

  arma::uword n_n = X.n_rows;
  arma::uword n_m = rank_y_in_Y.n_elem;
  double cv = 0;

  for (size_t i = 0; i < n_n; ++i){

    arma::vec Ki_h(n_n);
    const auto xrow_i = X.row(i);

    for (size_t j = 0; j < n_n; ++j){

      arma::vec Dji_h = arma::vectorise(X.row(j) - xrow_i) / h;
      Ki_h(j) = arma::prod(
        (arma::abs(Dji_h) < 1.0) %
          arma::square(1.0 - arma::square(Dji_h)) * 15.0 / 16.0
      ) * w(j);
    }

    double Dhat_i = arma::sum(Ki_h);
    arma::vec cKi_h = arma::cumsum(Ki_h);
    arma::vec Nhat_i = cKi_h.elem(rank_y_in_Y - 1);

    if (Dhat_i != 0){

      // arma::vec Fhat_i = arma::clamp(Nhat_i / Dhat_i, 0.0, 1.0);
      cv += arma::sum(arma::pow((rank_y_in_Y >= (i + 1)) - Nhat_i / Dhat_i, 2) % p_y) * w(i);

    }else{

      cv += arma::sum((rank_y_in_Y >= (i + 1)) % p_y) * w(i);
    }
  }

  cv /= n_n;

  return(cv);
}

// Using K4_Biweight kernel

// [[Rcpp::export]]
double SSDNWuniY_K4B_rcpp(arma::mat X,
                          arma::vec Y,
                          arma::vec h,
                          arma::uvec rank_y_in_Y,
                          arma::vec p_y){

  arma::uword n_n = X.n_rows;
  arma::uword n_m = rank_y_in_Y.n_elem;
  double cv = 0;

  for (size_t i = 0; i < n_n; ++i){

    arma::vec Ki_h(n_n);
    const auto xrow_i = X.row(i);

    for (size_t j = 0; j < n_n; ++j){

      arma::vec Dji_h_sq = arma::square(
        arma::vectorise(X.row(j) - xrow_i) / h
      );
      Ki_h(j) = arma::prod(
        (Dji_h_sq < 1) % (1.0 - Dji_h_sq * 3.0) %
          arma::square(1.0 - Dji_h_sq) * 105.0 / 64.0
      );
    }

    double Dhat_i = arma::sum(Ki_h);
    arma::vec cKi_h = arma::cumsum(Ki_h);
    arma::vec Nhat_i = cKi_h.elem(rank_y_in_Y - 1);

    if (Dhat_i != 0){

      arma::vec Fhat_i = arma::clamp(Nhat_i / Dhat_i, 0.0, 1.0);
      cv += arma::sum(arma::pow((rank_y_in_Y >= (i + 1)) - Fhat_i, 2) % p_y);

    }else{

      cv += arma::sum((rank_y_in_Y >= (i + 1)) % p_y);
    }
  }

  cv /= n_n;

  return(cv);
}

// [[Rcpp::export]]
double SSDNWuniY_K4B_w_rcpp(arma::mat X,
                            arma::vec Y,
                            arma::vec h,
                            arma::uvec rank_y_in_Y,
                            arma::vec p_y,
                            arma::vec w){

  arma::uword n_n = X.n_rows;
  arma::uword n_m = rank_y_in_Y.n_elem;
  double cv = 0;

  for (size_t i = 0; i < n_n; ++i){

    arma::vec Ki_h(n_n);
    const auto xrow_i = X.row(i);

    for (size_t j = 0; j < n_n; ++j){

      arma::vec Dji_h_sq = arma::square(
        arma::vectorise(X.row(j) - xrow_i) / h
      );
      Ki_h(j) = arma::prod(
        (Dji_h_sq < 1) % (1.0 - Dji_h_sq * 3.0) %
          arma::square(1.0 - Dji_h_sq) * 105.0 / 64.0
      ) * w(j);
    }

    double Dhat_i = arma::sum(Ki_h);
    arma::vec cKi_h = arma::cumsum(Ki_h);
    arma::vec Nhat_i = cKi_h.elem(rank_y_in_Y - 1);

    if (Dhat_i != 0){

      arma::vec Fhat_i = arma::clamp(Nhat_i / Dhat_i, 0.0, 1.0);
      cv += arma::sum(arma::pow((rank_y_in_Y >= (i + 1)) - Fhat_i, 2) % p_y) * w(i);

    }else{

      cv += arma::sum((rank_y_in_Y >= (i + 1)) % p_y) * w(i);
    }
  }

  cv /= n_n;

  return(cv);
}

// Using Gaussian kernel

// [[Rcpp::export]]
double SSDNWuniY_KG_rcpp(arma::mat X,
                         arma::vec Y,
                         arma::vec h,
                         arma::uvec rank_y_in_Y,
                         arma::vec p_y){

  arma::uword n_n = X.n_rows;
  arma::uword n_m = rank_y_in_Y.n_elem;
  double cv = 0;

  for (size_t i = 0; i < n_n; ++i){

    arma::vec Ki_h(n_n);
    const auto xrow_i = X.row(i);

    for (size_t j = 0; j < n_n; ++j){

      arma::vec Dji_h = arma::vectorise(X.row(j) - xrow_i) / h;
      Ki_h(j) = arma::prod(
        arma::exp(-arma::square(Dji_h) / 2.0) * 0.39894228
      );
    }

    double Dhat_i = arma::sum(Ki_h);
    arma::vec cKi_h = arma::cumsum(Ki_h);
    arma::vec Nhat_i = cKi_h.elem(rank_y_in_Y - 1);

    if (Dhat_i != 0){

      cv += arma::sum(arma::pow((rank_y_in_Y >= (i + 1)) - Nhat_i / Dhat_i, 2) % p_y);

    }else{

      cv += arma::sum((rank_y_in_Y >= (i + 1)) % p_y);
    }
  }

  cv /= n_n;

  return(cv);
}

// [[Rcpp::export]]
double SSDNWuniY_KG_w_rcpp(arma::mat X,
                           arma::vec Y,
                           arma::vec h,
                           arma::uvec rank_y_in_Y,
                           arma::vec p_y,
                           arma::vec w){

  arma::uword n_n = X.n_rows;
  arma::uword n_m = rank_y_in_Y.n_elem;
  double cv = 0;

  for (size_t i = 0; i < n_n; ++i){

    arma::vec Ki_h(n_n);
    const auto xrow_i = X.row(i);

    for (size_t j = 0; j < n_n; ++j){

      arma::vec Dji_h = arma::vectorise(X.row(j) - xrow_i) / h;
      Ki_h(j) = arma::prod(
        arma::exp(-arma::square(Dji_h) / 2.0) * 0.39894228
      ) * w(j);
    }

    double Dhat_i = arma::sum(Ki_h);
    arma::vec cKi_h = arma::cumsum(Ki_h);
    arma::vec Nhat_i = cKi_h.elem(rank_y_in_Y - 1);

    if (Dhat_i != 0){

      cv += arma::sum(arma::pow((rank_y_in_Y >= (i + 1)) - Nhat_i / Dhat_i, 2) % p_y) * w(i);

    }else{

      cv += arma::sum((rank_y_in_Y >= (i + 1)) % p_y) * w(i);
    }
  }

  cv /= n_n;

  return(cv);
}
