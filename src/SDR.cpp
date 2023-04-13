#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
double DRCV_K2B_rcpp(arma::mat X,
                     arma::mat Y_CP,
                     arma::vec h){

  arma::uword n_n = X.n_rows;
  arma::mat Nhat(n_n, n_n);
  arma::vec Dhat(n_n);
  double cv = 0.0;;

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

      for (size_t m = 0; m < n_n; ++m){

        Nhat(i, m) += Y_CP(j, m) * Kji_h;
        Nhat(j, m) += Y_CP(i, m) * Kji_h;
      }
    }
  }

  for (size_t m = 0; m < n_n; ++m){

    for (size_t i = 0; i < n_n; ++i){

      if (Dhat(i) != 0){

        cv += std::pow(Y_CP(i, m) - Nhat(i, m) / Dhat(i), 2);
      }
    }
  }

  return cv;
}



