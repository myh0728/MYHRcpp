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

// [[Rcpp::export]]
arma::vec KDE_K2B_rcpp_chatgpt(const arma::mat& X,
                               const arma::mat& x,
                               const arma::vec& h){

  const arma::uword n_n = X.n_rows;
  const arma::uword n_k = x.n_rows;
  arma::vec Dhat(n_k, arma::fill::zeros);

  // Precompute constants for the kernel
  const double c1 = 15.0 / 16.0;

  // Loop over the points x
  for (arma::uword k = 0; k < n_k; ++k) {

    // Loop over the points X
    for (arma::uword i = 0; i < n_n; ++i) {

      // Compute the kernel value
      const arma::vec Dik = (X.row(i) - x.row(k)).t() / h;
      const arma::vec Kik_h = c1 * arma::square((1.0 - arma::square(Dik))) %
        (arma::abs(Dik) < 1) / h;

      // Compute the contribution of this point to the estimate
      Dhat(k) += arma::prod(Kik_h);
    }
  }

  // Normalize the estimate
  Dhat /= n_n;

  return Dhat;
}

// [[Rcpp::export]]
arma::vec KDEcv_K2B_rcpp_o1(arma::mat X,
                            arma::vec h){

  arma::uword n_n = X.n_rows;
  arma::vec Dhat(n_n);

  for (size_t i = 0; i < n_n; ++i){

    for (size_t j = 0; j < n_n; ++j){

      if (j != i){

        arma::vec Dji_h = arma::vectorise(X.row(j) - X.row(i)) / h;
        double Kji_h = arma::prod(
          (arma::abs(Dji_h) < 1) %
            arma::square(1.0 - arma::square(Dji_h)) * 15.0 / 16.0 / h
        );
        Dhat(i) += Kji_h;
      }
    }
  }

  Dhat /= n_n-1;
  return(Dhat);
}

// [[Rcpp::export]]
arma::vec KDEcv_K2B_w_rcpp_o1(arma::mat X,
                              arma::vec h,
                              arma::vec w){
  arma::uword n_n = X.n_rows;
  arma::vec Dhat(n_n);

  for (size_t i = 0; i < n_n; ++i){

    for (size_t j = 0; j < n_n; ++j){

      if (j != i){

        arma::vec Dji_h = arma::vectorise(X.row(j) - X.row(i)) / h;
        double Kji_h = arma::prod(
          (arma::abs(Dji_h) < 1) %
            arma::square(1.0 - arma::square(Dji_h)) * 15.0 / 16.0 / h
        ) * w(j);
        Dhat(i) += Kji_h;
      }
    }
  }

  Dhat /= n_n-1;
  return(Dhat);
}

// [[Rcpp::export]]
arma::vec KDE_K4B_rcpp_o1(arma::mat X,
                          arma::mat x,
                          arma::vec h){

  arma::uword n_n = X.n_rows;
  arma::uword n_k = x.n_rows;
  arma::vec Dhat(n_k);

  for (size_t k = 0; k < n_k; ++k){

    for (size_t i = 0; i < n_n; ++i){

      arma::vec Dik_h = arma::vectorise(X.row(i) - x.row(k)) / h;
      double Kik_h = arma::prod(
        (arma::abs(Dik_h) < 1) % (1.0 - 3.0 * arma::square(Dik_h)) %
          arma::square(1.0 - arma::square(Dik_h)) * 105.0 / 64.0 / h
      );
      Dhat(k) += Kik_h;
    }
  }

  Dhat = Dhat/n_n;
  Dhat(find(Dhat < 0)).fill(0);
  return(Dhat);
}

// [[Rcpp::export]]
arma::vec KDEcv_K4B_rcpp_o1(arma::mat X,
                            arma::vec h){

  arma::uword n_n = X.n_rows;
  arma::vec Dhat(n_n);

  for (size_t i = 0; i < n_n; ++i){

    for (size_t j = 0; j < n_n; ++j){

      if (j != i){

        arma::vec Dji_h = arma::vectorise(X.row(j) - X.row(i)) / h;
        double Kji_h = arma::prod(
          (arma::abs(Dji_h) < 1) % (1.0 - 3.0 * arma::square(Dji_h)) %
            arma::square(1.0 - arma::square(Dji_h)) * 105.0 / 64.0 / h
        );
        Dhat(i) += Kji_h;
      }
    }
  }

  Dhat /= n_n-1;
  Dhat(find(Dhat < 0)).fill(0);
  return(Dhat);
}







// Nadaraya-Watson estimation for conditional distribution

// [[Rcpp::export]]
arma::mat NWF_K2B_rcpp(arma::mat X,
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
  for (arma::uword k = 0; k < n_k; ++k){
    for (arma::uword i = 0; i < n_n; ++i){
      arma::vec Dik_h = vectorise(X.row(i)-x.row(k))/h;
      double Kik_h = prod((abs(Dik_h) < 1)%pow(1-pow(Dik_h, 2), 2)*15/16);
      Dhat(k) += Kik_h;
      for (arma::uword m = 0; m < n_m; ++m){
        Nhat(k, m) += prod(Y.row(i) <= y.row(m))*Kik_h;
      }
    }
    if (Dhat(k) != 0){
      Fhat.row(k) = Nhat.row(k)/Dhat(k);
    }
  }
  return(Fhat);
}

// [[Rcpp::export]]
arma::mat NWcv_K2B_rcpp_o1(arma::mat X,
                           arma::mat Y,
                           arma::vec h){
  arma::uword n_n = X.n_rows;
  arma::uword n_m = Y.n_cols;
  arma::mat Nhat(n_n, n_m);
  arma::vec Dhat(n_n);
  arma::mat Yhat(n_n, n_m);
  for (arma::uword i = 0; i < n_n; ++i){
    for (arma::uword j = 0; j < n_n; ++j){
      if (j != i){
        arma::vec Dji_h = vectorise(X.row(j)-X.row(i))/h;
        double Kji_h = prod((abs(Dji_h) < 1)%pow(1-pow(Dji_h, 2), 2)*15/16);
        Dhat(i) += Kji_h;
        for (arma::uword m = 0; m < n_m; ++m){
          Nhat(i, m) += Y(j, m)*Kji_h;
        }
      }
    }
    if (Dhat(i) != 0){
      Yhat.row(i) = Nhat.row(i)/Dhat(i);
    }
  }
  return(Yhat);
}

