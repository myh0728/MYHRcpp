#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::vec KDE_rcpp(arma::mat X,
                   arma::mat x,
                   Rcpp::Function K,
                   arma::vec h){
  arma::uword n_n = X.n_rows;
  arma::uword n_k = x.n_rows;
  arma::vec Dhat(n_k);
  for (arma::uword k = 0; k < n_k; ++k){
    for (arma::uword i = 0; i < n_n; ++i){
      arma::vec Dik_h = vectorise(X.row(i)-x.row(k))/h;
      double Kik_h = prod(Rcpp::as<arma::vec>(K(Dik_h)));
      //prod((abs(Dik_h) < 1)%pow(1-pow(Dik_h, 2), 2)*15/16);
      Dhat(k) += Kik_h;
    }
  }
  Dhat = Dhat/n_n;
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
arma::vec KDE_K2B_rcpp_chatgpt(const arma::mat& X,
                               const arma::mat& x,
                               const arma::vec& h) {
  const arma::uword n_n = X.n_rows;
  const arma::uword n_k = x.n_rows;
  arma::vec Dhat(n_k, arma::fill::zeros);

  // Precompute constants for the kernel
  const double c1 = 15.0 / 16.0;
  const double c2 = 1.0 / 5.0;

  // Loop over the points x
  for (arma::uword k = 0; k < n_k; ++k) {
    // Loop over the points X
    for (arma::uword i = 0; i < n_n; ++i) {
      // Compute the kernel value
      const arma::vec Dik = (X.row(i) - x.row(k)).t() / h;
      const arma::vec Kik_h = c1 * pow((1.0 - pow(Dik, 2)), 2) % (abs(Dik) < 1.0);

      // Compute the contribution of this point to the estimate
      Dhat(k) += arma::prod(Kik_h);
    }
  }

  // Normalize the estimate
  Dhat /= n_n;

  return Dhat;
}
