#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat eXsq_rcpp(arma::mat data_X){

  arma::uword n_p = data_X.n_cols;
  arma::mat eXsq(n_p, n_p);

  for (size_t i = 0; i < n_p; ++i){

    for (size_t j = 0; j <= i; ++j){

      eXsq(i, j) = eXsq(j, i) = mean(data_X.col(i) % data_X.col(j));
    }
  }

  return eXsq;
}

// [[Rcpp::export]]
arma::mat eXsq_w_rcpp(arma::mat data_X,
                      arma::vec weight){

  arma::uword n_p = data_X.n_cols;
  arma::mat eXsq(n_p, n_p);

  for (size_t i = 0; i < n_p; ++i){

    for (size_t j = 0; j <= i; ++j){

      eXsq(i, j) = eXsq(j, i) = mean(data_X.col(i) % data_X.col(j) % weight);
    }
  }

  return eXsq;
}








// [[Rcpp::export]]
arma::mat Xsq_lowtri_rcpp(arma::mat data_X){
  arma::uword n_n = data_X.n_rows;
  arma::uword n_p = data_X.n_cols;
  arma::mat Xsq_lowtri(n_n, (n_p+1)*n_p/2);
  arma::uword k = 0;
  for (arma::uword i = 0; i < n_p; ++i){
    for (arma::uword j = i; j < n_p; ++j){
      Xsq_lowtri.col(k) = data_X.col(i) % data_X.col(j);
      k = k+1;
    }
  }
  return Xsq_lowtri;
}

// [[Rcpp::export]]
arma::mat twoXYsym_lowtri_rcpp(arma::mat data_X,
                               arma::mat data_Y){
  arma::uword n_n = data_X.n_rows;
  arma::uword n_p = data_X.n_cols;
  arma::mat twoXYsym_lowtri(n_n, (n_p+1)*n_p/2);
  arma::uword k = 0;
  for (arma::uword i = 0; i < n_p; ++i){
    for (arma::uword j = i; j < n_p; ++j){
      if (i == j){
        twoXYsym_lowtri.col(k) = data_X.col(i) % data_Y.col(j)*2;
      }else{
        twoXYsym_lowtri.col(k) = data_X.col(i) % data_Y.col(j)\
        +data_X.col(j) % data_Y.col(i);
      }
      k = k+1;
    }
  }
  return twoXYsym_lowtri;
}






