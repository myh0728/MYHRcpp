#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// kernel density estimation

// Using Gaussian kernel

// [[Rcpp::export]]
arma::vec KDE_KG_rcpp(arma::mat X,
                      arma::mat x,
                      arma::vec h){

  arma::uword n_n = X.n_rows;
  arma::uword n_k = x.n_rows;
  arma::vec Dhat(n_k);

  for (size_t k = 0; k < n_k; ++k){

    for (size_t i = 0; i < n_n; ++i){

      arma::vec Dik_h = arma::vectorise(X.row(i) - x.row(k)) / h;
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

    for (size_t i = 0; i < n_n; ++i){

      arma::vec Dik_h = arma::vectorise(X.row(i) - x.row(k)) / h;
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

  for (size_t i = 0; i < n_n; ++i)
  {

    const auto xrow_i = X.row(i);

    for (size_t j = 0; j < i; ++j)
    {

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

  for (size_t i = 0; i < n_n; ++i)
  {

    const auto xrow_i = X.row(i);

    for (size_t j = 0; j < i; ++j)
    {

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

// Using K2_Biweight kernel

// [[Rcpp::export]]
arma::vec KDE_K2B_rcpp(arma::mat X,
                       arma::mat x,
                       arma::vec h){

  arma::uword n_n = X.n_rows;
  arma::uword n_k = x.n_rows;
  arma::vec Dhat(n_k);

  for (size_t k = 0; k < n_k; ++k){

    for (size_t i = 0; i < n_n; ++i){

      arma::vec Dik_h = arma::vectorise(X.row(i) - x.row(k)) / h;
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

    for (size_t i = 0; i < n_n; ++i){

      arma::vec Dik_h = arma::vectorise(X.row(i) - x.row(k)) / h;
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

  for (size_t i = 0; i < n_n; ++i)
  {

    const auto xrow_i = X.row(i);

    for (size_t j = 0; j < i; ++j)
    {

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

  for (size_t i = 0; i < n_n; ++i)
  {

    const auto xrow_i = X.row(i);

    for (size_t j = 0; j < i; ++j)
    {

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

    for (size_t i = 0; i < n_n; ++i){

      arma::vec Dik_h_sq = arma::square(
        arma::vectorise(X.row(i) - x.row(k)) / h
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

    for (size_t i = 0; i < n_n; ++i){

      arma::vec Dik_h_sq = arma::square(
        arma::vectorise(X.row(i) - x.row(k)) / h
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

  for (size_t i = 0; i < n_n; ++i)
  {

    const auto xrow_i = X.row(i);

    for (size_t j = 0; j < i; ++j)
    {

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

  for (size_t i = 0; i < n_n; ++i)
  {

    const auto xrow_i = X.row(i);

    for (size_t j = 0; j < i; ++j)
    {

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










// Nadaraya-Watson estimation

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
  for (arma::uword k = 0; k < n_k; ++k){
    for (arma::uword i = 0; i < n_n; ++i){
      arma::vec Dik_h = vectorise(X.row(i)-x.row(k))/h;
      double Kik_h = prod((abs(Dik_h) < 1)%pow(1-pow(Dik_h, 2), 2)*15/16);
      Dhat(k) += Kik_h;
      for (arma::uword m = 0; m < n_m; ++m){
        Nhat(k, m) += Y(i, m)*Kik_h;
      }
    }
    if (Dhat(k) != 0){
      yhat.row(k) = Nhat.row(k)/Dhat(k);
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
  for (arma::uword k = 0; k < n_k; ++k){
    for (arma::uword i = 0; i < n_n; ++i){
      arma::vec Dik_h = vectorise(X.row(i)-x.row(k))/h;
      double Kik_h = prod((abs(Dik_h) < 1)%pow(1-pow(Dik_h, 2), 2)*15/16)*w(i);
      Dhat(k) += Kik_h;
      for (arma::uword m = 0; m < n_m; ++m){
        Nhat(k, m) += Y(i, m)*Kik_h;
      }
    }
    if (Dhat(k) != 0){
      yhat.row(k) = Nhat.row(k)/Dhat(k);
    }
  }
  return(yhat);
}

// [[Rcpp::export]]
arma::mat NWcv_K2B_rcpp(
    arma::mat X,
    arma::mat Y,
    arma::vec h){

  arma::uword n_n = X.n_rows;
  arma::uword n_m = Y.n_cols;
  arma::mat Nhat(n_n, n_m);
  arma::vec Dhat(n_n);
  arma::mat Yhat(n_n, n_m);

  for (size_t i = 0; i < n_n; ++i)
  {

    const auto xrow_i = X.row(i);

    for (size_t j = 0; j < i; ++j)
    {

      arma::vec Dji_h = arma::vectorise(X.row(j) - xrow_i) / h;
      double Kji_h = prod(
        (arma::abs(Dji_h) < 1.0) %
          arma::pow(1.0 - arma::pow(Dji_h, 2.0), 2.0) * 15.0 / 16.0
      );

      Dhat(i) += Kji_h;
      Dhat(j) += Kji_h;

      Nhat.row(i) += Y.row(j) * Kji_h;
      Nhat.row(j) += Y.row(i) * Kji_h;

    }

  }

  for (size_t i = 0; i < n_n; ++i)
  {

    if (Dhat(i) != 0){

      Yhat.row(i) = Nhat.row(i) / Dhat(i);

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
  for (arma::uword i = 0; i < n_n; ++i){
    for (arma::uword j = 0; j < n_n; ++j){
      if (j != i){
        arma::vec Dji_h = vectorise(X.row(j)-X.row(i))/h;
        double Kji_h = prod((abs(Dji_h) < 1)%pow(1-pow(Dji_h, 2), 2)*15/16)*w(j);
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
  for (arma::uword k = 0; k < n_k; ++k){
    for (arma::uword i = 0; i < n_n; ++i){
      arma::vec Dik_h = vectorise(X.row(i)-x.row(k))/h;
      double Kik_h = prod((abs(Dik_h) < 1)%(1-3*pow(Dik_h, 2))%
                          pow(1-pow(Dik_h, 2), 2)*105/64);
      Dhat(k) += Kik_h;
      for (arma::uword m = 0; m < n_m; ++m){
        Nhat(k, m) += Y(i, m)*Kik_h;
      }
    }
    if (Dhat(k) != 0){
      yhat.row(k) = Nhat.row(k)/Dhat(k);
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
  for (arma::uword k = 0; k < n_k; ++k){
    for (arma::uword i = 0; i < n_n; ++i){
      arma::vec Dik_h = vectorise(X.row(i)-x.row(k))/h;
      double Kik_h = prod((abs(Dik_h) < 1)%(1-3*pow(Dik_h, 2))%
                          pow(1-pow(Dik_h, 2), 2)*105/64)*w(i);
      Dhat(k) += Kik_h;
      for (arma::uword m = 0; m < n_m; ++m){
        Nhat(k, m) += Y(i, m)*Kik_h;
      }
    }
    if (Dhat(k) != 0){
      yhat.row(k) = Nhat.row(k)/Dhat(k);
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
  for (arma::uword i = 0; i < n_n; ++i){
    for (arma::uword j = 0; j < n_n; ++j){
      if (j != i){
        arma::vec Dji_h = vectorise(X.row(j)-X.row(i))/h;
        double Kji_h = prod((abs(Dji_h) < 1)%(1-3*pow(Dji_h, 2))%
                            pow(1-pow(Dji_h, 2), 2)*105/64);
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
  for (arma::uword i = 0; i < n_n; ++i){
    for (arma::uword j = 0; j < n_n; ++j){
      if (j != i){
        arma::vec Dji_h = vectorise(X.row(j)-X.row(i))/h;
        double Kji_h = prod((abs(Dji_h) < 1)%(1-3*pow(Dji_h, 2))%
                            pow(1-pow(Dji_h, 2), 2)*105/64)*w(j);
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



