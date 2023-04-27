#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::vec dG1_normal_rcpp(arma::mat X,
                          arma::vec Y,
                          arma::uword n1,
                          double alpha,
                          arma::vec beta,
                          double sigma,
                          arma::uword iter_max,
                          double stop_tol){

  arma::uword n_n = X.n_rows;

  arma::vec SI_all = alpha + X * beta;
  arma::mat f_yx_2_all(n_n - n1, n_n);

  for (size_t i = 0; i < n_n; ++i){

    for (size_t i2 = n1; i2 < n_n; ++i2){

      f_yx_2_all(i2 - n1, i) = exp(-pow(Y(i2) - SI_all(i), 2) /
        (2.0 * pow(sigma, 2))) * 0.3989423;
    }
  }

  arma::vec dG1 = arma::ones(n_n) / n_n;

  for (size_t iter = 1; iter <= iter_max; ++iter){

    arma::vec dF1_2 = f_yx_2_all * dG1;
    arma::vec dG1_new = 1.0 /
      (arma::sum(f_yx_2_all.each_col() / dF1_2, 0).t() + n1);
    double sum_dG1_new = arma::sum(dG1_new);

    if (sum_dG1_new == 0){

      break;

    }else{

      dG1_new /= sum_dG1_new;
    }

    if (arma::sum(arma::abs(dG1_new - dG1)) > stop_tol){

      dG1 = dG1_new;
    }else{

      break;
    }
  }

  return dG1;
}

// [[Rcpp::export]]
double lpL_normal_rcpp(arma::mat X,
                       arma::vec Y,
                       arma::uword n1,
                       double alpha,
                       arma::vec beta,
                       double sigma,
                       arma::uword iter_max,
                       double stop_tol){

  double lL = 0;
  arma::uword n_n = X.n_rows;

  arma::vec SI_all = alpha + X * beta;
  arma::mat f_yx_2_all(n_n - n1, n_n);
  arma::vec f_yx_all(n_n);

  for (size_t i = 0; i < n_n; ++i){

    for (size_t i2 = n1; i2 < n_n; ++i2){

      f_yx_2_all(i2 - n1, i) = exp(-pow(Y(i2) - SI_all(i), 2) /
        (2.0 * pow(sigma, 2))) * 0.3989423;

      if (i2 == i){

        f_yx_all(i) = f_yx_2_all(i - n1, i);
      }
    }
  }

  for (size_t i1 = 0; i1 < n1; ++i1){

    f_yx_all(i1) = exp(-pow(Y(i1) - SI_all(i1), 2) /
      (2.0 * pow(sigma, 2))) * 0.3989423;
  }

  arma::vec dG1 = arma::ones(n_n) / n_n;

  for (size_t iter = 1; iter <= iter_max; ++iter){

    arma::vec dF1_2 = f_yx_2_all * dG1;
    arma::vec dG1_new = 1.0 /
      (arma::sum(f_yx_2_all.each_col() / dF1_2, 0).t() + n1);
    double sum_dG1_new = arma::sum(dG1_new);

    if (sum_dG1_new == 0){

      break;

    }else{

      dG1_new /= sum_dG1_new;
    }

    if (arma::sum(arma::abs(dG1_new - dG1)) > stop_tol){

      dG1 = dG1_new;
    }else{

      break;
    }
  }

  arma::vec dF1_2 = f_yx_2_all * dG1;
  lL = arma::accu(arma::log(f_yx_all)) + arma::accu(arma::log(dG1)) -
    arma::accu(arma::log(dF1_2));

  return lL;
}

// [[Rcpp::export]]
arma::vec dG1alt_normal_rcpp(arma::mat X,
                             arma::vec Y,
                             arma::uword n1,
                             double alpha,
                             arma::vec beta,
                             double sigma,
                             arma::vec gamma,
                             arma::uword iter_max,
                             double stop_tol){

  arma::uword n_n = X.n_rows;

  arma::vec SI_all = alpha + X * beta;
  arma::vec extLI = X * gamma;
  arma::mat f_yx_2_all(n_n - n1, n_n);

  for (size_t i = 0; i < n_n; ++i){

    for (size_t i2 = n1; i2 < n_n; ++i2){

      f_yx_2_all(i2 - n1, i) = exp(-pow(Y(i2) - SI_all(i), 2) /
        (2.0 * pow(sigma, 2))) * 0.3989423 * exp(extLI(i));
    }
  }

  arma::vec dG1 = arma::ones(n_n) / n_n;

  for (size_t iter = 1; iter <= iter_max; ++iter){

    arma::vec dF1_2 = f_yx_2_all * dG1;
    arma::vec dG1_new = 1.0 /
      (arma::sum(f_yx_2_all.each_col() / dF1_2, 0).t() + n1);
    double sum_dG1_new = arma::sum(dG1_new);

    if (sum_dG1_new == 0){

      break;

    }else{

      dG1_new /= sum_dG1_new;
    }

    if (arma::sum(arma::abs(dG1_new - dG1)) > stop_tol){

      dG1 = dG1_new;
    }else{

      break;
    }
  }

  return dG1;
}

// [[Rcpp::export]]
double lpLalt_normal_rcpp(arma::mat X,
                          arma::vec Y,
                          arma::uword n1,
                          double alpha,
                          arma::vec beta,
                          double sigma,
                          arma::vec gamma,
                          arma::uword iter_max,
                          double stop_tol){

  double lL = 0;
  arma::uword n_n = X.n_rows;

  arma::vec SI_all = alpha + X * beta;
  arma::vec extLI = X * gamma;
  arma::mat f_yx_2_all(n_n - n1, n_n);
  arma::vec f_yx_all(n_n);

  for (size_t i = 0; i < n_n; ++i){

    for (size_t i2 = n1; i2 < n_n; ++i2){

      f_yx_2_all(i2 - n1, i) = exp(-pow(Y(i2) - SI_all(i), 2) /
        (2.0 * pow(sigma, 2))) * 0.3989423 * exp(extLI(i));

      if (i2 == i){

        f_yx_all(i) = f_yx_2_all(i - n1, i);
      }
    }
  }

  for (size_t i1 = 0; i1 < n1; ++i1){

    f_yx_all(i1) = exp(-pow(Y(i1) - SI_all(i1), 2) /
      (2.0 * pow(sigma, 2))) * 0.3989423 * exp(extLI(i1));
  }

  arma::vec dG1 = arma::ones(n_n) / n_n;

  for (size_t iter = 1; iter <= iter_max; ++iter){

    arma::vec dF1_2 = f_yx_2_all * dG1;
    arma::vec dG1_new = 1.0 /
      (arma::sum(f_yx_2_all.each_col() / dF1_2, 0).t() + n1);
    double sum_dG1_new = arma::sum(dG1_new);

    if (sum_dG1_new == 0){

      break;

    }else{

      dG1_new /= sum_dG1_new;
    }

    if (arma::sum(arma::abs(dG1_new - dG1)) > stop_tol){

      dG1 = dG1_new;
    }else{

      break;
    }
  }

  arma::vec dF1_2 = f_yx_2_all * dG1;
  lL = arma::accu(arma::log(f_yx_all)) + arma::accu(arma::log(dG1)) -
    arma::accu(arma::log(dF1_2));

  return lL;
}

// [[Rcpp::export]]
arma::vec dG1_logistic_rcpp(arma::mat X,
                            arma::vec Y,
                            arma::uword n1,
                            double alpha,
                            arma::vec beta,
                            arma::uword iter_max,
                            double stop_tol){

  arma::uword n_n = X.n_rows;

  arma::vec logiSI_all = 1.0 / (1.0 + 1.0 / arma::exp(alpha + X * beta));
  arma::mat f_yx_2_all(n_n - n1, n_n);

  for (size_t i = 0; i < n_n; ++i){

    for (size_t i2 = n1; i2 < n_n; ++i2){

      if (Y(i2) == 1){

        f_yx_2_all(i2 - n1, i) = logiSI_all(i);

      }else{

        f_yx_2_all(i2 - n1, i) = 1 - logiSI_all(i);
      }
    }
  }

  arma::vec dG1 = arma::ones(n_n) / n_n;

  for (size_t iter = 1; iter <= iter_max; ++iter){

    arma::vec dF1_2 = f_yx_2_all * dG1;
    arma::vec dG1_new = 1.0 /
      (arma::sum(f_yx_2_all.each_col() / dF1_2, 0).t() + n1);
    double sum_dG1_new = arma::sum(dG1_new);

    if (sum_dG1_new == 0){

      break;

    }else{

      dG1_new /= sum_dG1_new;
    }

    if (arma::sum(arma::abs(dG1_new - dG1)) > stop_tol){

      dG1 = dG1_new;
    }else{

      break;
    }
  }

  return dG1;
}

// [[Rcpp::export]]
double lpL_logistic_rcpp(arma::mat X,
                       arma::vec Y,
                       arma::uword n1,
                       double alpha,
                       arma::vec beta,
                       arma::uword iter_max,
                       double stop_tol){

  double lL = 0;
  arma::uword n_n = X.n_rows;

  arma::vec logiSI_all = 1.0 / (1.0 + 1.0 / arma::exp(alpha + X * beta));
  arma::mat f_yx_2_all(n_n - n1, n_n);
  arma::vec f_yx_all(n_n);

  for (size_t i = 0; i < n_n; ++i){

    for (size_t i2 = n1; i2 < n_n; ++i2){

      if (Y(i2) == 1){

        f_yx_2_all(i2 - n1, i) = logiSI_all(i);

      }else{

        f_yx_2_all(i2 - n1, i) = 1 - logiSI_all(i);
      }

      if (i2 == i){

        f_yx_all(i) = f_yx_2_all(i - n1, i);
      }
    }
  }

  for (size_t i1 = 0; i1 < n1; ++i1){

    if (Y(i1) == 1){

      f_yx_all(i1) = logiSI_all(i1);

    }else{

      f_yx_all(i1) = 1 - logiSI_all(i1);
    }
  }

  arma::vec dG1 = arma::ones(n_n) / n_n;

  for (size_t iter = 1; iter <= iter_max; ++iter){

    arma::vec dF1_2 = f_yx_2_all * dG1;
    arma::vec dG1_new = 1.0 /
      (arma::sum(f_yx_2_all.each_col() / dF1_2, 0).t() + n1);
    double sum_dG1_new = arma::sum(dG1_new);

    if (sum_dG1_new == 0){

      break;

    }else{

      dG1_new /= sum_dG1_new;
    }

    if (arma::sum(arma::abs(dG1_new - dG1)) > stop_tol){

      dG1 = dG1_new;
    }else{

      break;
    }
  }

  arma::vec dF1_2 = f_yx_2_all * dG1;
  lL = arma::accu(arma::log(f_yx_all)) + arma::accu(arma::log(dG1)) -
    arma::accu(arma::log(dF1_2));

  return lL;
}

// [[Rcpp::export]]
arma::vec dG1alt_logistic_rcpp(arma::mat X,
                               arma::vec Y,
                               arma::uword n1,
                               double alpha,
                               arma::vec beta,
                               arma::vec gamma,
                               arma::uword iter_max,
                               double stop_tol){

  arma::uword n_n = X.n_rows;

  arma::vec logiSI_all = 1.0 / (1.0 + 1.0 / arma::exp(alpha + X * beta));
  arma::vec extLI = X * gamma;
  arma::mat f_yx_2_all(n_n - n1, n_n);

  for (size_t i = 0; i < n_n; ++i){

    for (size_t i2 = n1; i2 < n_n; ++i2){

      if (Y(i2) == 1){

        f_yx_2_all(i2 - n1, i) = logiSI_all(i) * exp(extLI(i));

      }else{

        f_yx_2_all(i2 - n1, i) = (1 - logiSI_all(i)) * exp(extLI(i));
      }
    }
  }

  arma::vec dG1 = arma::ones(n_n) / n_n;

  for (size_t iter = 1; iter <= iter_max; ++iter){

    arma::vec dF1_2 = f_yx_2_all * dG1;
    arma::vec dG1_new = 1.0 /
      (arma::sum(f_yx_2_all.each_col() / dF1_2, 0).t() + n1);
    double sum_dG1_new = arma::sum(dG1_new);

    if (sum_dG1_new == 0){

      break;

    }else{

      dG1_new /= sum_dG1_new;
    }

    if (arma::sum(arma::abs(dG1_new - dG1)) > stop_tol){

      dG1 = dG1_new;
    }else{

      break;
    }
  }

  return dG1;
}

// [[Rcpp::export]]
double lpLalt_logistic_rcpp(arma::mat X,
                            arma::vec Y,
                            arma::uword n1,
                            double alpha,
                            arma::vec beta,
                            arma::vec gamma,
                            arma::uword iter_max,
                            double stop_tol){

  double lL = 0;
  arma::uword n_n = X.n_rows;

  arma::vec logiSI_all = 1.0 / (1.0 + 1.0 / arma::exp(alpha + X * beta));
  arma::vec extLI = X * gamma;
  arma::mat f_yx_2_all(n_n - n1, n_n);
  arma::vec f_yx_all(n_n);

  for (size_t i = 0; i < n_n; ++i){

    for (size_t i2 = n1; i2 < n_n; ++i2){

      if (Y(i2) == 1){

        f_yx_2_all(i2 - n1, i) = logiSI_all(i) * exp(extLI(i));

      }else{

        f_yx_2_all(i2 - n1, i) = (1 - logiSI_all(i)) * exp(extLI(i));
      }

      if (i2 == i){

        f_yx_all(i) = f_yx_2_all(i - n1, i);
      }
    }
  }

  for (size_t i1 = 0; i1 < n1; ++i1){

    if (Y(i1) == 1){

      f_yx_all(i1) = logiSI_all(i1) * exp(extLI(i1));

    }else{

      f_yx_all(i1) = (1 - logiSI_all(i1)) * exp(extLI(i1));
    }
  }

  arma::vec dG1 = arma::ones(n_n) / n_n;

  for (size_t iter = 1; iter <= iter_max; ++iter){

    arma::vec dF1_2 = f_yx_2_all * dG1;
    arma::vec dG1_new = 1.0 /
      (arma::sum(f_yx_2_all.each_col() / dF1_2, 0).t() + n1);
    double sum_dG1_new = arma::sum(dG1_new);

    if (sum_dG1_new == 0){

      break;

    }else{

      dG1_new /= sum_dG1_new;
    }

    if (arma::sum(arma::abs(dG1_new - dG1)) > stop_tol){

      dG1 = dG1_new;
    }else{

      break;
    }
  }

  arma::vec dF1_2 = f_yx_2_all * dG1;
  lL = arma::accu(arma::log(f_yx_all)) + arma::accu(arma::log(dG1)) -
    arma::accu(arma::log(dF1_2));

  return lL;
}

