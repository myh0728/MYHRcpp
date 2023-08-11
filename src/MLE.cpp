#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
List diff_lL_normal(arma::mat X,
                    arma::vec Y,
                    double alpha,
                    arma::vec beta,
                    double sigma){

  arma::uword n_n = X.n_rows;
  arma::uword n_p = X.n_cols;
  arma::vec gradient(n_p + 2);
  arma::mat hessian(n_p + 2, n_p + 2);
  List result;

  arma::vec X_intercept = arma::vec({1});

  for (size_t i = 0; i < n_n; ++i){

    const arma::vec Xrow_i = X.row(i).t();
    const arma::vec extXrow_i = arma::join_vert(X_intercept, Xrow_i);
    const double res_i = Y(i) - alpha - arma::dot(Xrow_i, beta);
    const double res_i_square = pow(res_i, 2);
    const arma::vec cross = extXrow_i * res_i;
    gradient.subvec(0, n_p) += cross;
    gradient(n_p + 1) += res_i_square;
    hessian.submat(0, 0, n_p, n_p) += -extXrow_i * extXrow_i.t();
    hessian(n_p + 1, n_p + 1) += -res_i_square;
    hessian.submat(0, n_p + 1, n_p, n_p + 1) += -cross;
    hessian.submat(n_p + 1, 0, n_p + 1, n_p) += -cross.t();
  }

  gradient.subvec(0, n_p) = gradient.subvec(0, n_p) / pow(sigma, 2) / n_n;
  gradient(n_p + 1) = gradient(n_p + 1) / pow(sigma, 3) / n_n - pow(sigma, -1);
  hessian.submat(0, 0, n_p, n_p) = hessian.submat(0, 0, n_p, n_p) / pow(sigma, 2) / n_n;
  hessian(n_p + 1, n_p + 1) = hessian(n_p + 1, n_p + 1) * 3 / pow(sigma, 4) / n_n + pow(sigma, -2);
  hessian.submat(0, n_p + 1, n_p, n_p + 1) = hessian.submat(0, n_p + 1, n_p, n_p + 1) * 2 / pow(sigma, 3) / n_n;
  hessian.submat(n_p + 1, 0, n_p + 1, n_p) = hessian.submat(n_p + 1, 0, n_p + 1, n_p) * 2 / pow(sigma, 3) / n_n;

  result["gradient"] = gradient;
  result["hessian"] = hessian;
  return result;
}


