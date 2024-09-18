#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
List aux_Lagrange_EY_normal_rcpp(arma::mat X,
                                 double alpha,
                                 arma::vec beta,
                                 double sigma,
                                 double phi,
                                 double eta){

  arma::uword n_n = X.n_rows;
  double value = 0.0;
  double gradient = 0.0;
  double hessian = 0.0;
  List result;

  arma::vec X_intercept = arma::vec({1});

  for (size_t i = 0; i < n_n; ++i){

    const arma::vec Xrow_i = X.row(i).t();
    const arma::vec extXrow_i = arma::join_vert(X_intercept, Xrow_i);
    const double SI_i = alpha + arma::dot(Xrow_i, beta);
    const double Psi_i = SI_i - phi;

    const double denominator = 1.0 + eta * Psi_i;
    if (denominator > 0.0){

      value += log(denominator);
      gradient += Psi_i / denominator;
      hessian += -pow(Psi_i, 2) / pow(denominator, 2);

    } else {

      value += -10000.0;
    }
  }

  result["value"] = value;
  result["gradient"] = gradient;
  result["hessian"] = hessian;
  return result;
}

// [[Rcpp::export]]
List aux_Lagrange_EXsubY_normal_rcpp(arma::mat X,
                                     double alpha,
                                     arma::vec beta,
                                     double sigma,
                                     arma::mat phi,
                                     arma::mat y_pts,
                                     arma::vec eta){

  arma::uword n_n = X.n_rows;
  arma::uword n_p = X.n_cols;
  arma::uword n_k = phi.n_rows;
  arma::uword n_m = n_p * n_k;
  double value = 0.0;
  arma::vec gradient(n_m);
  arma::mat hessian(n_m, n_m);
  List result;

  arma::vec X_intercept = arma::vec({1});

  for (size_t i = 0; i < n_n; ++i){

    const arma::vec Xrow_i = X.row(i).t();
    const arma::vec extXrow_i = arma::join_vert(X_intercept, Xrow_i);
    const double SI_i = alpha + arma::dot(Xrow_i, beta);
    const arma::mat q_pts = (y_pts - SI_i) / sigma;
    const arma::mat cdf_pts = arma::normcdf(q_pts);
    const arma::vec cdf_dist = cdf_pts.col(1) - cdf_pts.col(0);
    const arma::mat dist_e_f = (arma::repmat(Xrow_i.t(), n_k, 1) - phi);
    const arma::vec Psi_i = arma::reshape(
      (dist_e_f % arma::repmat(cdf_dist, 1, n_p)).t(), n_m, 1);

    const double denominator = 1.0 + arma::dot(eta, Psi_i);
    if (denominator > 0.0){

      value += log(denominator);
      gradient += Psi_i / denominator;
      hessian += -Psi_i * Psi_i.t() / pow(denominator, 2);

    } else {

      value += -10000.0;
    }
  }

  result["value"] = value;
  result["gradient"] = gradient;
  result["hessian"] = hessian;
  return result;
}

// [[Rcpp::export]]
List aux_Lagrange_EYsubX_normal_rcpp(arma::mat X,
                                     double alpha,
                                     arma::vec beta,
                                     double sigma,
                                     arma::vec phi,
                                     arma::umat index,
                                     arma::vec eta){

  arma::uword n_n = X.n_rows;
  arma::uword n_m = phi.n_elem;
  double value = 0.0;
  arma::vec gradient(n_m);
  arma::mat hessian(n_m, n_m);
  List result;

  arma::vec X_intercept = arma::vec({1});

  for (size_t i = 0; i < n_n; ++i){

    const arma::vec Xrow_i = X.row(i).t();
    const arma::vec extXrow_i = arma::join_vert(X_intercept, Xrow_i);
    const double SI_i = alpha + arma::dot(Xrow_i, beta);
    const arma::vec c_f = SI_i - phi;
    arma::vec Psi_i(n_m);

    for (size_t m = 0; m < n_m; ++m){

      if (index(i, m) == 1){

        Psi_i(m) += c_f(m);
      }
    }

    const double denominator = 1.0 + arma::dot(eta, Psi_i);
    if (denominator > 0.0){

      value += log(denominator);
      gradient += Psi_i / denominator;
      hessian += -Psi_i * Psi_i.t() / pow(denominator, 2);

    } else {

      value += -10000.0;
    }
  }

  result["value"] = value;
  result["gradient"] = gradient;
  result["hessian"] = hessian;
  return result;
}

// [[Rcpp::export]]
List aux_EY_normal_rcpp(arma::mat X,
                        double alpha,
                        arma::vec beta,
                        double sigma,
                        double phi){

  arma::uword n_n = X.n_rows;
  arma::uword n_p = X.n_cols;
  arma::vec Psi(1);
  arma::mat Psi_square(1, 1);
  arma::mat Psi_gradient(1, n_p + 3);
  List result;

  arma::vec X_intercept = arma::vec({1});

  for (size_t i = 0; i < n_n; ++i){

    const arma::vec Xrow_i = X.row(i).t();
    const arma::vec extXrow_i = arma::join_vert(X_intercept, Xrow_i);
    const double SI_i = alpha + arma::dot(Xrow_i, beta);
    const double c_f = SI_i - phi;

    Psi += c_f;
    Psi_square += Psi * Psi.t();
    Psi_gradient.submat(0, 0, 0, n_p) += extXrow_i.t();
    Psi_gradient(0, n_p + 2) += -1;
  }

  Psi /= n_n;
  Psi_square /= n_n;
  Psi_gradient /= n_n;

  result["score"] = Psi;
  result["score_square"] = Psi_square;
  result["score_gradient"] = Psi_gradient;
  return result;
}

// [[Rcpp::export]]
List aux_EXsubY_normal_rcpp(arma::mat X,
                            double alpha,
                            arma::vec beta,
                            double sigma,
                            arma::mat phi,
                            arma::mat y_pts){

  arma::uword n_n = X.n_rows;
  arma::uword n_p = X.n_cols;
  arma::uword n_k = phi.n_rows;
  arma::uword n_m = n_p * n_k;
  arma::vec Psi(n_m);
  arma::mat Psi_square(n_m, n_m);
  arma::mat Psi_gradient(n_m, n_p + n_m + 2);
  List result;

  arma::vec X_intercept = arma::vec({1});

  for (size_t i = 0; i < n_n; ++i){

    const arma::vec Xrow_i = X.row(i).t();
    const arma::vec extXrow_i = arma::join_vert(X_intercept, Xrow_i);
    const double SI_i = alpha + arma::dot(Xrow_i, beta);
    const arma::mat q_pts = (y_pts - SI_i) / sigma;
    const arma::mat cdf_pts = arma::normcdf(q_pts);
    const arma::mat pdf_pts = arma::normpdf(q_pts);
    const arma::vec cdf_dist = cdf_pts.col(1) - cdf_pts.col(0);
    const arma::vec pdf_dist = pdf_pts.col(1) - pdf_pts.col(0);
    const arma::mat dist = (arma::repmat(Xrow_i.t(), n_k, 1) - phi);
    const arma::mat Psi_i = dist % arma::repmat(cdf_dist, 1, n_p);
    Psi += arma::reshape(Psi_i.t(), n_m, 1);
    Psi_square += arma::reshape(Psi_i.t(), n_m, 1) *
      arma::reshape(Psi_i.t(), 1, n_m);
    Psi_gradient.cols(0, n_p) += arma::reshape((dist %
      arma::repmat(-pdf_dist / sigma, 1, n_p)).t(), n_m, 1) *
      extXrow_i.t();
    Psi_gradient.col(n_p + 1) += arma::reshape((dist %
      arma::repmat(-(q_pts.col(1) / sigma) % pdf_pts.col(1) +
      (q_pts.col(0) / sigma) % pdf_pts.col(0), 1, n_p)).t(), n_m, 1);
    Psi_gradient.cols(n_p + 2, n_p + n_m + 1).diag() +=
      -arma::reshape(arma::repmat(cdf_dist, 1, n_p).t(), n_m, 1);
  }

  Psi /= n_n;
  Psi_square /= n_n;
  Psi_gradient /= n_n;

  result["score"] = Psi;
  result["score_square"] = Psi_square;
  result["score_gradient"] = Psi_gradient;
  return result;
}

// [[Rcpp::export]]
List aux_EYsubX_normal_rcpp(arma::mat X,
                            double alpha,
                            arma::vec beta,
                            double sigma,
                            arma::vec phi,
                            arma::umat index){

  arma::uword n_n = X.n_rows;
  arma::uword n_p = X.n_cols;
  arma::uword n_m = phi.n_elem;
  arma::vec Psi(n_m);
  arma::mat Psi_square(n_m, n_m);
  arma::mat Psi_gradient(n_m, n_p + n_m + 2);
  List result;

  arma::vec X_intercept = arma::vec({1});

  for (size_t i = 0; i < n_n; ++i){

    const arma::vec Xrow_i = X.row(i).t();
    const arma::vec extXrow_i = arma::join_vert(X_intercept, Xrow_i);
    const double SI_i = alpha + arma::dot(Xrow_i, beta);
    arma::vec Psi_i(n_m);

    for (size_t m = 0; m < n_m; ++m){

      if (index(i, m) == 1){

        Psi_i(m) += SI_i - phi(m);
        Psi_gradient.submat(m, 0, m, n_p) += extXrow_i.t();
        Psi_gradient(m, n_p + 2 + m) += -1;
      }
    }

    Psi += Psi_i;
    Psi_square += Psi * Psi.t();
  }

  Psi /= n_n;
  Psi_square /= n_n;
  Psi_gradient /= n_n;

  result["score"] = Psi;
  result["score_square"] = Psi_square;
  result["score_gradient"] = Psi_gradient;
  return result;
}



// [[Rcpp::export]]
List auxLS_Lagrange_EX_normal_rcpp(arma::mat X,
                                   double alpha,
                                   arma::vec beta,
                                   double sigma,
                                   arma::vec phi,
                                   double LS_beta,
                                   arma::vec eta){

  arma::uword n_n = X.n_rows;
  arma::uword n_p = X.n_cols;
  double value = 0.0;
  arma::vec gradient(n_p);
  arma::mat hessian(n_p, n_p);
  List result;

  arma::vec X_intercept = arma::vec({1});

  for (size_t i = 0; i < n_n; ++i){

    const arma::vec Xrow_i = X.row(i).t();
    const arma::vec extXrow_i = arma::join_vert(X_intercept, Xrow_i);
    const double SI_i = alpha + arma::dot(Xrow_i, beta);
    const double e_f = exp((SI_i * 2 + pow(sigma, 2) * LS_beta) * LS_beta / 2.0);
    const arma::vec Psi_i = (Xrow_i - phi) * e_f;

    const double denominator = 1.0 + arma::dot(eta, Psi_i);
    if (denominator > 0.0){

      value += log(denominator);
      gradient += Psi_i / denominator;
      hessian += -Psi_i * Psi_i.t() / pow(denominator, 2);

    } else {

      value += -10000.0;
    }
  }

  result["value"] = value;
  result["gradient"] = gradient;
  result["hessian"] = hessian;
  return result;
}

// [[Rcpp::export]]
List auxLS_Lagrange_EY_normal_rcpp(arma::mat X,
                                   double alpha,
                                   arma::vec beta,
                                   double sigma,
                                   double phi,
                                   double LS_beta,
                                   double eta){

  arma::uword n_n = X.n_rows;
  double value = 0.0;
  double gradient = 0.0;
  double hessian = 0.0;
  List result;

  arma::vec X_intercept = arma::vec({1});

  for (size_t i = 0; i < n_n; ++i){

    const arma::vec Xrow_i = X.row(i).t();
    const arma::vec extXrow_i = arma::join_vert(X_intercept, Xrow_i);
    const double SI_i = alpha + arma::dot(Xrow_i, beta);
    const double c_f = SI_i + pow(sigma, 2) * LS_beta - phi;
    const double e_f = exp((SI_i * 2 + pow(sigma, 2) * LS_beta) * LS_beta / 2.0);
    const double Psi_i = c_f * e_f;

    const double denominator = 1.0 + eta * Psi_i;
    if (denominator > 0.0){

      value += log(denominator);
      gradient += Psi_i / denominator;
      hessian += -pow(Psi_i, 2) / pow(denominator, 2);

    } else {

      value += -10000.0;
    }
  }

  result["value"] = value;
  result["gradient"] = gradient;
  result["hessian"] = hessian;
  return result;
}

// [[Rcpp::export]]
List auxLS_Lagrange_EXsubY_normal_rcpp(arma::mat X,
                                       double alpha,
                                       arma::vec beta,
                                       double sigma,
                                       arma::mat phi,
                                       double LS_beta,
                                       arma::mat y_pts,
                                       arma::vec eta){

  arma::uword n_n = X.n_rows;
  arma::uword n_p = X.n_cols;
  arma::uword n_k = phi.n_rows;
  arma::uword n_m = n_p * n_k;
  double value = 0.0;
  arma::vec gradient(n_m);
  arma::mat hessian(n_m, n_m);
  List result;

  arma::vec X_intercept = arma::vec({1});

  for (size_t i = 0; i < n_n; ++i){

    const arma::vec Xrow_i = X.row(i).t();
    const arma::vec extXrow_i = arma::join_vert(X_intercept, Xrow_i);
    const double SI_i = alpha + arma::dot(Xrow_i, beta);
    const double e_f = exp((SI_i * 2 + pow(sigma, 2) * LS_beta) * LS_beta / 2.0);
    const arma::mat q_pts = (y_pts - SI_i - pow(sigma, 2) * LS_beta) / sigma;
    const arma::mat cdf_pts = arma::normcdf(q_pts);
    const arma::vec cdf_dist = cdf_pts.col(1) - cdf_pts.col(0);
    const arma::mat dist_e_f = (arma::repmat(Xrow_i.t(), n_k, 1) - phi) * e_f;
    const arma::vec Psi_i = arma::reshape(
      (dist_e_f % arma::repmat(cdf_dist, 1, n_p)).t(), n_m, 1);

    const double denominator = 1.0 + arma::dot(eta, Psi_i);
    if (denominator > 0.0){

      value += log(denominator);
      gradient += Psi_i / denominator;
      hessian += -Psi_i * Psi_i.t() / pow(denominator, 2);

    } else {

      value += -10000.0;
    }
  }

  result["value"] = value;
  result["gradient"] = gradient;
  result["hessian"] = hessian;
  return result;
}

// [[Rcpp::export]]
List auxLS_Lagrange_EYsubX_normal_rcpp(arma::mat X,
                                       double alpha,
                                       arma::vec beta,
                                       double sigma,
                                       arma::vec phi,
                                       double LS_beta,
                                       arma::umat index,
                                       arma::vec eta){

  arma::uword n_n = X.n_rows;
  arma::uword n_m = phi.n_elem;
  double value = 0.0;
  arma::vec gradient(n_m);
  arma::mat hessian(n_m, n_m);
  List result;

  arma::vec X_intercept = arma::vec({1});

  for (size_t i = 0; i < n_n; ++i){

    const arma::vec Xrow_i = X.row(i).t();
    const arma::vec extXrow_i = arma::join_vert(X_intercept, Xrow_i);
    const double SI_i = alpha + arma::dot(Xrow_i, beta);
    const arma::vec c_f = SI_i + pow(sigma, 2) * LS_beta - phi;
    const double e_f = exp((SI_i * 2 + pow(sigma, 2) * LS_beta) * LS_beta / 2.0);
    arma::vec Psi_i(n_m);

    for (size_t m = 0; m < n_m; ++m){

      if (index(i, m) == 1){

        Psi_i(m) += c_f(m) * e_f;
      }
    }

    const double denominator = 1.0 + arma::dot(eta, Psi_i);
    if (denominator > 0.0){

      value += log(denominator);
      gradient += Psi_i / denominator;
      hessian += -Psi_i * Psi_i.t() / pow(denominator, 2);

    } else {

      value += -10000.0;
    }
  }

  result["value"] = value;
  result["gradient"] = gradient;
  result["hessian"] = hessian;
  return result;
}

// [[Rcpp::export]]
List auxLS_EX_normal_rcpp(arma::mat X,
                          double alpha,
                          arma::vec beta,
                          double sigma,
                          arma::vec phi,
                          double LS_beta){

  arma::uword n_n = X.n_rows;
  arma::uword n_p = X.n_cols;
  arma::vec Psi(n_p);
  arma::mat Psi_square(n_p, n_p);
  arma::mat Psi_gradient(n_p, n_p * 2 + 3);
  List result;

  arma::vec X_intercept = arma::vec({1});

  for (size_t i = 0; i < n_n; ++i){

    const arma::vec Xrow_i = X.row(i).t();
    const arma::vec extXrow_i = arma::join_vert(X_intercept, Xrow_i);
    const double SI_i = alpha + arma::dot(Xrow_i, beta);
    const double e_f = exp((SI_i * 2 + pow(sigma, 2) * LS_beta) * LS_beta / 2.0);
    const arma::vec Psi_i = (Xrow_i - phi) * e_f;
    Psi += Psi_i;
    Psi_square += Psi_i * Psi_i.t();
    Psi_gradient.cols(0, n_p) += Psi_i * extXrow_i.t() * LS_beta;
    Psi_gradient.col(n_p + 1) += Psi_i * pow(LS_beta, 2) * sigma;
    Psi_gradient.cols(n_p + 2, n_p * 2 + 1).diag() += -e_f;
    Psi_gradient.col(n_p * 2 + 2) += Psi_i * pow(sigma, 2) * LS_beta;
  }

  Psi /= n_n;
  Psi_square /= n_n;
  Psi_gradient /= n_n;

  result["score"] = Psi;
  result["score_square"] = Psi_square;
  result["score_gradient"] = Psi_gradient;
  return result;
}

// [[Rcpp::export]]
List auxLS_EY_normal_rcpp(arma::mat X,
                          double alpha,
                          arma::vec beta,
                          double sigma,
                          double phi,
                          double LS_beta){

  arma::uword n_n = X.n_rows;
  arma::uword n_p = X.n_cols;
  arma::vec Psi(1);
  arma::mat Psi_square(1, 1);
  arma::mat Psi_gradient(1, n_p + 4);
  List result;

  arma::vec X_intercept = arma::vec({1});

  for (size_t i = 0; i < n_n; ++i){

    const arma::vec Xrow_i = X.row(i).t();
    const arma::vec extXrow_i = arma::join_vert(X_intercept, Xrow_i);
    const double SI_i = alpha + arma::dot(Xrow_i, beta);
    const double c_f = SI_i + pow(sigma, 2) * LS_beta - phi;
    const double e_f = exp((SI_i * 2 + pow(sigma, 2) * LS_beta) * LS_beta / 2.0);

    Psi += c_f * e_f;
    Psi_square += Psi * Psi.t();
    Psi_gradient.submat(0, 0, 0, n_p) += extXrow_i.t() * (c_f * LS_beta + 1) * e_f;
    Psi_gradient(0, n_p + 1) += (c_f * LS_beta + 2) * sigma * LS_beta * e_f;
    Psi_gradient(0, n_p + 2) += -e_f;
    Psi_gradient(0, n_p + 3) += (c_f * LS_beta + 1) * pow(sigma, 2) * e_f;
  }

  Psi /= n_n;
  Psi_square /= n_n;
  Psi_gradient /= n_n;

  result["score"] = Psi;
  result["score_square"] = Psi_square;
  result["score_gradient"] = Psi_gradient;
  return result;
}

// [[Rcpp::export]]
List auxLS_EXsubY_normal_rcpp(arma::mat X,
                              double alpha,
                              arma::vec beta,
                              double sigma,
                              arma::mat phi,
                              double LS_beta,
                              arma::mat y_pts){

  arma::uword n_n = X.n_rows;
  arma::uword n_p = X.n_cols;
  arma::uword n_k = phi.n_rows;
  arma::uword n_m = n_p * n_k;
  arma::vec Psi(n_m);
  arma::mat Psi_square(n_m, n_m);
  arma::mat Psi_gradient(n_m, n_p + n_m + 3);
  List result;

  arma::vec X_intercept = arma::vec({1});

  for (size_t i = 0; i < n_n; ++i){

    const arma::vec Xrow_i = X.row(i).t();
    const arma::vec extXrow_i = arma::join_vert(X_intercept, Xrow_i);
    const double SI_i = alpha + arma::dot(Xrow_i, beta);
    const double e_f = exp((SI_i * 2 + pow(sigma, 2) * LS_beta) * LS_beta / 2.0);
    const arma::mat q_pts = (y_pts - SI_i - pow(sigma, 2) * LS_beta) / sigma;
    const arma::mat cdf_pts = arma::normcdf(q_pts);
    const arma::mat pdf_pts = arma::normpdf(q_pts);
    const arma::vec cdf_dist = cdf_pts.col(1) - cdf_pts.col(0);
    const arma::vec pdf_dist = pdf_pts.col(1) - pdf_pts.col(0);
    const arma::mat dist_e_f = (arma::repmat(Xrow_i.t(), n_k, 1) - phi) * e_f;
    const arma::mat Psi_i = dist_e_f % arma::repmat(cdf_dist, 1, n_p);
    Psi += arma::reshape(Psi_i.t(), n_m, 1);
    Psi_square += arma::reshape(Psi_i.t(), n_m, 1) *
      arma::reshape(Psi_i.t(), 1, n_m);
    Psi_gradient.cols(0, n_p) += arma::reshape((dist_e_f %
      arma::repmat(cdf_dist * LS_beta - pdf_dist / sigma, 1, n_p)).t(), n_m, 1) *
      extXrow_i.t();
    Psi_gradient.col(n_p + 1) += arma::reshape((dist_e_f %
      arma::repmat(cdf_dist * sigma * pow(LS_beta, 2) -
      (q_pts.col(1) / sigma + LS_beta * 2) % pdf_pts.col(1) +
      (q_pts.col(0) / sigma + LS_beta * 2) % pdf_pts.col(0), 1, n_p)).t(), n_m, 1);
    Psi_gradient.cols(n_p + 2, n_p + n_m + 1).diag() +=
      arma::reshape(arma::repmat(cdf_dist, 1, n_p).t(), n_m, 1) * (-e_f);
    Psi_gradient.col(n_p + n_m + 2) +=
      arma::reshape((dist_e_f % arma::repmat(cdf_dist * LS_beta * pow(sigma, 2) -
      pdf_dist * sigma, 1, n_p)).t(), n_m, 1);
  }

  Psi /= n_n;
  Psi_square /= n_n;
  Psi_gradient /= n_n;

  result["score"] = Psi;
  result["score_square"] = Psi_square;
  result["score_gradient"] = Psi_gradient;
  return result;
}

// [[Rcpp::export]]
List auxLS_EYsubX_normal_rcpp(arma::mat X,
                              double alpha,
                              arma::vec beta,
                              double sigma,
                              arma::vec phi,
                              double LS_beta,
                              arma::umat index){

  arma::uword n_n = X.n_rows;
  arma::uword n_p = X.n_cols;
  arma::uword n_m = phi.n_elem;
  arma::vec Psi(n_m);
  arma::mat Psi_square(n_m, n_m);
  arma::mat Psi_gradient(n_m, n_p + n_m + 3);
  List result;

  arma::vec X_intercept = arma::vec({1});

  for (size_t i = 0; i < n_n; ++i){

    const arma::vec Xrow_i = X.row(i).t();
    const arma::vec extXrow_i = arma::join_vert(X_intercept, Xrow_i);
    const double SI_i = alpha + arma::dot(Xrow_i, beta);
    const arma::vec c_f = SI_i + pow(sigma, 2) * LS_beta - phi;
    const double e_f = exp((SI_i * 2 + pow(sigma, 2) * LS_beta) * LS_beta / 2.0);
    arma::vec Psi_i(n_m);

    for (size_t m = 0; m < n_m; ++m){

      if (index(i, m) == 1){

        Psi_i(m) += c_f(m) * e_f;
        Psi_gradient.submat(m, 0, m, n_p) += extXrow_i.t() * (c_f(m) * LS_beta + 1) * e_f;
        Psi_gradient(m, n_p + 1) += (c_f(m) * LS_beta + 2) * sigma * LS_beta * e_f;
        Psi_gradient(m, n_p + 2 + m) += -e_f;
        Psi_gradient(m, n_p + n_m + 2) += (c_f(m) * LS_beta + 1) * pow(sigma, 2) * e_f;
      }
    }

    Psi += Psi_i;
    Psi_square += Psi * Psi.t();
  }

  Psi /= n_n;
  Psi_square /= n_n;
  Psi_gradient /= n_n;

  result["score"] = Psi;
  result["score_square"] = Psi_square;
  result["score_gradient"] = Psi_gradient;
  return result;
}



// [[Rcpp::export]]
List auxCS_Lagrange_EY_normal_rcpp(arma::mat X,
                                   double alpha,
                                   arma::vec beta,
                                   double sigma,
                                   double phi,
                                   arma::vec CS_beta,
                                   double eta){

  arma::uword n_n = X.n_rows;
  double value = 0.0;
  double gradient = 0.0;
  double hessian = 0.0;
  List result;

  arma::vec X_intercept = arma::vec({1});

  for (size_t i = 0; i < n_n; ++i){

    const arma::vec Xrow_i = X.row(i).t();
    const arma::vec extXrow_i = arma::join_vert(X_intercept, Xrow_i);
    const double SI_i = alpha + arma::dot(Xrow_i, beta);
    const double eSI_CS_i = exp(arma::dot(Xrow_i, CS_beta));
    const double Psi_i = eSI_CS_i * (SI_i - phi);

    const double denominator = 1.0 + eta * Psi_i;
    if (denominator > 0.0){

      value += log(denominator);
      gradient += Psi_i / denominator;
      hessian += -pow(Psi_i, 2) / pow(denominator, 2);

    } else {

      value += -10000.0;
    }
  }

  result["value"] = value;
  result["gradient"] = gradient;
  result["hessian"] = hessian;
  return result;
}

// [[Rcpp::export]]
List auxCS_Lagrange_EXsubY_normal_rcpp(arma::mat X,
                              double alpha,
                              arma::vec beta,
                              double sigma,
                              arma::mat phi,
                              arma::vec CS_beta,
                              arma::mat y_pts,
                              arma::vec eta){

  arma::uword n_n = X.n_rows;
  arma::uword n_p = X.n_cols;
  arma::uword n_k = phi.n_rows;
  arma::uword n_m = n_p * n_k;
  double value = 0.0;
  arma::vec gradient(n_m);
  arma::mat hessian(n_m, n_m);
  List result;

  arma::vec X_intercept = arma::vec({1});

  for (size_t i = 0; i < n_n; ++i){

    const arma::vec Xrow_i = X.row(i).t();
    const arma::vec extXrow_i = arma::join_vert(X_intercept, Xrow_i);
    const double SI_i = alpha + arma::dot(Xrow_i, beta);
    const double eSI_CS_i = exp(arma::dot(Xrow_i, CS_beta));
    const arma::mat q_pts = (y_pts - SI_i) / sigma;
    const arma::mat cdf_pts = arma::normcdf(q_pts);
    const arma::mat pdf_pts = arma::normpdf(q_pts);
    const arma::vec cdf_dist = cdf_pts.col(1) - cdf_pts.col(0);
    const arma::vec pdf_dist = pdf_pts.col(1) - pdf_pts.col(0);
    const arma::mat x_phi_dist = (arma::repmat(Xrow_i.t(), n_k, 1) - phi);
    const arma::vec Psi_i = arma::reshape(
      x_phi_dist % arma::repmat(cdf_dist, 1, n_p) * eSI_CS_i, n_m, 1);

    const double denominator = 1.0 + arma::dot(eta, Psi_i);
    if (denominator > 0.0){

      value += log(denominator);
      gradient += Psi_i / denominator;
      hessian += -Psi_i * Psi_i.t() / pow(denominator, 2);

    } else {

      value += -10000.0;
    }
  }

  result["value"] = value;
  result["gradient"] = gradient;
  result["hessian"] = hessian;
  return result;
}

// [[Rcpp::export]]
List auxCS_Lagrange_EYsubX_normal_rcpp(arma::mat X,
                                       double alpha,
                                       arma::vec beta,
                                       double sigma,
                                       arma::vec phi,
                                       arma::vec CS_beta,
                                       arma::umat index,
                                       arma::vec eta){

  arma::uword n_n = X.n_rows;
  arma::uword n_m = phi.n_elem;
  double value = 0.0;
  arma::vec gradient(n_m);
  arma::mat hessian(n_m, n_m);
  List result;

  arma::vec X_intercept = arma::vec({1});

  for (size_t i = 0; i < n_n; ++i){

    const arma::vec Xrow_i = X.row(i).t();
    const arma::vec extXrow_i = arma::join_vert(X_intercept, Xrow_i);
    const double SI_i = alpha + arma::dot(Xrow_i, beta);
    const double eSI_CS_i = exp(arma::dot(Xrow_i, CS_beta));
    arma::vec Psi_i(n_m);

    for (size_t m = 0; m < n_m; ++m){

      if (index(i, m) == 1){

        Psi_i(m) += eSI_CS_i * (SI_i - phi(m));
      }
    }

    const double denominator = 1.0 + arma::dot(eta, Psi_i);
    if (denominator > 0.0){

      value += log(denominator);
      gradient += Psi_i / denominator;
      hessian += -Psi_i * Psi_i.t() / pow(denominator, 2);

    } else {

      value += -10000.0;
    }
  }

  result["value"] = value;
  result["gradient"] = gradient;
  result["hessian"] = hessian;
  return result;
}

// [[Rcpp::export]]
List auxCS_EY_normal_rcpp(arma::mat X,
                          double alpha,
                          arma::vec beta,
                          double sigma,
                          double phi,
                          arma::vec CS_beta){

  arma::uword n_n = X.n_rows;
  arma::uword n_p = X.n_cols;
  arma::vec Psi(1);
  arma::mat Psi_square(1, 1);
  arma::mat Psi_gradient(1, n_p * 2 + 3);
  List result;

  arma::vec X_intercept = arma::vec({1});

  for (size_t i = 0; i < n_n; ++i){

    const arma::vec Xrow_i = X.row(i).t();
    const arma::vec extXrow_i = arma::join_vert(X_intercept, Xrow_i);
    const double SI_i = alpha + arma::dot(Xrow_i, beta);
    const double eSI_CS_i = exp(arma::dot(Xrow_i, CS_beta));
    Psi += eSI_CS_i * (SI_i - phi);
    Psi_square += Psi * Psi.t();
    Psi_gradient.submat(0, 0, 0, n_p) += extXrow_i.t() * eSI_CS_i;
    Psi_gradient(0, n_p + 2) += -eSI_CS_i;
    Psi_gradient.submat(0, n_p + 3, 0, n_p * 2 + 2) += Xrow_i.t() * eSI_CS_i * (SI_i - phi);
  }

  Psi /= n_n;
  Psi_square /= n_n;
  Psi_gradient /= n_n;

  result["score"] = Psi;
  result["score_square"] = Psi_square;
  result["score_gradient"] = Psi_gradient;
  return result;
}

// [[Rcpp::export]]
List auxCS_EXsubY_normal_rcpp(arma::mat X,
                              double alpha,
                              arma::vec beta,
                              double sigma,
                              arma::mat phi,
                              arma::vec CS_beta,
                              arma::mat y_pts){

  arma::uword n_n = X.n_rows;
  arma::uword n_p = X.n_cols;
  arma::uword n_k = phi.n_rows;
  arma::uword n_m = n_p * n_k;
  arma::vec Psi(n_m);
  arma::mat Psi_square(n_m, n_m);
  arma::mat Psi_gradient(n_m, n_p * 2 + n_m + 2);
  List result;

  arma::vec X_intercept = arma::vec({1});

  for (size_t i = 0; i < n_n; ++i){

    const arma::vec Xrow_i = X.row(i).t();
    const arma::vec extXrow_i = arma::join_vert(X_intercept, Xrow_i);
    const double SI_i = alpha + arma::dot(Xrow_i, beta);
    const double eSI_CS_i = exp(arma::dot(Xrow_i, CS_beta));
    const arma::mat q_pts = (y_pts - SI_i) / sigma;
    const arma::mat cdf_pts = arma::normcdf(q_pts);
    const arma::mat pdf_pts = arma::normpdf(q_pts);
    const arma::vec cdf_dist = cdf_pts.col(1) - cdf_pts.col(0);
    const arma::vec pdf_dist = pdf_pts.col(1) - pdf_pts.col(0);
    const arma::mat x_phi_dist = (arma::repmat(Xrow_i.t(), n_k, 1) - phi);
    const arma::mat Psi_i = x_phi_dist % arma::repmat(cdf_dist, 1, n_p) * eSI_CS_i;
    Psi += arma::reshape(Psi_i.t(), n_m, 1);
    Psi_square += arma::reshape(Psi_i.t(), n_m, 1) *
      arma::reshape(Psi_i.t(), 1, n_m);
    Psi_gradient.cols(0, n_p) += -arma::reshape((x_phi_dist %
      arma::repmat(pdf_dist / sigma, 1, n_p)).t(), n_m, 1) *
      extXrow_i.t() * eSI_CS_i;
    Psi_gradient.col(n_p + 1) += -arma::reshape((x_phi_dist %
      arma::repmat((q_pts.col(1) / sigma) % pdf_pts.col(1) +
      (q_pts.col(0) / sigma) % pdf_pts.col(0), 1, n_p)).t(), n_m, 1) * eSI_CS_i;
    Psi_gradient.cols(n_p + 2, n_p + n_m + 1).diag() +=
      arma::reshape(arma::repmat(cdf_dist, 1, n_p).t(), n_m, 1) * (-eSI_CS_i);
    Psi_gradient.cols(n_p + n_m + 2, n_p * 2 + n_m + 1) += arma::reshape(Psi_i.t(), n_m, 1) * Xrow_i.t();
  }

  Psi /= n_n;
  Psi_square /= n_n;
  Psi_gradient /= n_n;

  result["score"] = Psi;
  result["score_square"] = Psi_square;
  result["score_gradient"] = Psi_gradient;
  return result;
}

// [[Rcpp::export]]
List auxCS_EYsubX_normal_rcpp(arma::mat X,
                              double alpha,
                              arma::vec beta,
                              double sigma,
                              arma::vec phi,
                              arma::vec CS_beta,
                              arma::umat index){

  arma::uword n_n = X.n_rows;
  arma::uword n_p = X.n_cols;
  arma::uword n_m = phi.n_elem;
  arma::vec Psi(n_m);
  arma::mat Psi_square(n_m, n_m);
  arma::mat Psi_gradient(n_m, n_p * 2 + n_m + 2);
  List result;

  arma::vec X_intercept = arma::vec({1});

  for (size_t i = 0; i < n_n; ++i){

    const arma::vec Xrow_i = X.row(i).t();
    const arma::vec extXrow_i = arma::join_vert(X_intercept, Xrow_i);
    const double SI_i = alpha + arma::dot(Xrow_i, beta);
    const double eSI_CS_i = exp(arma::dot(Xrow_i, CS_beta));
    arma::vec Psi_i(n_m);

    for (size_t m = 0; m < n_m; ++m){

      if (index(i, m) == 1){

        Psi_i(m) += eSI_CS_i * (SI_i - phi(m));
        Psi_gradient.submat(m, 0, m, n_p) += extXrow_i.t() * eSI_CS_i;
        Psi_gradient(m, n_p + 2 + m) += -eSI_CS_i;
        Psi_gradient.submat(m, n_p + n_m + 2, m, n_p * 2 + n_m + 1) += Xrow_i.t() *
          eSI_CS_i * (SI_i - phi(m));
      }
    }

    Psi += Psi_i;
    Psi_square += Psi * Psi.t();
  }

  Psi /= n_n;
  Psi_square /= n_n;
  Psi_gradient /= n_n;

  result["score"] = Psi;
  result["score_square"] = Psi_square;
  result["score_gradient"] = Psi_gradient;
  return result;
}



// [[Rcpp::export]]
List aux_Lagrange_EY_gamma_rcpp(arma::mat X,
                                double alpha,
                                arma::vec beta,
                                double nu,
                                double phi,
                                double eta){

  arma::uword n_n = X.n_rows;
  double value = 0.0;
  double gradient = 0.0;
  double hessian = 0.0;
  List result;

  arma::vec X_intercept = arma::vec({1});

  for (size_t i = 0; i < n_n; ++i){

    const arma::vec Xrow_i = X.row(i).t();
    const arma::vec extXrow_i = arma::join_vert(X_intercept, Xrow_i);
    const double eSI_i = exp(alpha + arma::dot(Xrow_i, beta));
    const double Psi_i = eSI_i - phi;

    const double denominator = 1.0 + eta * Psi_i;
    if (denominator > 0.0){

      value += log(denominator);
      gradient += Psi_i / denominator;
      hessian += -pow(Psi_i, 2) / pow(denominator, 2);

    } else {

      value += -10000.0;
    }
  }

  result["value"] = value;
  result["gradient"] = gradient;
  result["hessian"] = hessian;
  return result;
}

// [[Rcpp::export]]
List aux_Lagrange_EXsubY_gamma_rcpp(arma::mat X,
                                    double alpha,
                                    arma::vec beta,
                                    double nu,
                                    arma::mat phi,
                                    arma::mat y_pts,
                                    arma::vec eta){

  arma::uword n_n = X.n_rows;
  arma::uword n_p = X.n_cols;
  arma::uword n_k = phi.n_rows;
  arma::uword n_m = n_p * n_k;
  double value = 0.0;
  arma::vec gradient(n_m);
  arma::mat hessian(n_m, n_m);
  List result;

  arma::vec X_intercept = arma::vec({1});

  for (size_t i = 0; i < n_n; ++i){

    const arma::vec Xrow_i = X.row(i).t();
    const arma::vec extXrow_i = arma::join_vert(X_intercept, Xrow_i);
    const double eSI_i = exp(alpha + arma::dot(Xrow_i, beta));
    const double lambda_i = nu / eSI_i;
    const double lambda_i_inv = 1 / lambda_i;
    arma::mat Psi_i_mat(n_k, n_p);

    for (size_t k = 0; k < n_k; ++k){

      Psi_i_mat.row(k) = (Xrow_i.t() - phi.row(k)) *
        (R::pgamma(y_pts(k, 1), nu, lambda_i_inv, TRUE, FALSE) -
        R::pgamma(y_pts(k, 0), nu, lambda_i_inv, TRUE, FALSE));
    }

    const arma::vec Psi_i = arma::reshape(
      Psi_i_mat.t(), n_m, 1);

    const double denominator = 1.0 + arma::dot(eta, Psi_i);
    if (denominator > 0.0){

      value += log(denominator);
      gradient += Psi_i / denominator;
      hessian += -Psi_i * Psi_i.t() / pow(denominator, 2);

    } else {

      value += -10000.0;
    }
  }

  result["value"] = value;
  result["gradient"] = gradient;
  result["hessian"] = hessian;
  return result;
}

// [[Rcpp::export]]
List aux_Lagrange_EYsubX_gamma_rcpp(arma::mat X,
                                    double alpha,
                                    arma::vec beta,
                                    double nu,
                                    arma::vec phi,
                                    arma::umat index,
                                    arma::vec eta){

  arma::uword n_n = X.n_rows;
  arma::uword n_m = phi.n_elem;
  double value = 0.0;
  arma::vec gradient(n_m);
  arma::mat hessian(n_m, n_m);
  List result;

  arma::vec X_intercept = arma::vec({1});

  for (size_t i = 0; i < n_n; ++i){

    const arma::vec Xrow_i = X.row(i).t();
    const arma::vec extXrow_i = arma::join_vert(X_intercept, Xrow_i);
    const double eSI_i = exp(alpha + arma::dot(Xrow_i, beta));
    arma::vec Psi_i(n_m);

    for (size_t m = 0; m < n_m; ++m){

      if (index(i, m) == 1){

        Psi_i(m) += eSI_i - phi(m);
      }
    }

    const double denominator = 1.0 + arma::dot(eta, Psi_i);
    if (denominator > 0.0){

      value += log(denominator);
      gradient += Psi_i / denominator;
      hessian += -Psi_i * Psi_i.t() / pow(denominator, 2);

    } else {

      value += -10000.0;
    }
  }

  result["value"] = value;
  result["gradient"] = gradient;
  result["hessian"] = hessian;
  return result;
}



// [[Rcpp::export]]
List auxLS_Lagrange_EX_gamma_rcpp(arma::mat X,
                                  double alpha,
                                  arma::vec beta,
                                  double nu,
                                  arma::vec phi,
                                  double LS_beta,
                                  arma::vec eta){

  arma::uword n_n = X.n_rows;
  arma::uword n_p = X.n_cols;
  double value = 0.0;
  arma::vec gradient(n_p);
  arma::mat hessian(n_p, n_p);
  List result;

  arma::vec X_intercept = arma::vec({1});

  for (size_t i = 0; i < n_n; ++i){

    const arma::vec Xrow_i = X.row(i).t();
    const arma::vec extXrow_i = arma::join_vert(X_intercept, Xrow_i);
    const double eSI_i = exp(alpha + arma::dot(Xrow_i, beta));
    const double lambda_i = nu / eSI_i;
    const arma::vec Psi_i = (Xrow_i - phi) *
      pow(lambda_i / (lambda_i - LS_beta), nu);

    const double denominator = 1.0 + arma::dot(eta, Psi_i);
    if (denominator > 0.0){

      value += log(denominator);
      gradient += Psi_i / denominator;
      hessian += -Psi_i * Psi_i.t() / pow(denominator, 2);

    } else {

      value += -10000.0;
    }
  }

  result["value"] = value;
  result["gradient"] = gradient;
  result["hessian"] = hessian;
  return result;
}

// [[Rcpp::export]]
List auxLS_Lagrange_EY_gamma_rcpp(arma::mat X,
                                  double alpha,
                                  arma::vec beta,
                                  double nu,
                                  double phi,
                                  double LS_beta,
                                  double eta){

  arma::uword n_n = X.n_rows;
  double value = 0.0;
  double gradient = 0.0;
  double hessian = 0.0;
  List result;

  arma::vec X_intercept = arma::vec({1});

  for (size_t i = 0; i < n_n; ++i){

    const arma::vec Xrow_i = X.row(i).t();
    const arma::vec extXrow_i = arma::join_vert(X_intercept, Xrow_i);
    const double eSI_i = exp(alpha + arma::dot(Xrow_i, beta));
    const double lambda_i = nu / eSI_i;
    const double Psi_i = (nu / (lambda_i - LS_beta) - phi) *
      pow(lambda_i / (lambda_i - LS_beta), nu);

    const double denominator = 1.0 + eta * Psi_i;
    if (denominator > 0.0){

      value += log(denominator);
      gradient += Psi_i / denominator;
      hessian += -pow(Psi_i, 2) / pow(denominator, 2);

    } else {

      value += -10000.0;
    }
  }

  result["value"] = value;
  result["gradient"] = gradient;
  result["hessian"] = hessian;
  return result;
}

// [[Rcpp::export]]
List auxLS_Lagrange_EXsubY_gamma_rcpp(arma::mat X,
                                      double alpha,
                                      arma::vec beta,
                                      double nu,
                                      arma::mat phi,
                                      double LS_beta,
                                      arma::mat y_pts,
                                      arma::vec eta){

  arma::uword n_n = X.n_rows;
  arma::uword n_p = X.n_cols;
  arma::uword n_k = phi.n_rows;
  arma::uword n_m = n_p * n_k;
  double value = 0.0;
  arma::vec gradient(n_m);
  arma::mat hessian(n_m, n_m);
  List result;

  arma::vec X_intercept = arma::vec({1});

  for (size_t i = 0; i < n_n; ++i){

    const arma::vec Xrow_i = X.row(i).t();
    const arma::vec extXrow_i = arma::join_vert(X_intercept, Xrow_i);
    const double eSI_i = exp(alpha + arma::dot(Xrow_i, beta));
    const double lambda_i = nu / eSI_i;
    const double lambda_i_LS_beta = lambda_i - LS_beta;
    const double lambda_i_LS_beta_inv = 1 / lambda_i_LS_beta;
    const double ratio_i = pow(lambda_i / lambda_i_LS_beta, nu);
    arma::mat Psi_i_mat(n_k, n_p);

    for (size_t k = 0; k < n_k; ++k){

      Psi_i_mat.row(k) = (Xrow_i.t() - phi.row(k)) * ratio_i *
        (R::pgamma(y_pts(k, 1), nu, lambda_i_LS_beta_inv, TRUE, FALSE) -
        R::pgamma(y_pts(k, 0), nu, lambda_i_LS_beta_inv, TRUE, FALSE));
    }

    const arma::vec Psi_i = arma::reshape(
      Psi_i_mat.t(), n_m, 1);

    const double denominator = 1.0 + arma::dot(eta, Psi_i);
    if (denominator > 0.0){

      value += log(denominator);
      gradient += Psi_i / denominator;
      hessian += -Psi_i * Psi_i.t() / pow(denominator, 2);

    } else {

      value += -10000.0;
    }
  }

  result["value"] = value;
  result["gradient"] = gradient;
  result["hessian"] = hessian;
  return result;
}

// [[Rcpp::export]]
List auxLS_Lagrange_EYsubX_gamma_rcpp(arma::mat X,
                                      double alpha,
                                      arma::vec beta,
                                      double nu,
                                      arma::vec phi,
                                      double LS_beta,
                                      arma::umat index,
                                      arma::vec eta){

  arma::uword n_n = X.n_rows;
  arma::uword n_m = phi.n_elem;
  double value = 0.0;
  arma::vec gradient(n_m);
  arma::mat hessian(n_m, n_m);
  List result;

  arma::vec X_intercept = arma::vec({1});

  for (size_t i = 0; i < n_n; ++i){

    const arma::vec Xrow_i = X.row(i).t();
    const arma::vec extXrow_i = arma::join_vert(X_intercept, Xrow_i);
    const double eSI_i = exp(alpha + arma::dot(Xrow_i, beta));
    const double lambda_i = nu / eSI_i;
    arma::vec Psi_i(n_m);

    for (size_t m = 0; m < n_m; ++m){

      if (index(i, m) == 1){

        Psi_i(m) += (nu / (lambda_i - LS_beta) - phi(m)) *
          pow(lambda_i / (lambda_i - LS_beta), nu);
      }
    }

    const double denominator = 1.0 + arma::dot(eta, Psi_i);
    if (denominator > 0.0){

      value += log(denominator);
      gradient += Psi_i / denominator;
      hessian += -Psi_i * Psi_i.t() / pow(denominator, 2);

    } else {

      value += -10000.0;
    }
  }

  result["value"] = value;
  result["gradient"] = gradient;
  result["hessian"] = hessian;
  return result;
}



// [[Rcpp::export]]
List auxCS_Lagrange_EY_gamma_rcpp(arma::mat X,
                                   double alpha,
                                   arma::vec beta,
                                   double nu,
                                   double phi,
                                   arma::vec CS_beta,
                                   double eta){

  arma::uword n_n = X.n_rows;
  double value = 0.0;
  double gradient = 0.0;
  double hessian = 0.0;
  List result;

  arma::vec X_intercept = arma::vec({1});

  for (size_t i = 0; i < n_n; ++i){

    const arma::vec Xrow_i = X.row(i).t();
    const arma::vec extXrow_i = arma::join_vert(X_intercept, Xrow_i);
    const double eSI_i = exp(alpha + arma::dot(Xrow_i, beta));
    const double eSI_CS_i = exp(arma::dot(Xrow_i, CS_beta));
    const double Psi_i = eSI_CS_i * (eSI_i - phi);

    const double denominator = 1.0 + eta * Psi_i;
    if (denominator > 0.0){

      value += log(denominator);
      gradient += Psi_i / denominator;
      hessian += -pow(Psi_i, 2) / pow(denominator, 2);

    } else {

      value += -10000.0;
    }
  }

  result["value"] = value;
  result["gradient"] = gradient;
  result["hessian"] = hessian;
  return result;
}

// [[Rcpp::export]]
List auxCS_Lagrange_EXsubY_gamma_rcpp(arma::mat X,
                                      double alpha,
                                      arma::vec beta,
                                      double nu,
                                      arma::mat phi,
                                      arma::vec CS_beta,
                                      arma::mat y_pts,
                                      arma::vec eta){

  arma::uword n_n = X.n_rows;
  arma::uword n_p = X.n_cols;
  arma::uword n_k = phi.n_rows;
  arma::uword n_m = n_p * n_k;
  double value = 0.0;
  arma::vec gradient(n_m);
  arma::mat hessian(n_m, n_m);
  List result;

  arma::vec X_intercept = arma::vec({1});

  for (size_t i = 0; i < n_n; ++i){

    const arma::vec Xrow_i = X.row(i).t();
    const arma::vec extXrow_i = arma::join_vert(X_intercept, Xrow_i);
    const double eSI_i = exp(alpha + arma::dot(Xrow_i, beta));
    const double lambda_i = nu / eSI_i;
    const double lambda_i_inv = 1 / lambda_i;
    const double eSI_CS_i = exp(arma::dot(Xrow_i, CS_beta));
    arma::mat Psi_i_mat(n_k, n_p);

    for (size_t k = 0; k < n_k; ++k){

      Psi_i_mat.row(k) = (Xrow_i.t() - phi.row(k)) * eSI_CS_i *
        (R::pgamma(y_pts(k, 1), nu, lambda_i_inv, TRUE, FALSE) -
        R::pgamma(y_pts(k, 0), nu, lambda_i_inv, TRUE, FALSE));
    }

    const arma::vec Psi_i = arma::reshape(
      Psi_i_mat.t(), n_m, 1);

    const double denominator = 1.0 + arma::dot(eta, Psi_i);
    if (denominator > 0.0){

      value += log(denominator);
      gradient += Psi_i / denominator;
      hessian += -Psi_i * Psi_i.t() / pow(denominator, 2);

    } else {

      value += -10000.0;
    }
  }

  result["value"] = value;
  result["gradient"] = gradient;
  result["hessian"] = hessian;
  return result;
}

// [[Rcpp::export]]
List auxCS_Lagrange_EYsubX_gamma_rcpp(arma::mat X,
                                      double alpha,
                                      arma::vec beta,
                                      double nu,
                                      arma::vec phi,
                                      arma::vec CS_beta,
                                      arma::umat index,
                                      arma::vec eta){

  arma::uword n_n = X.n_rows;
  arma::uword n_m = phi.n_elem;
  double value = 0.0;
  arma::vec gradient(n_m);
  arma::mat hessian(n_m, n_m);
  List result;

  arma::vec X_intercept = arma::vec({1});

  for (size_t i = 0; i < n_n; ++i){

    const arma::vec Xrow_i = X.row(i).t();
    const arma::vec extXrow_i = arma::join_vert(X_intercept, Xrow_i);
    const double eSI_i = exp(alpha + arma::dot(Xrow_i, beta));
    const double eSI_CS_i = exp(arma::dot(Xrow_i, CS_beta));
    arma::vec Psi_i(n_m);

    for (size_t m = 0; m < n_m; ++m){

      if (index(i, m) == 1){

        Psi_i(m) += eSI_CS_i * (eSI_i - phi(m));
      }
    }

    const double denominator = 1.0 + arma::dot(eta, Psi_i);
    if (denominator > 0.0){

      value += log(denominator);
      gradient += Psi_i / denominator;
      hessian += -Psi_i * Psi_i.t() / pow(denominator, 2);

    } else {

      value += -10000.0;
    }
  }

  result["value"] = value;
  result["gradient"] = gradient;
  result["hessian"] = hessian;
  return result;
}








// [[Rcpp::export]]
double pgamma_rcpp(double c,
                   double shape,
                   double rate)
{
  double cdf = R::pgamma(c, shape, 1 / rate, TRUE, FALSE);
  return(cdf);
}












// [[Rcpp::export]]
List aux_EXsubgroupY_logistic_rcpp(arma::mat X,
                                   double alpha,
                                   arma::vec beta,
                                   arma::mat phi){

  arma::uword n_n = X.n_rows;
  arma::uword n_p = X.n_cols;
  arma::uword n_m = n_p * 2;
  arma::vec Psi(n_m);
  arma::mat Psi_square(n_m, n_m);
  arma::mat Psi_gradient(n_m, n_p + n_m + 1);
  List result;

  arma::vec X_intercept = arma::vec({1});

  for (size_t i = 0; i < n_n; ++i){

    const arma::vec Xrow_i = X.row(i).t();
    const arma::vec extXrow_i = arma::join_vert(X_intercept, Xrow_i);
    const double SI_i = alpha + arma::dot(Xrow_i, beta);
    const double e_f = exp(SI_i);
    const arma::vec p_vec = {e_f, 1};
    const arma::vec pn_vec = {1, -1};
    const arma::mat dist_X_phi = (arma::repmat(Xrow_i.t(), 2, 1) - phi);
    const arma::mat Psi_i = dist_X_phi % arma::repmat(p_vec, 1, n_p) / (e_f + 1.0);
    Psi += arma::reshape(Psi_i.t(), n_m, 1);
    Psi_square += arma::reshape(Psi_i.t(), n_m, 1) *
      arma::reshape(Psi_i.t(), 1, n_m);
    Psi_gradient.cols(0, n_p) += arma::reshape(
      (dist_X_phi % arma::repmat(pn_vec, 1, n_p)).t(), n_m, 1
    ) * e_f / pow(e_f + 1, 2) * extXrow_i.t();
    Psi_gradient.cols(n_p + 1, n_p + n_m).diag() += -arma::reshape(
      arma::repmat(p_vec, 1, n_p).t(), n_m, 1) / (e_f + 1.0);
  }

  Psi /= n_n;
  Psi_square /= n_n;
  Psi_gradient /= n_n;

  result["score"] = Psi;
  result["score_square"] = Psi_square;
  result["score_gradient"] = Psi_gradient;
  return result;
}

// [[Rcpp::export]]
List auxLS_EYsubgroupX_logistic_rcpp(arma::mat X,
                                     double alpha,
                                     arma::vec beta,
                                     arma::vec phi,
                                     double LS_beta,
                                     arma::umat index){

  arma::uword n_n = X.n_rows;
  arma::uword n_p = X.n_cols;
  arma::uword n_m = phi.n_elem;
  arma::vec Psi(n_m);
  arma::mat Psi_square(n_m, n_m);
  arma::mat Psi_gradient(n_m, n_p + n_m + 2);
  List result;

  arma::vec X_intercept = arma::vec({1});

  for (size_t i = 0; i < n_n; ++i){

    const arma::vec Xrow_i = X.row(i).t();
    const arma::vec extXrow_i = arma::join_vert(X_intercept, Xrow_i);
    const double SI_i = alpha + arma::dot(Xrow_i, beta);
    const double eSI_i = exp(SI_i);
    const double eSIb_i = eSI_i * exp(LS_beta);
    arma::vec Psi_i(n_m);

    for (size_t m = 0; m < n_m; ++m){

      if (index(i, m) == 1){

        Psi_i(m) += ((1 - phi(m)) * eSIb_i - phi(m)) / (1 + eSI_i);
        Psi_gradient.submat(m, 0, m, n_p) += extXrow_i.t() * eSI_i /
          pow(1 + eSI_i, 2) * ((1 - phi(m)) * exp(LS_beta) + phi(m));
        Psi_gradient(m, n_p + 1 + m) += -(1 + eSIb_i) / (1 + eSI_i);
        Psi_gradient(m, n_p + n_m + 1) += (1 - phi(m)) * eSIb_i / (1 + eSI_i);
      }
    }

    Psi += Psi_i;
    Psi_square += Psi * Psi.t();
  }

  Psi /= n_n;
  Psi_square /= n_n;
  Psi_gradient /= n_n;

  result["score"] = Psi;
  result["score_square"] = Psi_square;
  result["score_gradient"] = Psi_gradient;
  return result;
}

// [[Rcpp::export]]
List aux_EYsubgroupX_logistic_rcpp(arma::mat X,
                                   double alpha,
                                   arma::vec beta,
                                   arma::vec phi,
                                   arma::umat index){

  arma::uword n_n = X.n_rows;
  arma::uword n_p = X.n_cols;
  arma::uword n_m = phi.n_elem;
  arma::vec Psi(n_m);
  arma::mat Psi_square(n_m, n_m);
  arma::mat Psi_gradient(n_m, n_p + n_m + 1);
  List result;

  arma::vec X_intercept = arma::vec({1});

  for (size_t i = 0; i < n_n; ++i){

    const arma::vec Xrow_i = X.row(i).t();
    const arma::vec extXrow_i = arma::join_vert(X_intercept, Xrow_i);
    const double eSI_i = exp(alpha + arma::dot(Xrow_i, beta));
    arma::vec Psi_i(n_m);

    for (size_t m = 0; m < n_m; ++m){

      if (index(i, m) == 1){

        Psi_i(m) += ((1 - phi(m)) * eSI_i - phi(m)) / (1 + eSI_i);
        Psi_gradient.submat(m, 0, m, n_p) += extXrow_i.t() * eSI_i / pow(1 + eSI_i, 2);
        Psi_gradient(m, n_p + 1 + m) += -1;
      }
    }

    Psi += Psi_i;
    Psi_square += Psi * Psi.t();
  }

  Psi /= n_n;
  Psi_square /= n_n;
  Psi_gradient /= n_n;

  result["score"] = Psi;
  result["score_square"] = Psi_square;
  result["score_gradient"] = Psi_gradient;
  return result;
}
