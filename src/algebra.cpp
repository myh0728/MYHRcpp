#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat pinv_rcpp(arma::mat M_A){

  arma::mat M_X = arma::pinv(M_A);
  return M_X;
}

// [[Rcpp::export]]
arma::mat solve_rcpp(arma::mat M_A, arma::mat M_B){

  arma::mat M_X = arma::solve(M_A, M_B);
  return M_X;
}

// [[Rcpp::export]]
arma::mat inv_sympd_rcpp(arma::mat M_S){

  arma::mat M_S_inv = arma::inv_sympd(M_S);
  return M_S_inv;
}

// [[Rcpp::export]]
List eigen_rcpp(arma::mat M_S){

  List result;
  arma::vec eigval;
  arma::mat eigvec;

  eig_sym(eigval, eigvec, M_S);
  eigval = arma::reverse(eigval);
  eigvec = arma::reverse(eigvec, 1);

  result["value"] = eigval;
  result["vector"] = eigvec;
  return result;
}

// [[Rcpp::export]]
arma::mat GroupSum_rcpp(arma::mat MM,
                        arma::uvec id){

  arma::uvec id_unique = arma::unique(id);
  arma::uword n_group = id_unique.n_elem;
  arma::uword n_obs = MM.n_rows;
  arma::uword n_p = MM.n_cols;
  arma::mat result(n_group, n_p);

  for (size_t i = 0; i < n_obs; ++i){

    for (size_t j = 0; j < n_group; ++j){

      if (id(i) == id_unique(j)){

        result.row(j) += MM.row(i);
      }
    }
  }

  return result;
}

// [[Rcpp::export]]
arma::vec countAinB_rcpp(arma::vec A,
                         arma::vec B){

  arma::uword l_A = A.n_elem;
  arma::uword l_B = B.n_elem;
  arma::vec result(l_A);

  for (size_t i = 0; i < l_A; ++i){

    for (size_t j = 0; j < l_B; ++j){

      if (A(i) == B(j)){

        result(i) += 1;
      }
    }
  }

  return result;
}

// [[Rcpp::export]]
arma::vec countAinB_W_rcpp(arma::vec A,
                           arma::vec B,
                           arma::vec W){

  arma::uword l_A = A.n_elem;
  arma::uword l_B = B.n_elem;
  arma::vec result(l_A);

  for (size_t i = 0; i < l_A; ++i){

    for (size_t j = 0; j < l_B; ++j){

      if (A(i) == B(j)){

        result(i) += W(j);
      }
    }
  }

  return result;
}

// [[Rcpp::export]]
arma::vec rankAinB_rcpp(arma::vec A,
                        arma::vec B){

  arma::uword l_A = A.n_elem;
  arma::uword l_B = B.n_elem;
  arma::vec result(l_A);

  for (size_t i = 0; i < l_A; ++i){

    for (size_t j = 0; j < l_B; ++j){

      if (A(i) >= B(j)){

        result(i) += 1;

      }else{

        break;
      }
    }
  }

  return result;
}

// [[Rcpp::export]]
arma::mat outer_minus_rcpp(arma::vec a,
                           arma::vec b){

  arma::uword l_a = a.n_elem;
  arma::uword l_b = b.n_elem;
  arma::mat outer(l_a, l_b);

  for (size_t i = 0; i < l_a; ++i){

    for (size_t j = 0; j < l_b; ++j){

      outer(i, j) = a(i) - b(j);
    }
  }

  return outer;
}

// [[Rcpp::export]]
arma::mat outer_plus_rcpp(arma::vec a,
                          arma::vec b){

  arma::uword l_a = a.n_elem;
  arma::uword l_b = b.n_elem;
  arma::mat outer(l_a, l_b);

  for (size_t i = 0; i < l_a; ++i){

    for (size_t j = 0; j < l_b; ++j){

      outer(i, j) = a(i) + b(j);
    }
  }

  return outer;
}

// [[Rcpp::export]]
arma::mat outer_leq_rcpp(arma::vec a,
                         arma::vec b){

  arma::uword l_a = a.n_elem;
  arma::uword l_b = b.n_elem;
  arma::mat outer(l_a, l_b);

  for (size_t i = 0; i < l_a; ++i){

    for (size_t j = 0; j < l_b; ++j){

      outer(i, j) = (a(i) <= b(j));
    }
  }

  return outer;
}

// [[Rcpp::export]]
arma::mat outer_geq_rcpp(arma::vec a,
                         arma::vec b){

  arma::uword l_a = a.n_elem;
  arma::uword l_b = b.n_elem;
  arma::mat outer(l_a, l_b);

  for (size_t i = 0; i < l_a; ++i){

    for (size_t j = 0; j < l_b; ++j){

      outer(i, j) = (a(i) >= b(j));
    }
  }

  return outer;
}


