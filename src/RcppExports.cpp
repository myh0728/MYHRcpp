// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// eXsq_rcpp
arma::mat eXsq_rcpp(arma::mat data_X);
RcppExport SEXP _MYHRcpp_eXsq_rcpp(SEXP data_XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type data_X(data_XSEXP);
    rcpp_result_gen = Rcpp::wrap(eXsq_rcpp(data_X));
    return rcpp_result_gen;
END_RCPP
}
// eXsq_w_rcpp
arma::mat eXsq_w_rcpp(arma::mat data_X, arma::vec weight);
RcppExport SEXP _MYHRcpp_eXsq_w_rcpp(SEXP data_XSEXP, SEXP weightSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type data_X(data_XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type weight(weightSEXP);
    rcpp_result_gen = Rcpp::wrap(eXsq_w_rcpp(data_X, weight));
    return rcpp_result_gen;
END_RCPP
}
// Xsq_lowtri_rcpp
arma::mat Xsq_lowtri_rcpp(arma::mat data_X);
RcppExport SEXP _MYHRcpp_Xsq_lowtri_rcpp(SEXP data_XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type data_X(data_XSEXP);
    rcpp_result_gen = Rcpp::wrap(Xsq_lowtri_rcpp(data_X));
    return rcpp_result_gen;
END_RCPP
}
// twoXYsym_lowtri_rcpp
arma::mat twoXYsym_lowtri_rcpp(arma::mat data_X, arma::mat data_Y);
RcppExport SEXP _MYHRcpp_twoXYsym_lowtri_rcpp(SEXP data_XSEXP, SEXP data_YSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type data_X(data_XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type data_Y(data_YSEXP);
    rcpp_result_gen = Rcpp::wrap(twoXYsym_lowtri_rcpp(data_X, data_Y));
    return rcpp_result_gen;
END_RCPP
}
// ctingP_rcpp
arma::mat ctingP_rcpp(arma::mat Y, arma::mat y);
RcppExport SEXP _MYHRcpp_ctingP_rcpp(SEXP YSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(ctingP_rcpp(Y, y));
    return rcpp_result_gen;
END_RCPP
}
// ctingP_uni_rcpp
arma::mat ctingP_uni_rcpp(arma::vec Y, arma::vec y);
RcppExport SEXP _MYHRcpp_ctingP_uni_rcpp(SEXP YSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(ctingP_uni_rcpp(Y, y));
    return rcpp_result_gen;
END_RCPP
}
// pinv_rcpp
arma::mat pinv_rcpp(arma::mat M_A);
RcppExport SEXP _MYHRcpp_pinv_rcpp(SEXP M_ASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type M_A(M_ASEXP);
    rcpp_result_gen = Rcpp::wrap(pinv_rcpp(M_A));
    return rcpp_result_gen;
END_RCPP
}
// solve_rcpp
arma::mat solve_rcpp(arma::mat M_A, arma::mat M_B);
RcppExport SEXP _MYHRcpp_solve_rcpp(SEXP M_ASEXP, SEXP M_BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type M_A(M_ASEXP);
    Rcpp::traits::input_parameter< arma::mat >::type M_B(M_BSEXP);
    rcpp_result_gen = Rcpp::wrap(solve_rcpp(M_A, M_B));
    return rcpp_result_gen;
END_RCPP
}
// inv_sympd_rcpp
arma::mat inv_sympd_rcpp(arma::mat M_S);
RcppExport SEXP _MYHRcpp_inv_sympd_rcpp(SEXP M_SSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type M_S(M_SSEXP);
    rcpp_result_gen = Rcpp::wrap(inv_sympd_rcpp(M_S));
    return rcpp_result_gen;
END_RCPP
}
// eigen_rcpp
List eigen_rcpp(arma::mat M_S);
RcppExport SEXP _MYHRcpp_eigen_rcpp(SEXP M_SSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type M_S(M_SSEXP);
    rcpp_result_gen = Rcpp::wrap(eigen_rcpp(M_S));
    return rcpp_result_gen;
END_RCPP
}
// GroupSum_rcpp
arma::mat GroupSum_rcpp(arma::mat MM, arma::uvec id);
RcppExport SEXP _MYHRcpp_GroupSum_rcpp(SEXP MMSEXP, SEXP idSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type MM(MMSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type id(idSEXP);
    rcpp_result_gen = Rcpp::wrap(GroupSum_rcpp(MM, id));
    return rcpp_result_gen;
END_RCPP
}
// countAinB_rcpp
arma::vec countAinB_rcpp(arma::vec A, arma::vec B);
RcppExport SEXP _MYHRcpp_countAinB_rcpp(SEXP ASEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::vec >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(countAinB_rcpp(A, B));
    return rcpp_result_gen;
END_RCPP
}
// countAinB_W_rcpp
arma::vec countAinB_W_rcpp(arma::vec A, arma::vec B, arma::vec W);
RcppExport SEXP _MYHRcpp_countAinB_W_rcpp(SEXP ASEXP, SEXP BSEXP, SEXP WSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::vec >::type B(BSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type W(WSEXP);
    rcpp_result_gen = Rcpp::wrap(countAinB_W_rcpp(A, B, W));
    return rcpp_result_gen;
END_RCPP
}
// rankAinB_rcpp
arma::vec rankAinB_rcpp(arma::vec A, arma::vec B);
RcppExport SEXP _MYHRcpp_rankAinB_rcpp(SEXP ASEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::vec >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(rankAinB_rcpp(A, B));
    return rcpp_result_gen;
END_RCPP
}
// atRisk_integral_rcpp
arma::mat atRisk_integral_rcpp(arma::mat integrand, arma::vec t_start, arma::vec t_stop, arma::vec t_event);
RcppExport SEXP _MYHRcpp_atRisk_integral_rcpp(SEXP integrandSEXP, SEXP t_startSEXP, SEXP t_stopSEXP, SEXP t_eventSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type integrand(integrandSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type t_start(t_startSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type t_stop(t_stopSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type t_event(t_eventSEXP);
    rcpp_result_gen = Rcpp::wrap(atRisk_integral_rcpp(integrand, t_start, t_stop, t_event));
    return rcpp_result_gen;
END_RCPP
}
// sum_atRisk_rcpp
arma::mat sum_atRisk_rcpp(arma::mat summand, arma::vec t_start, arma::vec t_stop, arma::vec t_event);
RcppExport SEXP _MYHRcpp_sum_atRisk_rcpp(SEXP summandSEXP, SEXP t_startSEXP, SEXP t_stopSEXP, SEXP t_eventSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type summand(summandSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type t_start(t_startSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type t_stop(t_stopSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type t_event(t_eventSEXP);
    rcpp_result_gen = Rcpp::wrap(sum_atRisk_rcpp(summand, t_start, t_stop, t_event));
    return rcpp_result_gen;
END_RCPP
}
// KDE_KG_rcpp
arma::vec KDE_KG_rcpp(arma::mat X, arma::mat x, arma::vec h);
RcppExport SEXP _MYHRcpp_KDE_KG_rcpp(SEXP XSEXP, SEXP xSEXP, SEXP hSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type h(hSEXP);
    rcpp_result_gen = Rcpp::wrap(KDE_KG_rcpp(X, x, h));
    return rcpp_result_gen;
END_RCPP
}
// KDE_KG_w_rcpp
arma::vec KDE_KG_w_rcpp(arma::mat X, arma::mat x, arma::vec h, arma::vec w);
RcppExport SEXP _MYHRcpp_KDE_KG_w_rcpp(SEXP XSEXP, SEXP xSEXP, SEXP hSEXP, SEXP wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type h(hSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type w(wSEXP);
    rcpp_result_gen = Rcpp::wrap(KDE_KG_w_rcpp(X, x, h, w));
    return rcpp_result_gen;
END_RCPP
}
// KDEcv_KG_rcpp
arma::vec KDEcv_KG_rcpp(arma::mat X, arma::vec h);
RcppExport SEXP _MYHRcpp_KDEcv_KG_rcpp(SEXP XSEXP, SEXP hSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type h(hSEXP);
    rcpp_result_gen = Rcpp::wrap(KDEcv_KG_rcpp(X, h));
    return rcpp_result_gen;
END_RCPP
}
// KDEcv_KG_w_rcpp
arma::vec KDEcv_KG_w_rcpp(arma::mat X, arma::vec h, arma::vec w);
RcppExport SEXP _MYHRcpp_KDEcv_KG_w_rcpp(SEXP XSEXP, SEXP hSEXP, SEXP wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type h(hSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type w(wSEXP);
    rcpp_result_gen = Rcpp::wrap(KDEcv_KG_w_rcpp(X, h, w));
    return rcpp_result_gen;
END_RCPP
}
// KDE_K2B_rcpp
arma::vec KDE_K2B_rcpp(arma::mat X, arma::mat x, arma::vec h);
RcppExport SEXP _MYHRcpp_KDE_K2B_rcpp(SEXP XSEXP, SEXP xSEXP, SEXP hSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type h(hSEXP);
    rcpp_result_gen = Rcpp::wrap(KDE_K2B_rcpp(X, x, h));
    return rcpp_result_gen;
END_RCPP
}
// KDE_K2B_w_rcpp
arma::vec KDE_K2B_w_rcpp(arma::mat X, arma::mat x, arma::vec h, arma::vec w);
RcppExport SEXP _MYHRcpp_KDE_K2B_w_rcpp(SEXP XSEXP, SEXP xSEXP, SEXP hSEXP, SEXP wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type h(hSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type w(wSEXP);
    rcpp_result_gen = Rcpp::wrap(KDE_K2B_w_rcpp(X, x, h, w));
    return rcpp_result_gen;
END_RCPP
}
// KDEcv_K2B_rcpp
arma::vec KDEcv_K2B_rcpp(arma::mat X, arma::vec h);
RcppExport SEXP _MYHRcpp_KDEcv_K2B_rcpp(SEXP XSEXP, SEXP hSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type h(hSEXP);
    rcpp_result_gen = Rcpp::wrap(KDEcv_K2B_rcpp(X, h));
    return rcpp_result_gen;
END_RCPP
}
// KDEcv_K2B_w_rcpp
arma::vec KDEcv_K2B_w_rcpp(arma::mat X, arma::vec h, arma::vec w);
RcppExport SEXP _MYHRcpp_KDEcv_K2B_w_rcpp(SEXP XSEXP, SEXP hSEXP, SEXP wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type h(hSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type w(wSEXP);
    rcpp_result_gen = Rcpp::wrap(KDEcv_K2B_w_rcpp(X, h, w));
    return rcpp_result_gen;
END_RCPP
}
// KDE_K4B_rcpp
arma::vec KDE_K4B_rcpp(arma::mat X, arma::mat x, arma::vec h);
RcppExport SEXP _MYHRcpp_KDE_K4B_rcpp(SEXP XSEXP, SEXP xSEXP, SEXP hSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type h(hSEXP);
    rcpp_result_gen = Rcpp::wrap(KDE_K4B_rcpp(X, x, h));
    return rcpp_result_gen;
END_RCPP
}
// KDE_K4B_w_rcpp
arma::vec KDE_K4B_w_rcpp(arma::mat X, arma::mat x, arma::vec h, arma::vec w);
RcppExport SEXP _MYHRcpp_KDE_K4B_w_rcpp(SEXP XSEXP, SEXP xSEXP, SEXP hSEXP, SEXP wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type h(hSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type w(wSEXP);
    rcpp_result_gen = Rcpp::wrap(KDE_K4B_w_rcpp(X, x, h, w));
    return rcpp_result_gen;
END_RCPP
}
// KDEcv_K4B_rcpp
arma::vec KDEcv_K4B_rcpp(arma::mat X, arma::vec h);
RcppExport SEXP _MYHRcpp_KDEcv_K4B_rcpp(SEXP XSEXP, SEXP hSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type h(hSEXP);
    rcpp_result_gen = Rcpp::wrap(KDEcv_K4B_rcpp(X, h));
    return rcpp_result_gen;
END_RCPP
}
// KDEcv_K4B_w_rcpp
arma::vec KDEcv_K4B_w_rcpp(arma::mat X, arma::vec h, arma::vec w);
RcppExport SEXP _MYHRcpp_KDEcv_K4B_w_rcpp(SEXP XSEXP, SEXP hSEXP, SEXP wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type h(hSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type w(wSEXP);
    rcpp_result_gen = Rcpp::wrap(KDEcv_K4B_w_rcpp(X, h, w));
    return rcpp_result_gen;
END_RCPP
}
// NW_K2B_rcpp
arma::mat NW_K2B_rcpp(arma::mat X, arma::mat Y, arma::mat x, arma::vec h);
RcppExport SEXP _MYHRcpp_NW_K2B_rcpp(SEXP XSEXP, SEXP YSEXP, SEXP xSEXP, SEXP hSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type h(hSEXP);
    rcpp_result_gen = Rcpp::wrap(NW_K2B_rcpp(X, Y, x, h));
    return rcpp_result_gen;
END_RCPP
}
// NW_K2B_w_rcpp
arma::mat NW_K2B_w_rcpp(arma::mat X, arma::mat Y, arma::mat x, arma::vec h, arma::vec w);
RcppExport SEXP _MYHRcpp_NW_K2B_w_rcpp(SEXP XSEXP, SEXP YSEXP, SEXP xSEXP, SEXP hSEXP, SEXP wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type h(hSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type w(wSEXP);
    rcpp_result_gen = Rcpp::wrap(NW_K2B_w_rcpp(X, Y, x, h, w));
    return rcpp_result_gen;
END_RCPP
}
// NWcv_K2B_rcpp
arma::mat NWcv_K2B_rcpp(arma::mat X, arma::mat Y, arma::vec h);
RcppExport SEXP _MYHRcpp_NWcv_K2B_rcpp(SEXP XSEXP, SEXP YSEXP, SEXP hSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type h(hSEXP);
    rcpp_result_gen = Rcpp::wrap(NWcv_K2B_rcpp(X, Y, h));
    return rcpp_result_gen;
END_RCPP
}
// NWcv_K2B_w_rcpp
arma::mat NWcv_K2B_w_rcpp(arma::mat X, arma::mat Y, arma::vec h, arma::vec w);
RcppExport SEXP _MYHRcpp_NWcv_K2B_w_rcpp(SEXP XSEXP, SEXP YSEXP, SEXP hSEXP, SEXP wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type h(hSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type w(wSEXP);
    rcpp_result_gen = Rcpp::wrap(NWcv_K2B_w_rcpp(X, Y, h, w));
    return rcpp_result_gen;
END_RCPP
}
// NW_K4B_rcpp
arma::mat NW_K4B_rcpp(arma::mat X, arma::mat Y, arma::mat x, arma::vec h);
RcppExport SEXP _MYHRcpp_NW_K4B_rcpp(SEXP XSEXP, SEXP YSEXP, SEXP xSEXP, SEXP hSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type h(hSEXP);
    rcpp_result_gen = Rcpp::wrap(NW_K4B_rcpp(X, Y, x, h));
    return rcpp_result_gen;
END_RCPP
}
// NW_K4B_w_rcpp
arma::mat NW_K4B_w_rcpp(arma::mat X, arma::mat Y, arma::mat x, arma::vec h, arma::vec w);
RcppExport SEXP _MYHRcpp_NW_K4B_w_rcpp(SEXP XSEXP, SEXP YSEXP, SEXP xSEXP, SEXP hSEXP, SEXP wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type h(hSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type w(wSEXP);
    rcpp_result_gen = Rcpp::wrap(NW_K4B_w_rcpp(X, Y, x, h, w));
    return rcpp_result_gen;
END_RCPP
}
// NWcv_K4B_rcpp
arma::mat NWcv_K4B_rcpp(arma::mat X, arma::mat Y, arma::vec h);
RcppExport SEXP _MYHRcpp_NWcv_K4B_rcpp(SEXP XSEXP, SEXP YSEXP, SEXP hSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type h(hSEXP);
    rcpp_result_gen = Rcpp::wrap(NWcv_K4B_rcpp(X, Y, h));
    return rcpp_result_gen;
END_RCPP
}
// NWcv_K4B_w_rcpp
arma::mat NWcv_K4B_w_rcpp(arma::mat X, arma::mat Y, arma::vec h, arma::vec w);
RcppExport SEXP _MYHRcpp_NWcv_K4B_w_rcpp(SEXP XSEXP, SEXP YSEXP, SEXP hSEXP, SEXP wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type h(hSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type w(wSEXP);
    rcpp_result_gen = Rcpp::wrap(NWcv_K4B_w_rcpp(X, Y, h, w));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_hello
List rcpp_hello();
RcppExport SEXP _MYHRcpp_rcpp_hello() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(rcpp_hello());
    return rcpp_result_gen;
END_RCPP
}
// KME_rcpp
arma::vec KME_rcpp(arma::vec t_stop, arma::uvec is_event, arma::vec t_event);
RcppExport SEXP _MYHRcpp_KME_rcpp(SEXP t_stopSEXP, SEXP is_eventSEXP, SEXP t_eventSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type t_stop(t_stopSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type is_event(is_eventSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type t_event(t_eventSEXP);
    rcpp_result_gen = Rcpp::wrap(KME_rcpp(t_stop, is_event, t_event));
    return rcpp_result_gen;
END_RCPP
}
// KME_w_rcpp
arma::vec KME_w_rcpp(arma::vec t_stop, arma::uvec is_event, arma::vec t_event, arma::vec w);
RcppExport SEXP _MYHRcpp_KME_w_rcpp(SEXP t_stopSEXP, SEXP is_eventSEXP, SEXP t_eventSEXP, SEXP wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type t_stop(t_stopSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type is_event(is_eventSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type t_event(t_eventSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type w(wSEXP);
    rcpp_result_gen = Rcpp::wrap(KME_w_rcpp(t_stop, is_event, t_event, w));
    return rcpp_result_gen;
END_RCPP
}
// SKME_K2B_rcpp
arma::mat SKME_K2B_rcpp(arma::vec t_stop, arma::uvec is_event, arma::vec t_event, arma::mat X, arma::mat x, arma::vec h);
RcppExport SEXP _MYHRcpp_SKME_K2B_rcpp(SEXP t_stopSEXP, SEXP is_eventSEXP, SEXP t_eventSEXP, SEXP XSEXP, SEXP xSEXP, SEXP hSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type t_stop(t_stopSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type is_event(is_eventSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type t_event(t_eventSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type h(hSEXP);
    rcpp_result_gen = Rcpp::wrap(SKME_K2B_rcpp(t_stop, is_event, t_event, X, x, h));
    return rcpp_result_gen;
END_RCPP
}
// SKME_K2B_w_rcpp
arma::mat SKME_K2B_w_rcpp(arma::vec t_stop, arma::uvec is_event, arma::vec t_event, arma::mat X, arma::mat x, arma::vec h, arma::vec w);
RcppExport SEXP _MYHRcpp_SKME_K2B_w_rcpp(SEXP t_stopSEXP, SEXP is_eventSEXP, SEXP t_eventSEXP, SEXP XSEXP, SEXP xSEXP, SEXP hSEXP, SEXP wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type t_stop(t_stopSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type is_event(is_eventSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type t_event(t_eventSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type h(hSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type w(wSEXP);
    rcpp_result_gen = Rcpp::wrap(SKME_K2B_w_rcpp(t_stop, is_event, t_event, X, x, h, w));
    return rcpp_result_gen;
END_RCPP
}
// SKME_K4B_rcpp
arma::mat SKME_K4B_rcpp(arma::vec t_stop, arma::uvec is_event, arma::vec t_event, arma::mat X, arma::mat x, arma::vec h);
RcppExport SEXP _MYHRcpp_SKME_K4B_rcpp(SEXP t_stopSEXP, SEXP is_eventSEXP, SEXP t_eventSEXP, SEXP XSEXP, SEXP xSEXP, SEXP hSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type t_stop(t_stopSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type is_event(is_eventSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type t_event(t_eventSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type h(hSEXP);
    rcpp_result_gen = Rcpp::wrap(SKME_K4B_rcpp(t_stop, is_event, t_event, X, x, h));
    return rcpp_result_gen;
END_RCPP
}
// SKME_K4B_w_rcpp
arma::mat SKME_K4B_w_rcpp(arma::vec t_stop, arma::uvec is_event, arma::vec t_event, arma::mat X, arma::mat x, arma::vec h, arma::vec w);
RcppExport SEXP _MYHRcpp_SKME_K4B_w_rcpp(SEXP t_stopSEXP, SEXP is_eventSEXP, SEXP t_eventSEXP, SEXP XSEXP, SEXP xSEXP, SEXP hSEXP, SEXP wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type t_stop(t_stopSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type is_event(is_eventSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type t_event(t_eventSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type h(hSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type w(wSEXP);
    rcpp_result_gen = Rcpp::wrap(SKME_K4B_w_rcpp(t_stop, is_event, t_event, X, x, h, w));
    return rcpp_result_gen;
END_RCPP
}
// testfunction
arma::vec testfunction(arma::vec a, Rcpp::Function my_r_func);
RcppExport SEXP _MYHRcpp_testfunction(SEXP aSEXP, SEXP my_r_funcSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type a(aSEXP);
    Rcpp::traits::input_parameter< Rcpp::Function >::type my_r_func(my_r_funcSEXP);
    rcpp_result_gen = Rcpp::wrap(testfunction(a, my_r_func));
    return rcpp_result_gen;
END_RCPP
}
// KDE_rcpp_kernel
arma::vec KDE_rcpp_kernel(arma::mat X, arma::mat x, Rcpp::Function K, arma::vec h);
RcppExport SEXP _MYHRcpp_KDE_rcpp_kernel(SEXP XSEXP, SEXP xSEXP, SEXP KSEXP, SEXP hSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< Rcpp::Function >::type K(KSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type h(hSEXP);
    rcpp_result_gen = Rcpp::wrap(KDE_rcpp_kernel(X, x, K, h));
    return rcpp_result_gen;
END_RCPP
}
// KDE_K2B_rcpp_chatgpt
arma::vec KDE_K2B_rcpp_chatgpt(const arma::mat& X, const arma::mat& x, const arma::vec& h);
RcppExport SEXP _MYHRcpp_KDE_K2B_rcpp_chatgpt(SEXP XSEXP, SEXP xSEXP, SEXP hSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type h(hSEXP);
    rcpp_result_gen = Rcpp::wrap(KDE_K2B_rcpp_chatgpt(X, x, h));
    return rcpp_result_gen;
END_RCPP
}
// KDEcv_K2B_rcpp_o1
arma::vec KDEcv_K2B_rcpp_o1(arma::mat X, arma::vec h);
RcppExport SEXP _MYHRcpp_KDEcv_K2B_rcpp_o1(SEXP XSEXP, SEXP hSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type h(hSEXP);
    rcpp_result_gen = Rcpp::wrap(KDEcv_K2B_rcpp_o1(X, h));
    return rcpp_result_gen;
END_RCPP
}
// KDEcv_K2B_w_rcpp_o1
arma::vec KDEcv_K2B_w_rcpp_o1(arma::mat X, arma::vec h, arma::vec w);
RcppExport SEXP _MYHRcpp_KDEcv_K2B_w_rcpp_o1(SEXP XSEXP, SEXP hSEXP, SEXP wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type h(hSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type w(wSEXP);
    rcpp_result_gen = Rcpp::wrap(KDEcv_K2B_w_rcpp_o1(X, h, w));
    return rcpp_result_gen;
END_RCPP
}
// KDE_K4B_rcpp_o1
arma::vec KDE_K4B_rcpp_o1(arma::mat X, arma::mat x, arma::vec h);
RcppExport SEXP _MYHRcpp_KDE_K4B_rcpp_o1(SEXP XSEXP, SEXP xSEXP, SEXP hSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type h(hSEXP);
    rcpp_result_gen = Rcpp::wrap(KDE_K4B_rcpp_o1(X, x, h));
    return rcpp_result_gen;
END_RCPP
}
// KDEcv_K4B_rcpp_o1
arma::vec KDEcv_K4B_rcpp_o1(arma::mat X, arma::vec h);
RcppExport SEXP _MYHRcpp_KDEcv_K4B_rcpp_o1(SEXP XSEXP, SEXP hSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type h(hSEXP);
    rcpp_result_gen = Rcpp::wrap(KDEcv_K4B_rcpp_o1(X, h));
    return rcpp_result_gen;
END_RCPP
}
// NWF_K2B_rcpp
arma::mat NWF_K2B_rcpp(arma::mat X, arma::mat Y, arma::mat x, arma::mat y, arma::vec h);
RcppExport SEXP _MYHRcpp_NWF_K2B_rcpp(SEXP XSEXP, SEXP YSEXP, SEXP xSEXP, SEXP ySEXP, SEXP hSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::vec >::type h(hSEXP);
    rcpp_result_gen = Rcpp::wrap(NWF_K2B_rcpp(X, Y, x, y, h));
    return rcpp_result_gen;
END_RCPP
}
// NWcv_K2B_rcpp_o1
arma::mat NWcv_K2B_rcpp_o1(arma::mat X, arma::mat Y, arma::vec h);
RcppExport SEXP _MYHRcpp_NWcv_K2B_rcpp_o1(SEXP XSEXP, SEXP YSEXP, SEXP hSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type h(hSEXP);
    rcpp_result_gen = Rcpp::wrap(NWcv_K2B_rcpp_o1(X, Y, h));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_MYHRcpp_eXsq_rcpp", (DL_FUNC) &_MYHRcpp_eXsq_rcpp, 1},
    {"_MYHRcpp_eXsq_w_rcpp", (DL_FUNC) &_MYHRcpp_eXsq_w_rcpp, 2},
    {"_MYHRcpp_Xsq_lowtri_rcpp", (DL_FUNC) &_MYHRcpp_Xsq_lowtri_rcpp, 1},
    {"_MYHRcpp_twoXYsym_lowtri_rcpp", (DL_FUNC) &_MYHRcpp_twoXYsym_lowtri_rcpp, 2},
    {"_MYHRcpp_ctingP_rcpp", (DL_FUNC) &_MYHRcpp_ctingP_rcpp, 2},
    {"_MYHRcpp_ctingP_uni_rcpp", (DL_FUNC) &_MYHRcpp_ctingP_uni_rcpp, 2},
    {"_MYHRcpp_pinv_rcpp", (DL_FUNC) &_MYHRcpp_pinv_rcpp, 1},
    {"_MYHRcpp_solve_rcpp", (DL_FUNC) &_MYHRcpp_solve_rcpp, 2},
    {"_MYHRcpp_inv_sympd_rcpp", (DL_FUNC) &_MYHRcpp_inv_sympd_rcpp, 1},
    {"_MYHRcpp_eigen_rcpp", (DL_FUNC) &_MYHRcpp_eigen_rcpp, 1},
    {"_MYHRcpp_GroupSum_rcpp", (DL_FUNC) &_MYHRcpp_GroupSum_rcpp, 2},
    {"_MYHRcpp_countAinB_rcpp", (DL_FUNC) &_MYHRcpp_countAinB_rcpp, 2},
    {"_MYHRcpp_countAinB_W_rcpp", (DL_FUNC) &_MYHRcpp_countAinB_W_rcpp, 3},
    {"_MYHRcpp_rankAinB_rcpp", (DL_FUNC) &_MYHRcpp_rankAinB_rcpp, 2},
    {"_MYHRcpp_atRisk_integral_rcpp", (DL_FUNC) &_MYHRcpp_atRisk_integral_rcpp, 4},
    {"_MYHRcpp_sum_atRisk_rcpp", (DL_FUNC) &_MYHRcpp_sum_atRisk_rcpp, 4},
    {"_MYHRcpp_KDE_KG_rcpp", (DL_FUNC) &_MYHRcpp_KDE_KG_rcpp, 3},
    {"_MYHRcpp_KDE_KG_w_rcpp", (DL_FUNC) &_MYHRcpp_KDE_KG_w_rcpp, 4},
    {"_MYHRcpp_KDEcv_KG_rcpp", (DL_FUNC) &_MYHRcpp_KDEcv_KG_rcpp, 2},
    {"_MYHRcpp_KDEcv_KG_w_rcpp", (DL_FUNC) &_MYHRcpp_KDEcv_KG_w_rcpp, 3},
    {"_MYHRcpp_KDE_K2B_rcpp", (DL_FUNC) &_MYHRcpp_KDE_K2B_rcpp, 3},
    {"_MYHRcpp_KDE_K2B_w_rcpp", (DL_FUNC) &_MYHRcpp_KDE_K2B_w_rcpp, 4},
    {"_MYHRcpp_KDEcv_K2B_rcpp", (DL_FUNC) &_MYHRcpp_KDEcv_K2B_rcpp, 2},
    {"_MYHRcpp_KDEcv_K2B_w_rcpp", (DL_FUNC) &_MYHRcpp_KDEcv_K2B_w_rcpp, 3},
    {"_MYHRcpp_KDE_K4B_rcpp", (DL_FUNC) &_MYHRcpp_KDE_K4B_rcpp, 3},
    {"_MYHRcpp_KDE_K4B_w_rcpp", (DL_FUNC) &_MYHRcpp_KDE_K4B_w_rcpp, 4},
    {"_MYHRcpp_KDEcv_K4B_rcpp", (DL_FUNC) &_MYHRcpp_KDEcv_K4B_rcpp, 2},
    {"_MYHRcpp_KDEcv_K4B_w_rcpp", (DL_FUNC) &_MYHRcpp_KDEcv_K4B_w_rcpp, 3},
    {"_MYHRcpp_NW_K2B_rcpp", (DL_FUNC) &_MYHRcpp_NW_K2B_rcpp, 4},
    {"_MYHRcpp_NW_K2B_w_rcpp", (DL_FUNC) &_MYHRcpp_NW_K2B_w_rcpp, 5},
    {"_MYHRcpp_NWcv_K2B_rcpp", (DL_FUNC) &_MYHRcpp_NWcv_K2B_rcpp, 3},
    {"_MYHRcpp_NWcv_K2B_w_rcpp", (DL_FUNC) &_MYHRcpp_NWcv_K2B_w_rcpp, 4},
    {"_MYHRcpp_NW_K4B_rcpp", (DL_FUNC) &_MYHRcpp_NW_K4B_rcpp, 4},
    {"_MYHRcpp_NW_K4B_w_rcpp", (DL_FUNC) &_MYHRcpp_NW_K4B_w_rcpp, 5},
    {"_MYHRcpp_NWcv_K4B_rcpp", (DL_FUNC) &_MYHRcpp_NWcv_K4B_rcpp, 3},
    {"_MYHRcpp_NWcv_K4B_w_rcpp", (DL_FUNC) &_MYHRcpp_NWcv_K4B_w_rcpp, 4},
    {"_MYHRcpp_rcpp_hello", (DL_FUNC) &_MYHRcpp_rcpp_hello, 0},
    {"_MYHRcpp_KME_rcpp", (DL_FUNC) &_MYHRcpp_KME_rcpp, 3},
    {"_MYHRcpp_KME_w_rcpp", (DL_FUNC) &_MYHRcpp_KME_w_rcpp, 4},
    {"_MYHRcpp_SKME_K2B_rcpp", (DL_FUNC) &_MYHRcpp_SKME_K2B_rcpp, 6},
    {"_MYHRcpp_SKME_K2B_w_rcpp", (DL_FUNC) &_MYHRcpp_SKME_K2B_w_rcpp, 7},
    {"_MYHRcpp_SKME_K4B_rcpp", (DL_FUNC) &_MYHRcpp_SKME_K4B_rcpp, 6},
    {"_MYHRcpp_SKME_K4B_w_rcpp", (DL_FUNC) &_MYHRcpp_SKME_K4B_w_rcpp, 7},
    {"_MYHRcpp_testfunction", (DL_FUNC) &_MYHRcpp_testfunction, 2},
    {"_MYHRcpp_KDE_rcpp_kernel", (DL_FUNC) &_MYHRcpp_KDE_rcpp_kernel, 4},
    {"_MYHRcpp_KDE_K2B_rcpp_chatgpt", (DL_FUNC) &_MYHRcpp_KDE_K2B_rcpp_chatgpt, 3},
    {"_MYHRcpp_KDEcv_K2B_rcpp_o1", (DL_FUNC) &_MYHRcpp_KDEcv_K2B_rcpp_o1, 2},
    {"_MYHRcpp_KDEcv_K2B_w_rcpp_o1", (DL_FUNC) &_MYHRcpp_KDEcv_K2B_w_rcpp_o1, 3},
    {"_MYHRcpp_KDE_K4B_rcpp_o1", (DL_FUNC) &_MYHRcpp_KDE_K4B_rcpp_o1, 3},
    {"_MYHRcpp_KDEcv_K4B_rcpp_o1", (DL_FUNC) &_MYHRcpp_KDEcv_K4B_rcpp_o1, 2},
    {"_MYHRcpp_NWF_K2B_rcpp", (DL_FUNC) &_MYHRcpp_NWF_K2B_rcpp, 5},
    {"_MYHRcpp_NWcv_K2B_rcpp_o1", (DL_FUNC) &_MYHRcpp_NWcv_K2B_rcpp_o1, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_MYHRcpp(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
