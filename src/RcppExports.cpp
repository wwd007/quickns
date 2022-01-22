// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// make_design2
NumericMatrix make_design2(int nx, int ncoef, int ord, NumericMatrix temp);
RcppExport SEXP _quickns_make_design2(SEXP nxSEXP, SEXP ncoefSEXP, SEXP ordSEXP, SEXP tempSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type nx(nxSEXP);
    Rcpp::traits::input_parameter< int >::type ncoef(ncoefSEXP);
    Rcpp::traits::input_parameter< int >::type ord(ordSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type temp(tempSEXP);
    rcpp_result_gen = Rcpp::wrap(make_design2(nx, ncoef, ord, temp));
    return rcpp_result_gen;
END_RCPP
}
// make_design_no_intcpt_cpp
NumericMatrix make_design_no_intcpt_cpp(int nx, int ncoef, int ord, NumericMatrix temp);
RcppExport SEXP _quickns_make_design_no_intcpt_cpp(SEXP nxSEXP, SEXP ncoefSEXP, SEXP ordSEXP, SEXP tempSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type nx(nxSEXP);
    Rcpp::traits::input_parameter< int >::type ncoef(ncoefSEXP);
    Rcpp::traits::input_parameter< int >::type ord(ordSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type temp(tempSEXP);
    rcpp_result_gen = Rcpp::wrap(make_design_no_intcpt_cpp(nx, ncoef, ord, temp));
    return rcpp_result_gen;
END_RCPP
}
// get_basis
Eigen::MatrixXd get_basis(const Eigen::Map<Eigen::MatrixXd> const_m, const Eigen::Map<Eigen::MatrixXd> basis_m);
RcppExport SEXP _quickns_get_basis(SEXP const_mSEXP, SEXP basis_mSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type const_m(const_mSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type basis_m(basis_mSEXP);
    rcpp_result_gen = Rcpp::wrap(get_basis(const_m, basis_m));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP spline_basis(SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_quickns_make_design2", (DL_FUNC) &_quickns_make_design2, 4},
    {"_quickns_make_design_no_intcpt_cpp", (DL_FUNC) &_quickns_make_design_no_intcpt_cpp, 4},
    {"_quickns_get_basis", (DL_FUNC) &_quickns_get_basis, 2},
    {"spline_basis", (DL_FUNC) &spline_basis, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_quickns(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
