// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>

// via the depends attribute we tell Rcpp to create hooks for
// RcppEigen so that the build process will know what to do
//
// [[Rcpp::depends(RcppEigen)]]

// three internal util functions
using Rcpp::NumericMatrix;
using Rcpp::IntegerMatrix;
using Rcpp::IntegerVector;
using Rcpp::seq_len;
using Rcpp::rep_each;
using Rcpp::outer;
// [[Rcpp::export]]
NumericMatrix make_design2(int nx, int ncoef, int ord, NumericMatrix temp) {
    IntegerVector offset = temp.attr("Offsets");
    IntegerVector ii = rep_each(seq_len(nx) - 1, ord);
    IntegerMatrix jj = outer(seq_len(ord) - 1, offset, std::plus<int>());
    NumericMatrix design(nx, ncoef);
    for (int i = 0; i < ii.size(); ++i) {
        design(ii[i], jj[i]) = temp[i];
    }
    return design;
}

// [[Rcpp::export]]
NumericMatrix make_design_no_intcpt_cpp(int nx, int ncoef, int ord, NumericMatrix temp) {
    IntegerVector offset = temp.attr("Offsets");
    IntegerVector ii = rep_each(seq_len(nx) - 1, ord);
    IntegerMatrix jj = outer(seq_len(ord) - 2, offset, std::plus<int>());
    NumericMatrix design(nx, ncoef - 1);
    for (int i = 0; i < ii.size(); ++i) {
        if (jj[i] >= 0) {
            design(ii[i], jj[i]) = temp[i];
        }
    }
    return design;
}

// [[Rcpp::export]]
Eigen::MatrixXd get_basis(const Eigen::Map<Eigen::MatrixXd> const_m, const Eigen::Map<Eigen::MatrixXd> basis_m) {
    const Eigen::HouseholderQR<Eigen::MatrixXd> QR(const_m.adjoint());
    const Eigen::MatrixXd A((QR.householderQ().adjoint() * basis_m.adjoint()).adjoint());
    return A;
}