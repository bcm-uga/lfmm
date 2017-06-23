// -*- mode: poly-c++r -*-

#include <cmath> 
#include <limits>
#include <RcppEigen.h>

using namespace Rcpp;
using namespace Eigen;

#include "helpers.h"

// [[Rcpp::export]]
Rcpp::List lassoLFMM_main(const Eigen::Map<Eigen::MatrixXd> Y,
                          const Eigen::Map<Eigen::MatrixXd> X,
                          const double gamma,
                          const double lambda,
                          const double relative_err_epsilon,
                          const int it_max,
                          const Eigen::Map<Eigen::MatrixXd> U0,
                          const Eigen::Map<Eigen::MatrixXd> V0,
                          const Eigen::Map<Eigen::MatrixXd> B0) {
  // constants
  const int n = Y.rows();
  const int p = Y.cols();

  // variables
  MatrixXd U = U0;
  MatrixXd V = V0;
  MatrixXd Yux = Y;
  MatrixXd B = B0;
  double err = 0.0;
  Yux = Y - U * V.transpose();
  Yux = Yux - X * B.transpose();
  double err_new = Yux.squaredNorm() / n / p;
  double relative_err = std::numeric_limits<double>::max();
  int it = 1;

  // main loop
  while ((it <= it_max) && (relative_err > relative_err_epsilon)) {
    err = err_new;
    Rcpp::Rcout << "It = " << it << "/" << it_max << ", err2 = " << err << std::endl;

    // step B
    Yux = Y - U * V.transpose();
    B = compute_B_lasso(Yux, X, lambda);

    // compute W = UV^T
    Yux = Y - X * B.transpose();
    compute_soft_SVD(Yux, gamma, U, V);


    // err
    Yux = Yux - U * V.transpose();
    err_new = Yux.squaredNorm() / n / p;
    relative_err = std::abs(err_new - err) / err;
    it++;
  }

  return Rcpp::List::create(
                            Named("U") = U,
                            Named("V") = V,
                            Named("B") = B
                            );
}

