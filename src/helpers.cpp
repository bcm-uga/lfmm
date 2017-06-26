// -*- mode: poly-c++r -*-

#include <algorithm>
#include <cmath>
#include <RcppEigen.h>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;
using namespace Eigen;

// [[Rcpp::export]]
Rcpp::List compute_P(const Eigen::Map<Eigen::MatrixXd> & X,
                          const double lambda) {

  // constants
  const int n = X.rows();
  const int d = X.cols();

  // compute oblique projection matrix
  MatrixXd In =  MatrixXd::Identity(n,n);
  MatrixXd D =  MatrixXd::Identity(d,d);
  LDLT<MatrixXd> lltX = (X.transpose() * X + lambda * D).ldlt();
  MatrixXd P = In - X * lltX.solve(X.transpose());

  // compute of eigen decomposition
  SelfAdjointEigenSolver<MatrixXd> es(P);

  return Rcpp::List::create(Named("sqrt.P") = es.operatorSqrt(),
                            Named("sqrt.P.inv") = es.operatorInverseSqrt()
                            // Named("P") = P
                            );
}


