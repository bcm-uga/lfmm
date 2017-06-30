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


// [[Rcpp::export]]
void impute_lfmm_cpp(Eigen::Map<Eigen::MatrixXd> Y,
                     const Eigen::Map<Eigen::MatrixXd> X,
                     const Eigen::Map<Eigen::MatrixXd> U,
                     const Eigen::Map<Eigen::MatrixXd> V,
                     const Eigen::Map<Eigen::MatrixXd> B,
                     NumericVector missingId
                     ) {

  // constants
  const int n = Y.rows();
  int k = 0;
  int i = 0;
  int j = 0;

  for (int t = 0; t < missingId.size(); t++) {
    k = missingId[t] - 1;
    i = k % n;
    j = k / n;
    Y(i, j) = U.row(i).dot(V.row(j)) +
      X.row(i).dot(B.row(j));
  }
}

// [[Rcpp::export]]
double err2_lfmm_cpp(const Eigen::Map<Eigen::MatrixXd> Y,
                     const Eigen::Map<Eigen::MatrixXd> X,
                     const Eigen::Map<Eigen::MatrixXd> U,
                     const Eigen::Map<Eigen::MatrixXd> V,
                     const Eigen::Map<Eigen::MatrixXd> B) {
  // constants
  const int n = Y.rows();
  const int p = Y.cols();
  double err2 = 0.0;
  double aux = 0.0;

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < p; j++) {
      aux = Y(i, j) -
        U.row(i).dot(V.row(j)) -
        X.row(i).dot(B.row(j));
      err2 += aux * aux;
    }
  }
  return(err2 / n / p);
}

