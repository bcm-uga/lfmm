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
Rcpp::List compute_eigen_svd(const Eigen::Map<Eigen::MatrixXd> & X) {

  // compute svd of X
  JacobiSVD<MatrixXd> svd(X, ComputeFullU | ComputeFullV);
  return Rcpp::List::create(Named("Q") = svd.matrixU(),
                            Named("R") = svd.matrixV(),
                            Named("sigma") = svd.singularValues()
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

// [[Rcpp::export]]
Eigen::VectorXd err2s_lfmm_cpp(const Eigen::Map<Eigen::MatrixXd> Y,
                              const Eigen::Map<Eigen::MatrixXd> X,
                              const Eigen::Map<Eigen::MatrixXd> U,
                              const Eigen::Map<Eigen::MatrixXd> V,
                              const Eigen::Map<Eigen::MatrixXd> B) {
  // constants
  const int n = Y.rows();
  const int p = Y.cols();
  Eigen::VectorXd err2s = VectorXd(p); 
  double aux = 0.0;

  for (int j = 0; j < p; j++) {
    err2s(j) = 0.0;
    for (int i = 0; i < n; i++) {
      aux = Y(i, j) -
        U.row(i).dot(V.row(j)) -
        X.row(i).dot(B.row(j));
      err2s(j) += aux * aux;
    }
  }
  return(err2s);
}

// [[Rcpp::export]]
Eigen::VectorXd sum2_lm_cpp(const Eigen::Map<Eigen::MatrixXd> Y,
                            const Eigen::Map<Eigen::MatrixXd> X,
                            const Eigen::Map<Eigen::MatrixXd> B) {
  // constants
  const int n = Y.rows();
  const int p = Y.cols();
  VectorXd err2 = VectorXd::Zero(p);
  double aux = 0.0;

  for (int j = 0; j < p; j++) {
    err2(j) = 0.0;
    for (int i = 0; i < n; i++) {
      aux = Y(i, j) -
        X.row(i).dot(B.row(j));
      err2(j) += aux * aux;
    }
  }
  return(err2);
}
