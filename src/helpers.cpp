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


// // [[Rcpp::export]]
// Eigen::MatrixXd compute_B_lasso(const Eigen::MatrixXd & Y,
//                                 const Eigen::MatrixXd & X,
//                                 const double lambda) {
//   // RMK: we assume that co variates are orthogonal
//   MatrixXd B_hat = compute_B_ridge(Y, X, 0.0);
//   MatrixXd B = B_hat;
//   for (int i = 0; i < B.rows(); ++i) {
//     for (int j = 0; j < B.cols(); ++j) {
//       B(i, j) = std::abs(B_hat(i, j)) - lambda;
//       B(i, j) = std::max(B(i, j), 0.0);
//       B(i, j) = copysign(B(i, j), B_hat(i, j));
//     }
//   }
//   return B;
// }


void compute_K_SVD(const Eigen::MatrixXd & Y,
                   const int K,
                   Eigen::MatrixXd & U,
                   Eigen::MatrixXd & V) {
  // svd
  JacobiSVD<MatrixXd> svd(Y, ComputeThinU | ComputeThinV);

  // compute of U
  VectorXd svs = svd.singularValues().head(K);
  U = svd.matrixU().leftCols(K) * svs.asDiagonal();

  // compute of V
  V = svd.matrixV().leftCols(K);
}

// [[Rcpp::export]]
void compute_soft_SVD(const Eigen::MatrixXd & Y,
                      const double gamma,
                      Eigen::MatrixXd & U,
                      Eigen::MatrixXd & V) {
  // svd
  JacobiSVD<MatrixXd> svd(Y, ComputeThinU | ComputeThinV);

  // compute K
  ArrayXd svs = svd.singularValues().array();
  int K = (svs >= gamma).cast<int>().sum(); // True give zero...
  VectorXd svs_thresholded = (svs.head(K) - gamma).matrix(); // BUG.... if i write un svs ...???

  // compute of U
  U = svd.matrixU().leftCols(K) * svs_thresholded.asDiagonal();

  // compute of V
  V = svd.matrixV().leftCols(K);
}
