#ifndef HELPERS_H
#define HELPERS_H

Eigen::MatrixXd compute_B_ridge(const Eigen::MatrixXd & Y,
                                const Eigen::MatrixXd & X,
                                const double lambda);

Eigen::MatrixXd compute_B_lasso(const Eigen::MatrixXd & Y,
                                const Eigen::MatrixXd & X,
                                const double lambda);

void compute_K_SVD(const Eigen::MatrixXd & Y,
                   const int K,
                   Eigen::MatrixXd & U,
                   Eigen::MatrixXd & V);

void compute_soft_SVD(const Eigen::MatrixXd & Y,
                      const double gamma,
                      Eigen::MatrixXd & U,
                      Eigen::MatrixXd & V);

#endif
