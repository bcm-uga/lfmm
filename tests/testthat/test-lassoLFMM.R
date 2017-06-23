library(testthat)
context("lassoLFMM")

test_that("lassoLFMM_heuristic_gamma and lassoLFMM_heuristic_lambda_range", {

  K <- 3
  dat <- lfmm_sampler(n = 100, p = 1000, K = K,
                      outlier.prop = 0.1,
                      cs = c(0.8),
                      sigma = 0.2,
                      B.sd = 1.0,
                      U.sd = 1.0,
                      V.sd = 1.0)


  m <- lassoLFMM(K = 3)

  params <- lassoLFMM_heuristic_gamma_lambda_range(m, dat)

  expect_equal(length(params$lambda.range), 100)

})

test_that("compute_soft_svd", {

  skip("To debug")
  ## Eigen::Map<Eigen::MatrixXd> U,
  ## Eigen::Map<Eigen::MatrixXd> V


  D_thau <- function(X, gamma) {
    m <- list()
    svd.res <- svd(X, nu = 0, nv = 0) # compute only singular value
    aux <- svd.res$d - gamma
    Sigma <- diag(aux[aux > 0.0])
    m$K <- ncol(Sigma)

    svd.res <- svd(X, nu = m$K, nv = m$K)
    m$U <- svd.res$u %*% Sigma
    m$V <- svd.res$v
    m
  }

  K <- 3
  dat <- lfmm_sampler(n = 100, p = 1000, K = K,
                      outlier.prop = 0.1,
                      cs = c(0.8),
                      sigma = 0.2,
                      B.sd = 1.0,
                      U.sd = 1.0,
                      V.sd = 1.0)

  ## compute gamma
  m <- lassoLFMM(K = 3)
  params <- lassoLFMM_heuristic_gamma_lambda_range(m, dat)

  ## c++
  U <- matrix(0, 100, 3)
  V <- matrix(0, 1000, 3)
  compute_soft_SVD(dat$Y,params$gamma,U,V)
  W.cpp <- tcrossprod(U,V)

  ## R
  res <- D_thau(dat$Y, params$gamma)
  dim(res$U)
  W.r <- tcrossprod(res$U,res$V)

  ## comp
  mean(abs(W.r - W.cpp))

})


test_that("lassoLFMM_main", {

  K <- 3
  dat <- lfmm_sampler(n = 100, p = 1000, K = K,
                      outlier.prop = 0.1,
                      cs = c(0.8),
                      sigma = 0.2,
                      B.sd = 1.0,
                      U.sd = 1.0,
                      V.sd = 1.0)
  dat$G <- dat$Y

  ## compute gamma
  m <- lassoLFMM(K = 3,
                 nozero.prop = 0.1,
                 lambda.K = 20,
                 lambda.eps = 0.001)
  params <- lassoLFMM_heuristic_gamma_lambda_range(m, dat)
  
  ## init B
  if (is.null(m$B)) {
    m$B <- matrix(0.0, ncol(dat$Y), ncol(dat$X))
  }

  ## init U and V
  if (is.null(m$U)) {
    m$U <- matrix(0.0, nrow(dat$Y), m$K)
  }
  if (is.null(m$V)) {
    m$V <- matrix(0.0, ncol(dat$Y), m$K)
  }

  ## main loop
  lambda <- params$lambda.range[20]
  relative.err.epsilon = 1e-6
  it.max <- 100
  
  ## c++
  res.cpp <- lassoLFMM_main(dat$Y, dat$X,
                            params$gamma, lambda,
                            relative.err.epsilon,
                            it.max,
                            m$U,
                            m$V,
                            m$B)
  ## why err decrease and then increase ??
  ## m <- res.cpp
  ## res.cpp <- lassoLFMM_main(dat$Y, dat$X,
  ##                           params$gamma, lambda,
  ##                           relative.err.epsilon,
  ##                           it.max,
  ##                           m$U,
  ##                           m$V,
  ##                           m$B)
  ## ok it not recompute all
  expect_equal(dim(res.cpp$U), c(100, 3))
  expect_equal(dim(res.cpp$V), c(1000, 3))
  expect_equal(dim(res.cpp$B), c(1000, 1))


  ## R implementation
  res.rr <- lassoLFMM_main_R(dat$Y, dat$X,
                            params$gamma, lambda,
                            relative.err.epsilon,
                            it.max,
                            m$U,
                            m$V,
                            m$B)

  ## ThesisRpackage implementation
  skip_if_not_installed("ThesisRpackage")
  futile.logger::flog.threshold(futile.logger::TRACE, name = "ThesisRpackage")
  res.r <- ThesisRpackage::LassoLFMMMethod(K = NULL,
                                           gamma = params$gamma,
                                           it.max = it.max,
                                           err.max = relative.err.epsilon,
                                           lambda = lambda,
                                           center = FALSE)
  res.r <- ThesisRpackage::fit(res.r , dat)

  ## comp
  W.r <- tcrossprod(res.r$U, res.r$V)
  W.cpp <- tcrossprod(res.cpp$U, res.cpp$V)
  expect_lte(mean(abs(W.r - W.cpp)), 1e-10)
  expect_lte(mean(abs(res.cpp$B - t(res.r$B))), 1e-10)
  expect_equal(mean(res.cpp$B != 0), mean(res.r$B != 0))

})

test_that("lassoLFMM", {

  K <- 3
  dat <- lfmm_sampler(n = 100, p = 1000, K = K,
                      outlier.prop = 0.1,
                      cs = c(0.8),
                      sigma = 0.2,
                      B.sd = 1.0,
                      U.sd = 1.0,
                      V.sd = 1.0)
  dat$G <- dat$Y

  ## lassoLFMM
  m <- lassoLFMM(K = 3,
                 nozero.prop = 0.1,
                 lambda.K = 20,
                 lambda.eps = 0.001)
  m <- MatrixFactorizationR_fit(m, dat,
                                it.max = 100, relative.err.epsilon = 1e-4)

  skip_if_not_installed("ThesisRpackage")
  futile.logger::flog.threshold(futile.logger::TRACE, name = "ThesisRpackage")
  m.ThesisRpackage <- ThesisRpackage::finalLfmmLassoMethod(K = 3, sparse.prop = 0.1,
                                                           lambda.K = 20, lambda.eps = 0.001)
  m.ThesisRpackage$center <- FALSE
  m.ThesisRpackage <- ThesisRpackage::fit(m.ThesisRpackage, dat)


  ## comp
  W <- tcrossprod(m$U, m$V)
  W.thesis <- tcrossprod(m.ThesisRpackage$U, m.ThesisRpackage$V)
  expect_lte(mean(abs(W - W.thesis)), 1e-10)
  mean(m$B != 0)
  mean(m.ThesisRpackage$B != 0)
  expect_lte(mean(abs(m$B - t(m.ThesisRpackage$B))), 1e-10)

})

test_that("lassoLFMM", {

  K <- 3
  n <- 100
  p <- 1000
  dat <- lfmm_sampler(n = n, p = p, K = K,
                      outlier.prop = 0.1,
                      cs = c(0.8),
                      sigma = 0.2,
                      B.sd = 1.0,
                      U.sd = 1.0,
                      V.sd = 1.0)

  ## no NA
  lfmm.noNA <- lassoLFMM(K = 3, nozero.prop = 0.1)
  lfmm.noNA <- MatrixFactorizationR_fit(lfmm.noNA, dat)

  ## add na
  na.ind <- sample.int(n * p, 0.01 * n * p)
  dat$Y[na.ind] <- NA

  ## lfmm with na
  lfmm.NA <- lassoLFMM(K = 3, nozero.prop = 0.1)
  lfmm.NA <- MatrixFactorizationR_fit(lfmm.NA, dat)

  ## impute by median first
  dat$Y <- impute_median(dat$Y)
  lfmm.NA.impute <- lassoLFMM(K = 3, nozero.prop = 0.1)
  lfmm.NA.impute <- MatrixFactorizationR_fit(lfmm.NA.impute, dat)

  ## comparison W
  W.NA <- tcrossprod(lfmm.NA$U, lfmm.NA$V)
  W.noNA <- tcrossprod(lfmm.noNA$U, lfmm.noNA$V)
  W.NA.impute <- tcrossprod(lfmm.NA.impute$U, lfmm.NA.impute$V)
  e1 <- sqrt(mean((W.NA - W.noNA) ^ 2))
  e2 <- sqrt(mean((W.NA.impute - W.noNA) ^ 2))
  expect_gt((e2 - e1) / e1, 0.1)

  ## comparison B
  e1 <- sqrt(mean((lfmm.noNA$B - lfmm.NA$B) ^ 2))
  e2 <- sqrt(mean((lfmm.noNA$B - lfmm.NA.impute$B) ^ 2))
  expect_gt((e2 - e1) / e1, 0.1)

})
