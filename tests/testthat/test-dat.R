library(testthat)
context("Dat wrapper")

test_that("Dat", {
  Y <- "../Data/1000Genomes/Phase3/European_Chrm22.maf.05.sample.10000.rds"
  skip_if_not(file.exists(Y))

  dat <- Dat(Y = Y)

  expect_equal(dim(dat$Y), c(503, 10000))

})

test_that("LfmmDat", {

  Y <- "../Data/1000Genomes/Phase3/European_Chrm22.maf.05.sample.10000.rds"
  skip_if_not(file.exists(Y))

  X <- matrix(rnorm(503, 503, 1))

  dat <- LfmmDat(Y = Y, X = X)

  expect_equal(dim(dat$Y), c(503, 10000))
  expect_equal(dim(dat$X), c(503, 1))
})

test_that("Dat and Rspectra", {

  n <- 100
  m <- 1000
  k <- 3
  dat <- lfmm_sampler(n = 100, p = 1000, K = k,
                      outlier.prop = 0.1,
                      cs = c(0.8),
                      sigma = 0.2,
                      B.sd = 1.0,
                      U.sd = 1.0,
                      V.sd = 1.0)

  Af <- function(x, args) {
    args$dat$productY(x)
  }
  Atransf <- function(x, args) {
    args$dat$productYt(x)
  }
  res.rspectra <- RSpectra::svds(A = Af,
                                 k, nu = k, nv = k,
                                 Atrans = Atransf, dim = c(n, m), opts = list(tol = 1e-10),
                                 args = list(dat = dat))
  res.svd <- svd(dat$Y, k, k)

  expect_lt(mean(abs(res.rspectra$u - res.svd$u)), 1)
  expect_lt(mean(abs(res.rspectra$d - res.svd$d[1:k])), 1e-10)
  expect_lt(mean(abs(res.rspectra$v - res.svd$v)), 1)
  W.svd <- tcrossprod(res.svd$u %*% diag(res.svd$d[1:k]), res.svd$v)
  W.rspectra <- tcrossprod(res.rspectra$u %*% diag(res.rspectra$d[1:k]), res.rspectra$v)
  expect_lt(mean(abs(W.rspectra - W.svd)), 1e-10)
  ## error because au PC get same variance

})

test_that("LfmmDat impute and err2", {

  dat <- lfmm_sampler(n = 100, p = 1000, K = 3,
                      outlier.prop = 0.1,
                      cs = c(0.8),
                      sigma = 0.2,
                      B.sd = 1.0,
                      U.sd = 1.0,
                      V.sd = 1.0)
  dat$U <- NULL
  dat$V <- NULL
  dat$B <- NULL
  ## run lfmm ridge
  m <- ridgeLFMM(K = 3, 1e-5)
  m <- MatrixFactorizationR_fit(m, dat)

  ## NA
  prop <- 0.01
  n <- nrow(dat$Y)
  p <- ncol(dat$Y)
  dat$missing.ind <- sample(n * p, prop * n * p)
  dat$Y[dat$missing.ind] <- NA

  ## impute
  Y <- dat$Y
  Yux <- tcrossprod(dat$X, m$B) + tcrossprod(m$U, m$V)
  Y[dat$missing.ind]  <- Yux[dat$missing.ind]
  dat$impute_lfmm(m$U, m$V, m$B)
  anyNA(dat$Y)
  expect_lte(mean(abs(dat$Y - Y)), 1e-18)

  ## err2
  Yux <- tcrossprod(dat$X, m$B) + tcrossprod(m$U, m$V)
  Yux <- Y - Yux
  err2.R <- mean(Yux ^ 2)
  err2.cpp <- dat$err2_lfmm(m$U, m$V, m$B)
  expect_lte(mean(abs(err2.cpp - err2.R)), 1e-10)

  ## sum2_lm
  E <- dat$Y - tcrossprod(dat$X, m$B) 
  effective.degree.freedom <- 99
  epsilon.sigma2 <- apply(E, 2, function(x) sum(x ^ 2) / effective.degree.freedom)
  s2.cpp <- dat$sigma2_lm(dat$X, m$B, effective.degree.freedom)
  expect_lte(mean(abs(epsilon.sigma2 - s2.cpp)), 1e-10)

})

test_that("Dat svd", {

  n <- 100
  m <- 1000
  k <- 3
  dat <- lfmm_sampler(n = 100, p = 1000, K = k,
                      outlier.prop = 0.1,
                      cs = c(0.8),
                      sigma = 0.2,
                      B.sd = 1.0,
                      U.sd = 1.0,
                      V.sd = 1.0)
  Af <- function(x, args) {
    dat$productY(x)
  }
  Atransf <- function(x, args) {
    dat$productYt(x)
  }
  res.rspectra <- compute_svd(Af, Atransf, k, k, k, dim = c(nrow(dat$Y), ncol(dat$Y)))
  res.svd <- svd(dat$Y, k, k)

  expect_lt(mean(abs(res.rspectra$u - res.svd$u)), 1)
  expect_lt(mean(abs(res.rspectra$d - res.svd$d[1:k])), 1e-10)
  expect_lt(mean(abs(res.rspectra$v - res.svd$v)), 1)
  W.svd <- tcrossprod(res.svd$u %*% diag(res.svd$d[1:k]), res.svd$v)
  W.rspectra <- tcrossprod(res.rspectra$u %*% diag(res.rspectra$d[1:k]), res.rspectra$v)
  ## error because au PC get same variance

  expect_lt(mean(abs(W.rspectra - W.svd)), 1e-10)

})

