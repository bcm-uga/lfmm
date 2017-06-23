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
    args$productY(x)
  }
  Atransf <- function(x, args) {
    args$productYt(x)
  }
  res.rspectra <- RSpectra::svds(A = Af,
                                 k, nu = k, nv = k,
                                 Atrans = Atransf, dim = c(n, m), opts = list(tol = 1e-10))
  res.svd <- svd(dat$Y, k, k)

  expect_lt(mean(abs(res.rspectra$u - res.svd$u)), 1)
  expect_lt(mean(abs(res.rspectra$d - res.svd$d[1:k])), 1e-10)
  expect_lt(mean(abs(res.rspectra$v - res.svd$v)), 1)
  W.svd <- tcrossprod(res.svd$u %*% diag(res.svd$d[1:k]), res.svd$v)
  W.rspectra <- tcrossprod(res.rspectra$u %*% diag(res.rspectra$d[1:k]), res.rspectra$v)
  expect_lt(mean(abs(W.rspectra - W.svd)), 1e-10)
  ## error because au PC get same variance

})
