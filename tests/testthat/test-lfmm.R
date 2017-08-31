library(testthat)
context("lfmm")

test_that("ridge_lfmm", {

  K <- 3
  dat <- lfmm_sampler(n = 100, p = 1000, K = K,
                      outlier.prop = 0.1,
                      cs = c(0.8),
                      sigma = 0.2,
                      B.sd = 1.0,
                      U.sd = 1.0,
                      V.sd = 1.0)

  lfmm.res <- ridge_lfmm(Y = dat$Y, X = dat$X, K = 3, lambda = 1e-5)

  skip("plot")
  id <- seq_along(lfmm.res$B)
  cols <- c('red', 'green')[as.numeric(id %in% dat$outlier) + 1]
  plot(id, lfmm.res$B, col = cols)

})

test_that("ridge_lasso", {

  K <- 3
  dat <- lfmm_sampler(n = 100, p = 1000, K = K,
                      outlier.prop = 0.1,
                      cs = c(0.6),
                      sigma = 0.2,
                      B.sd = 1.0,
                      U.sd = 1.0,
                      V.sd = 1.0)

  lfmm.res <- lasso_lfmm(Y = dat$Y, X = dat$X, K = 3, nozero.prop= 0.2)

  skip("plot")
  id <- seq_along(lfmm.res$B)
  cols <- c('red', 'green')[as.numeric(id %in% dat$outlier) + 1]
  plot(id, lfmm.res$B, col = cols)

})

test_that("hypothesis_test_lfmm", {
  
  K <- 3
  dat <- lfmm_sampler(n = 100, p = 1000, K = K,
                      outlier.prop = 0.1,
                      cs = c(0.8),
                      sigma = 0.2,
                      B.sd = 1.0,
                      U.sd = 1.0,
                      V.sd = 1.0)

  ## lfmm
  lfmm.res <- ridge_lfmm(Y = dat$Y, X = dat$X, K = 3, lambda = 1e-5)

  ## hp
  hp.res <- hypothesis_test_lfmm(Y = dat$Y, X = dat$X, lfmm = lfmm.res, calibrate = TRUE)

  skip("plot")

  ## plot score
  id <- seq_along(hp.res$calibrated.score)
  cols <- c('red', 'green')[as.numeric(id %in% dat$outlier) + 1]
  plot(id, hp.res$calibrated.score, col = cols)

  ## plot pvalue
  id <- seq_along(hp.res$calibrated.pvalue)
  cols <- c('red', 'green')[as.numeric(id %in% dat$outlier) + 1]
  plot(id, -log10(hp.res$calibrated.pvalue), col = cols)

})
