library(testthat)
context("lfmm")

test_that("lfmm_ridge", {

  K <- 3
  dat <- lfmm_sampler(n = 100, p = 1000, K = K,
                      outlier.prop = 0.1,
                      cs = c(0.8),
                      sigma = 0.2,
                      B.sd = 1.0,
                      U.sd = 1.0,
                      V.sd = 1.0)

  lfmm.res <- lfmm_ridge(Y = dat$Y, X = dat$X, K = 3, lambda = 1e-5)

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

  lfmm.res <- lfmm_lasso(Y = dat$Y, X = dat$X, K = 3, nozero.prop= 0.2)

  skip("plot")
  id <- seq_along(lfmm.res$B)
  cols <- c('red', 'green')[as.numeric(id %in% dat$outlier) + 1]
  plot(id, lfmm.res$B, col = cols)

})

test_that("lfmm_test", {

  K <- 3
  dat <- lfmm_sampler(n = 100, p = 1000, K = K,
                      outlier.prop = 0.01,
                      cs = c(0.8),
                      sigma = 0.2,
                      B.sd = 1.0,
                      U.sd = 1.0,
                      V.sd = 1.0)

  ## intercept
  dat$X <- cbind(1, dat$X)

  ## lfmm
  lfmm.res <- lfmm_ridge(Y = dat$Y, X = dat$X, K = 3, lambda = 1e-5)

  ## hp with gif
  hp.res.gif <- lfmm_test(Y = dat$Y, X = dat$X, lfmm = lfmm.res, calibrate = "gif")
  expect_equal(length(hp.res.gif$gif), 2)

  ## hp with median+MAD
  hp.res.mad <- lfmm_test(Y = dat$Y, X = dat$X, lfmm = lfmm.res, calibrate = "median+MAD")
  expect_equal(length(hp.res.mad$mad), 2)
  expect_equal(length(hp.res.mad$median), 2)

  ## hp with median+MAD
  hp.res <- lfmm_test(Y = dat$Y, X = dat$X, lfmm = lfmm.res, calibrate = NULL)

  skip("plot")

  ## with gif
  ## plot score
  d <- 2
  id <- seq_along(hp.res.gif$calibrated.score[,d])
  cols <- c('red', 'green')[as.numeric(id %in% dat$outlier) + 1]
  plot(id, hp.res.gif$calibrated.score2[,d], col = cols)
  ## plot pvalue
  cols <- c('red', 'green')[as.numeric(id %in% dat$outlier) + 1]
  plot(id, -log10(hp.res.gif$calibrated.pvalue[,d]), col = cols)
  hist(hp.res.gif$calibrated.pvalue[,d])
  hist(hp.res.gif$pvalue[,d])

  ## median+mad
  d <- 1
  id <- seq_along(hp.res.mad$calibrated.score[,d])
  cols <- c('red', 'green')[as.numeric(id %in% dat$outlier) + 1]
  plot(id, hp.res.mad$calibrated.score[,d], col = cols)
  ## plot pvalue
  cols <- c('red', 'green')[as.numeric(id %in% dat$outlier) + 1]
  plot(id, -log10(hp.res.mad$calibrated.pvalue[,d]), col = cols)
  hist(hp.res.mad$calibrated.pvalue[,d])
  hist(hp.res.mad$pvalue[,d])

})
