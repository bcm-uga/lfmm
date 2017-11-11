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


test_that("lfmm_ridge alternated", {

  ## rmk : the alternated algorithm is not guarenty to converge to the global minimum.
  K <- 3
  dat <- lfmm_sampler(n = 100, p = 1000, K = K,
                      outlier.prop = 0.1,
                      cs = c(0.8),
                      sigma = 0.2,
                      B.sd = 1.0,
                      U.sd = 1.0,
                      V.sd = 1.0)
  dat$X <- cbind(dat$X, matrix(rnorm(100), 100,1))

  lfmm.res <- lfmm_ridge(Y = dat$Y, X = dat$X, K = 3, lambda = 1e-5)
  lfmm.alt.res <- lfmm_ridge(Y = dat$Y, X = dat$X, K = 3,
                             lambda = 1e-5,
                             algorithm = "alternated",
                             it.max=100,
                             relative.err.min=1e-8)

  expect_lte(max((lfmm.res$B - lfmm.alt.res$B)^2), 1e-1)
  ## axes not in same order.
  ## cor(lfmm.res$U, lfmm.alt.res$U)

  ## hypothesis testing
  test.res <- lfmm_test(Y = dat$Y, X = dat$X, lfmm = lfmm.res)
  test.altr.res <- lfmm_test(Y = dat$Y, X = dat$X, lfmm = lfmm.alt.res)

  expect_lte(mean((test.altr.res$B - test.res$B)^2), 1e-2)
})

test_that("lfmm_ridge CV", {

  library(ggplot2)

  ## sample data
  K <- 3
  dat <- lfmm_sampler(n = 100, p = 1000, K = K,
                      outlier.prop = 0.1,
                      cs = c(0.8),
                      sigma = 0.2,
                      B.sd = 1.0,
                      U.sd = 1.0,
                      V.sd = 1.0)

  ## run cross validation
  errs <- lfmm_ridge_CV(Y = dat$Y,
                          X = dat$X,
                          n.fold.row = 5,
                          n.fold.col = 5,
                          lambdas = c(1e-10, 1, 1e20),
                          Ks = c(1,2,3,4,5,6))

  ## plot error
  skip("plot")
  ggplot(errs, aes(y = err, x = as.factor(K))) +
    geom_boxplot() +
    facet_grid(lambda ~ ., scale = "free")

  ggplot(errs, aes(y = err, x = as.factor(lambda))) +
    geom_boxplot() +
    facet_grid(K ~ ., scales = "free")

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
                      outlier.prop = 0.05,
                      cs = c(0.8),
                      sigma = 0.2,
                      B.sd = 1.0,
                      U.sd = 1.0,
                      V.sd = 1.0)

  ## random co variate
  dat$X <- cbind(dat$X, matrix(rnorm(100, 100, 1)))

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
