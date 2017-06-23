library(testthat)
context("Sampler")


test_that("lfmm_sampler", {

  dat <- lfmm_sampler(n = 100, p = 1000, K = 3,
                      outlier.prop = 0.1,
                      cs = c(0.8),
                      sigma = 0.2,
                      B.sd = 1.0,
                      U.sd = 1.0,
                      V.sd = 1.0)

  expect_equal(dim(dat$Y), c(100, 1000))
  expect_equal(dim(dat$X), c(100, 1))
  expect_equal(dim(dat$B), c(1000, 1))
  expect_equal(dim(dat$U), c(100, 3))
  expect_equal(dim(dat$V), c(1000, 3))
  expect_equal(dim(dat$Epsilon), c(100, 1000))

})

