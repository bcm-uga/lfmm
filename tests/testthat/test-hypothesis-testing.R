library(testthat)
context("hypothesis testing")


test_that("hypothesis_testing_lm", {

  dat <- lfmm_sampler(n = 100, p = 1000, K = 1,
                      outlier.prop = 0.1,
                      cs = c(0.8),
                      sigma = 0.2,
                      B.sd = 1.0,
                      U.sd = 1.0,
                      V.sd = 1.0)

  X <- cbind(dat$X, dat$U, rnorm(100))
  hp <- hypothesis_testing_lm(dat, X = X)

  lm.res <- lm(dat$Y ~ X - 1)
  B.lm <- lm.res$coefficients
  expect_lt(mean(abs(t(B.lm) - hp$B)), 1e-15)

  s <- summary(lm.res)

  ## E
  E.lm <- lm.res$residuals


  ## score
  score <- sapply(seq_along(s), function(i) s[[i]]$coefficients[,3])
  dim(score)
  expect_lt(mean(abs(t(score) - hp$score)), 1e-10)

  ## pvalue
  pvalue <- sapply(seq_along(s), function(i) s[[i]]$coefficients[,4])
  dim(pvalue)
  expect_lt(mean(abs(t(pvalue) - hp$pvalue)), 1e-14)
  hist(pvalue[3,]) ## ok 
  hist(hp$pvalue[,3]) ## ok

})

