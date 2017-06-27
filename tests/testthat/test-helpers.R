library(testthat)
context("Helpers")


test_that("compute_P", {

  set.seed(135435)
  dat <- lfmm_sampler(n = 100, p = 1000, K = 3,
                      outlier.prop = 0.1,
                      cs = c(0.8),
                      sigma = 0.2,
                      B.sd = 1.0,
                      U.sd = 1.0,
                      V.sd = 1.0)

  res <- compute_P(X = dat$X,
                   lambda = 1e-5)

  expect_lte(mean(abs(diag(1, nrow(dat$X), nrow(dat$X))) - res$sqrt.P %*% res$sqrt.P.inv), 1e-16)
  ## mean(abs(res$P - res$sqrt.P %*% res$sqrt.P))
  
})

test_that("compute_B_ridge", {

  dat <- lfmm_sampler(n = 100, p = 1000, K = 3,
                      outlier.prop = 0.1,
                      cs = c(0.8),
                      sigma = 0.2,
                      B.sd = 1.0,
                      U.sd = 1.0,
                      V.sd = 1.0)

  X <- cbind(matrix(1,100,1), dat$X)
  B <- compute_B_ridge(A = dat$Y,
                       X = X,
                       lambda = 1e-3)
  expect_equal(dim(B), c(1000, 2))
  ## hist(B[,1]) ## RMK when lambda -> infinity B -> 0 OK

})

test_that("compute_B_lasso", {

  skip("deprecated")
  dat <- lfmm_sampler(n = 100, p = 1000, K = 3,
                      outlier.prop = 0.1,
                      cs = c(0.8),
                      sigma = 0.2,
                      B.sd = 1.0,
                      U.sd = 1.0,
                      V.sd = 1.0)

  B_lasso <- function(A, X, lambda) {
    B_hat <- solve((crossprod(X,X) ), crossprod(X, A))
    sign(B_hat) * ((abs(B_hat) - lambda) %>% purrr::map_dbl(~ max(.x, 0)))
  }

  lambda <- 1.5e0

  X <- cbind(matrix(1,100,1), dat$X)
  X <- svd(X)$u ## to have othogonal value
  B.cpp <- compute_B_lasso(Y = dat$Y,
                       X = X,
                       lambda = lambda)

  expect_equal(dim(B.cpp), c(1000, 2))
  ## hist(B[,2]) 
  ## mean(B[,2] == 0)

  ## with R
  B.r <- B_lasso(dat$Y, X, lambda)

  ## comp
  expect_lt(mean(abs(t(B.r) - B.cpp)), 1e-10)
  expect_equal(mean(B.r != 0), mean(B.cpp != 0))
  
})
