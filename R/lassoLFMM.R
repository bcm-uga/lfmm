##' Lasso LFMM
##'
##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @author cayek
##' @export
lassoLFMM <- function(K, nozero.prop = 0.1,
                      lambda.K = 100, lambda.eps = 0.001) {
  m <- list(K = K,
            nozero.prop = nozero.prop,
            lambda.K = 100,
            lambda.eps = 0.001)
  class(m) <- "lassoLFMM"
  m
}


lassoLFMM_heuristic_gamma_lambda_range<- function(m, dat) {

  K <- m$K

  res <- list()
  svd.res <- svd(dat$Y, nu = K, nv = K) # compute only singular value
  res$gamma <- (svd.res$d[K] + svd.res$d[K + 1]) / 2


  Y <- dat$Y - tcrossprod(svd.res$u[,1:K] %*% diag(svd.res$d[1:K], K, K),
                          svd.res$v[,1:K])

  B <- compute_B_ridge(Y, dat$X, 0.0)

  ## lambda max and min
  lambda.max <- max(B)
  # lambda.min = lambda.eps * lambda.max like in Friedman et al. 2010
  lambda.min <- lambda.max * m$lambda.eps

  ## strategie presented in Friedman et al. 2010
  ##   log scaled sequence
  res$lambda.range <- exp(seq(log(lambda.max), log(lambda.min), length.out = m$lambda.K))

  res
}

lassoLFMM_init <- function(m, dat) {

  ## compute gamma
  m$params <- lassoLFMM_heuristic_gamma_lambda_range(m, dat)

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
  m
}

lassoLFMM_noNA <- function(m, dat, it.max = 100, relative.err.epsilon = 1e-6) {

  m <- lassoLFMM_init(m, dat)

  ## main loop
  for (lambda in m$params$lambda.range) {
    ## c++
    res <- lassoLFMM_main(dat$Y, dat$X,
                          m$params$gamma, lambda,
                          relative.err.epsilon,
                          it.max,
                          m$U,
                          m$V,
                          m$B)

    m[names(res)] <- res

    nozero.prop <- mean(m$B > 0.0)
    message("=== lambda = ", lambda, ", no zero B proportion = ", nozero.prop)
    if( nozero.prop > m$nozero.prop) {
      break
    }
  }
  m
}

lassoLFMM_withNA<- function(m, dat, it.max = 100, relative.err.epsilon = 1e-6) {

  ## NA and input by median
  missing.index <- which(is.na(dat$Y))
  dat$Y <- impute_median(dat$Y)

  ## lasso init
  m <- lassoLFMM_init(m, dat)

  ## main loop
  for (lambda in m$params$lambda.range) {
    ## c++
    res <- lassoLFMM_main_R(dat$Y, dat$X,
                            m$params$gamma, lambda,
                            relative.err.epsilon,
                            it.max,
                            m$U,
                            m$V,
                            m$B,
                            missing.index)

    m[names(res)] <- res

    nozero.prop <- mean(m$B > 0.0)
    message("=== lambda = ", lambda, ", no zero B proportion = ", nozero.prop)
    if( nozero.prop > m$nozero.prop) {
      break
    }
  }

  m
}


##' @export
MatrixFactorizationR_fit.lassoLFMM <- function(m, dat, it.max = 100, relative.err.epsilon = 1e-6) {

  ## test if there missing value in Y
  if (anyNA(dat$Y)) {
    res <- lassoLFMM_withNA(m, dat, it.max, relative.err.epsilon)
  } else {
    res <- lassoLFMM_noNA(m, dat, it.max, relative.err.epsilon)
  }
}

compute_soft_SVD_R <- function(X, gamma) {
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

lassoLFMM_main_R <- function(Y, X, gamma, lambda, relative_err_epsilon, it_max,
                             U0, V0, B0,
                             missing.index = NULL) {

  ## constants
  n = nrow(Y)
  p = ncol(Y)

  ## variables
  U = U0
  V = V0
  Yux = Y
  B = B0
  err = 0.0
  Yux = Y - tcrossprod(U, V)
  Yux = Yux - tcrossprod(X, B)
  err_new = mean(Yux ^ 2)
  relative_err = .Machine$double.xmax
  it = 1

  ## main loop
  while ((it <= it_max) && (relative_err > relative_err_epsilon)) {
    err = err_new;
    message("It = ", it , "/", it_max, ", err2 = " ,err)

    ## step B
    Yux = Y - tcrossprod(U, V)
    B = compute_B_lasso(Yux, X, lambda)

    ## compute W = UV^T
    Yux = Y - tcrossprod(X, B)
    res <- compute_soft_SVD_R(Yux, gamma)
    U <- res$U
    V <- res$V

    ## impute NA
    Yux = tcrossprod(X, B) + tcrossprod(U, V)
    if (!is.null(missing.index)) {
      Y[missing.index] <- Yux[missing.index]
    }

    ## err
    Yux = Y - Yux
    err_new = mean(Yux ^ 2)
    relative_err = abs(err_new - err) / err
    it = it + 1
  }

  list(U = U,
       V = V,
       B = B)
}
