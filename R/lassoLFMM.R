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

  ## compute gamma
  res <- list()
  Af <- function(x, args) {
    dat$productY(x)
  }
  Atransf <- function(x, args) {
    dat$productYt(x)
  }
  svd.res <- compute_svd(Af, Atransf, K + 1, K, K, dim = c(nrow(dat$Y), ncol(dat$Y)))
  res$gamma <- (svd.res$d[K] + svd.res$d[K + 1]) / 2
  U <-svd.res$u[,1:K] %*% diag(svd.res$d[1:K], K, K)
  V <- svd.res$v[,1:K]

  ## compute B
  Af <- function(x) {
    t(dat$productYt(x)) - tcrossprod(crossprod(x, U), V)
  }
  B <- compute_B_ridge(Af, dat$X, 0.0)

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

lassoLFMM_main <- function(m, dat, it.max = 100, relative.err.epsilon = 1e-6) {

  m <- lassoLFMM_init(m, dat)

  ## NA and input by median
  dat$missing.ind <- which(is.na(dat$Y))
  dat$Y <- impute_median(dat$Y)

  ## main loop
  for (lambda in m$params$lambda.range) {

    m <- lassoLFMM_loop(m, dat,
                        m$params$gamma, lambda,
                        relative.err.epsilon,
                        it.max)

    nozero.prop <- mean(m$B > 0.0)
    message("=== lambda = ", lambda, ", no zero B proportion = ", nozero.prop)
    if( nozero.prop > m$nozero.prop) {
      break
    }
  }

  ## to avoid side effect
  dat$Y[dat$missing.ind] <- NA

  m
}

##' @export
MatrixFactorizationR_fit.lassoLFMM <- function(m, dat, it.max = 100, relative.err.epsilon = 1e-6) {

  res <- lassoLFMM_main(m, dat, it.max, relative.err.epsilon)
  res
}

lassoLFMM_loop <- function(m, dat, gamma, lambda, relative_err_epsilon, it_max) {

  ## constants
  n = nrow(dat$Y)
  p = ncol(dat$Y)

  ## variables
  err = 0.0
  err_new = dat$err2_lfmm(m$U, m$V, m$B)
  relative_err = .Machine$double.xmax
  it = 1

  ## main loop
  while ((it <= it_max) && (relative_err > relative_err_epsilon)) {
    err = err_new;
    message("It = ", it , "/", it_max, ", err2 = " ,err)

    ## step B
    Af <- function(x) {
      t(dat$productYt(x)) - tcrossprod(crossprod(x, m$U), m$V)
    }
    m$B <- compute_B_lasso(Af, dat$X, lambda)

    ## compute W = UV^T
    Af <- function(x, args) {
      dat$productY(x)- dat$X %*% crossprod(m$B, x)
    }
    Atransf <- function(x, args) {
      dat$productYt(x) - m$B %*% crossprod(dat$X, x)
    }
    res.rspectra <- compute_svd_soft(Af, Atransf, gamma, m$K, dim = c(nrow(dat$Y), ncol(dat$Y)))
    m$U <- res.rspectra$u %*% diag(res.rspectra$d, length(res.rspectra$d), length(res.rspectra$d))
    m$V <- res.rspectra$v

    ## impute NA
    dat$impute_lfmm(m$U, m$V, m$B)

    ## err
    err_new = dat$err2_lfmm(m$U, m$V, m$B)
    relative_err = abs(err_new - err) / err
    it = it + 1
  }
  m
}
