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
  svd.res <- dat$svd(k = K + 1, K, K)# compute only singular value
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

    m[names(res)] <- res

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
  n = nrow(Y)
  p = ncol(Y)

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
    res.svd <- compute_soft_SVD_R(Yux, gamma)
    U <- res$U
    V <- res$V

    ## impute NA
    Yux <- tcrossprod(X, B) + tcrossprod(U, V)
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
