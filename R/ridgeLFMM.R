##' ridge LFMM
##'
##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @author cayek
##' @export
ridgeLFMM <- function(K, lambda) {
  m <- list( K = K,
            lambda = lambda)
  class(m) <- "ridgeLFMM"
  m
}

ridgeLFMM_noNA<- function(m, dat) {
  ## compute of P
  P.list <- compute_P(X = dat$X, lambda = m$lambda)

  ## c++ algorithm
  res <- ridgeLFMM_main(Y = dat$Y,
                        X = dat$X,
                        lambda = m$lambda,
                        K = m$K,
                        sqrtP = P.list$sqrt.P,
                        invSqrtP = P.list$sqrt.P.inv)
  m[names(res)] <- res
  m
}

ridgeLFMM_main <- function(m, dat, P.list) {

  d <- ncol(dat$X)
  n <- nrow(dat$Y)
  p <- ncol(dat$Y)

  ## UV
  Af <- function(x, args) {
    args$P %*% args$dat$productY(x)
  }
  Atransf <- function(x, args) {
    args$P %*% args$dat$productYt(x)
  }
  res.rspectra <- RSpectra::svds(A = dat$productY,
                                 Atrans = dat$productYt,
                                 k = m$K,
                                 nu = m$K, nv = m$K,
                                 opts = list(tol = 10e-10),
                                 dim = c(n, p),
                                 args = list(P = P.list$sqrt.P, dat = dat))
  m$U <- res.rspectra$u %*% diag(res.rspectra)
  m$U <- p.list$sqrt.P.inv %*% m$U
  m$V <- res.rspectra$v

  ## B
  Af <- function(x) {
    t(dat$productYt(x)) - tcrossprod(crossprod(x, m$U), m$V)
  }
  m$B <- compute_B_ridge(Af, dat$X, lambda)
  m
}

ridgeLFMM_withNA<- function(m, dat, relative.err.min = 1e-6, it.max = 100) {

  ## NA and input by median
  missing.index <- which(is.na(dat$Y))
  dat$Y <- impute_median(dat$Y)

  ## compute of P
  P.list <- compute_P(X = dat$X, lambda = m$lambda)

  ## main loop
  err2 <- .Machine$double.xmax
  it <- 1
  repeat{
    ## c++ algorithm
    res <- ridgeLFMM_main_R(Y = dat$Y,
                            X = dat$X,
                            lambda = m$lambda,
                            K = m$K,
                            sqrtP = P.list$sqrt.P,
                            invSqrtP = P.list$sqrt.P.inv)

    dat$Y[missing.index] <- NA
    dat <- MatrixFactorizationR_impute.ridgeLFMM(res, dat)
    err2.new <- MatrixFactorizationR_residual_error2.ridgeLFMM(res, dat)

    if(it > it.max || (abs(err2 - err2.new) / err2) < relative.err.min) {
      break
    }
    err2 <- err2.new
    message("It = ", it, "/", it.max, ", err2 = ", err2)
    it <- it + 1
  }
  m[names(res)] <- res
  m
}

##' @export
MatrixFactorizationR_fit.ridgeLFMM <- function(m, dat, it.max = 100, relative.err.min = 1e-6) {
  ## test if there missing value in Y
  if (anyNA(dat$Y)) {
    res <- ridgeLFMM_withNA(m, dat,
                     relative.err.min = relative.err.min,
                     it.max = it.max)
  } else {
    res <- ridgeLFMM_noNA(m, dat)
  }
}

##' Fit assuming V and B
##'
##' @export
MatrixFactorizationR_fit_knowing_loadings.ridgeLFMM <- function(m, dat) {
  m$U <- (dat$Y -  tcrossprod(dat$X, m$B)) %*% m$V
  m
}

##' @export
MatrixFactorizationR_impute.ridgeLFMM <- function(m, dat) {
  missing.index <- which(is.na(dat$Y))
  dat$Y[missing.index] <- tcrossprod(m$U, m$V)[missing.index]
  dat$Y[missing.index] <- dat$Y[missing.index] + tcrossprod(dat$X, m$B)[missing.index]
  dat
}

MatrixFactorizationR_CV.ridgeLFMM <- function(m, dat, kfold.row, kfold.col, lambdas , Ks) {

  params <- base::expand.grid(list(lambda = lambdas, K = Ks))
  CV(m = m,
     dat = dat,
     kfold.row = kfold.row,
     kfold.col = kfold.col,
     params = params)
 }

##' @export
MatrixFactorizationR_residual_error2.ridgeLFMM <- function(m, dat) {
  E <- dat$Y - tcrossprod(m$U, m$V)
  E <- E - tcrossprod(dat$X, m$B)
  mean(E ^ 2)
}
