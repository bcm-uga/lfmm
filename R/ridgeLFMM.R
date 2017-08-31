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

  ## main algorithm
  m <- ridgeLFMM_main(m, dat, P.list)
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
    args$dat$productYt(t(args$P) %*% x)
  }
  res.rspectra <- RSpectra::svds(A = Af,
                                 Atrans = Atransf,
                                 k = m$K,
                                 nu = m$K, nv = m$K,
                                 opts = list(tol = 10e-10),
                                 dim = c(n, p),
                                 args = list(P = P.list$sqrt.P, dat = dat))
  ## res.rspectra <- svd(P.list$sqrt.P %*% dat$Y, nu = m$K, nv = m$K) ## debug
  m$U <- res.rspectra$u %*% diag(res.rspectra$d[1:m$K], m$K, m$K)
  m$U <- P.list$sqrt.P.inv %*% m$U
  m$V <- res.rspectra$v

  ## B
  Af <- function(x) {
    t(dat$productYt(x)) - tcrossprod(crossprod(x, m$U), m$V)
  }
  m$B <- compute_B_ridge(Af, dat$X, m$lambda)
  ## m$B <- compute_B_ridge(dat$Y - tcrossprod(m$U, m$V), dat$X, m$lambda) ## debug
  m
}

ridgeLFMM_withNA <- function(m, dat, relative.err.min = 1e-6, it.max = 100) {

  ## NA and input by median
  dat$missing.ind <- which(is.na(dat$Y))
  dat$Y <- impute_median(dat$Y)

  ## compute of P
  P.list <- compute_P(X = dat$X, lambda = m$lambda)

  ## main loop
  err2 <- .Machine$double.xmax
  it <- 1
  repeat {
    ## main algorithm
    m <- ridgeLFMM_main(m, dat, P.list)

    dat$impute_lfmm(m$U, m$V, m$B)
    err2.new <- dat$err2_lfmm(m$U, m$V, m$B)
    if(it > it.max || (abs(err2 - err2.new) / err2) < relative.err.min) {
      break
    }
    err2 <- err2.new
    message("It = ", it, "/", it.max, ", err2 = ", err2)
    it <- it + 1
  }

  ## to avoid side effect
  dat$Y[dat$missing.ind] <- NA

  m
}

##' @export
lfmm_fit.ridgeLFMM <- function(m, dat, it.max = 100, relative.err.min = 1e-6) {
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
lfmm_fit_knowing_loadings.ridgeLFMM <- function(m, dat) {
  m$U <- (dat$Y -  tcrossprod(dat$X, m$B)) %*% m$V
  m
}

##' @export
lfmm_CV.ridgeLFMM <- function(m, dat, n.fold.row, n.fold.col, lambdas , Ks,
                                              col.prop = 1.0) {

  params <- base::expand.grid(list(lambda = lambdas, K = Ks))
  CV(m = m,
     dat = dat,
     n.fold.row = n.fold.row,
     n.fold.col = n.fold.col,
     params = params,
     col.prop = col.prop)

}

