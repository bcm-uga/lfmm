##' LFMM \eqn{L_2} regularized estimator
##'
##' This function compute the \eqnL{_2} regularized least squares estimator of LFMM.
##'
##' The algorithm optimize the loss function
##' \deqn{ Lridge(U, V, B) = \frac{1}{2} ||Y - U V^{T} - X B^T||_{F}^2 + \frac{\lambda}{2} norm{B}^{2}_{2}}.
##'
##' @param Y Explained variables matrix. Each column is a variable.
##' @param X Explaining variables matrix. Each column is a variable.
##' @param K Number of latent factor.
##' @param lambda The value of regularization term.
##' @return A lfmm class object with following attribute: 
##'  - U the latent variable score matrix.
##'  - V the latent variable axes matrix.
##'  - B the effect size matrix.
##'
##' @export
##' @author cayek
##' @examples
##' library(lfmm)
##'
##' ## sample data
##' K <- 3
##' dat <- lfmm_sampler(n = 100, p = 1000, K = K,
##'                     outlier.prop = 0.1,
##'                     cs = c(0.8),
##'                     sigma = 0.2,
##'                     B.sd = 1.0,
##'                     U.sd = 1.0,
##'                     V.sd = 1.0)
##' ## run lfmm
##' lfmm.res <- lfmm_ridge(Y = dat$Y, X = dat$X, K = 3, lambda = 1e-5)
##'
##' ## plot size effect matrix
##' id <- seq_along(lfmm.res$B)
##' cols <- c('red', 'green')[as.numeric(id %in% dat$outlier) + 1]
##' plot(id, lfmm.res$B, col = cols)
lfmm_ridge <- function(Y, X, K, lambda = 1e-5) {

  ## init
  m <- ridgeLFMM(K = K, lambda = lambda)
  dat <- LfmmDat(Y = Y, X = X)

  ## run
  m <- lfmm_fit(m, dat)

  ## return
  m
}

##' LFMM \eqn{L_1} regularized estimator
##'
##' This function compute the $L_1$ regularized least squares estimator of LFMM.
##' The algorithm optimize the loss function
##'
##' \deqn{ Llasso(U, V, B) =
##' frac{1}{2} ||Y - U V^{T} - X B^T||_{F}^2 + \frac{\lambda}{2}
##' \norm{B}^{2}_{2} }.
##'
##' @param Y Explained variables matrix. Each column is a variable.
##' @param X Explaining variables matrix. Each column is a variable.
##' @param K Number of latent factor.
##' @param nozero.prop The proportion of line of B expected to be different from
##'   zero.
##' @param lambda.num The number of 'lambda' values.
##' @param lambda.min.ratio Smallest value for `lambda`, as a fraction of
##'   `lambda.max`, the (data derived) entry value (i.e. the smallest value for
##'   which all coefficients are zero).
##' @param lambda Smallest value for `lambda`, as a fraction of 'lambda.max',
##'   the (data derived) entry value (i.e. the smallest value for which all
##'   coefficients are zero).
##' @param it.max The number of iteration of the algorithm.
##' @param relative.err.epsilon The relative error used to determine if the
##'   algorithm converged.
##' @return A lfmm class object with following attribute: - U the latent
##'   variable score matrix. - V the latent variable axes matrix. - B the effect
##'   size matrix.
##'
##' @export
##' @author cayek
##' @examples
##' library(lfmm)
##'
##' ## sample data
##' K <- 3
##' dat <- lfmm_sampler(n = 100, p = 1000, K = K,
##'                     outlier.prop = 0.1,
##'                     cs = c(0.6),
##'                     sigma = 0.2,
##'                     B.sd = 1.0,
##'                     U.sd = 1.0,
##'                     V.sd = 1.0)
##' ## run lfmm
##' lfmm.res <- lfmm_lasso(Y = dat$Y, X = dat$X, K = 3, nozero.prop= 0.2)
##'
##' ## plot size effect matrix
##' id <- seq_along(lfmm.res$B)
##' cols <- c('red', 'green')[as.numeric(id %in% dat$outlier) + 1]
##' plot(id, lfmm.res$B, col = cols)
lfmm_lasso <- function(Y, X, K,
                       nozero.prop = 0.1,
                       lambda.num = 100,
                       lambda.min.ratio = 0.01,
                       lambda = NULL,
                       it.max = 100, relative.err.epsilon = 1e-6) {
  ## init
  m <- lassoLFMM(K = K,
                 nozero.prop = nozero.prop,
                 lambda.num = lambda.num,
                 lambda.min.ratio = lambda.min.ratio,
                 lambda = lambda)
  dat <- LfmmDat(Y = Y, X = X)

  ## run
  m <- lfmm_fit(m, dat, it.max = it.max, relative.err.epsilon = relative.err.epsilon)

  ## return
  m
}

##' Hypothesis testing of association of Y with X with correction by latent factor.
##' This function compute the pvalue of the association test of each column of Y
##' with X. The hypothesis testing take into account latent variables computed
##' by lfmm_lasso or lfmm_ridge.
##'
##'
##' @param Y Explained variables matrix. Each column is a variable.
##' @param X Explaining variables matrix. Each column is a variable.
##' @param lfmm Object returned by lfmm_lasso or lfmm_ridge
##' @param calibrate If "gif", pvalue are calibrated by computing the genomic
##'   inflation factor of the zscore. If "median+MAD", pvalue are calibrated by
##'   computing the median and MAD of the zscore. If NULL, pvalue are not
##'   calibrated.
##' @return A list with pvalue, zscore and effect.
##'
##' @export
##' @author cayek
##' @examples
##' library(lfmm)
##'
##' K <- 3
##' dat <- lfmm_sampler(n = 100, p = 1000, K = K,
##'                     outlier.prop = 0.1,
##'                     cs = c(0.8),
##'                     sigma = 0.2,
##'                     B.sd = 1.0,
##'                     U.sd = 1.0,
##'                     V.sd = 1.0)
##' ## lfmm
##' lfmm.res <- lfmm_ridge(Y = dat$Y, X = dat$X, K = 3, lambda = 1e-5)
##'
##' ## hp
##' hp.res <- lfmm_test(Y = dat$Y, X = dat$X, lfmm = lfmm.res)
##'
##' ## plot score
##' id <- seq_along(hp.res$calibrated.score)
##' cols <- c('red', 'green')[as.numeric(id %in% dat$outlier) + 1]
##' plot(id, hp.res$calibrated.score2, col = cols)
##'
##' ## plot pvalue
##' id <- seq_along(hp.res$calibrated.pvalue)
##' cols <- c('red', 'green')[as.numeric(id %in% dat$outlier) + 1]
##' plot(id, -log10(hp.res$calibrated.pvalue), col = cols)
lfmm_test <- function(Y, X, lfmm, calibrate = c("gif", "median+MAD")) {

  ## init
  dat <- LfmmDat(Y = Y, X = X)

  ## hp
  X <- cbind(dat$X, lfmm$U)
  d <- ncol(dat$X)
  hp <- hypothesis_testing_lm(dat, X)
  hp$score <- hp$score[,1:d, drop = FALSE]
  hp$pvalue <- hp$pvalue[,1:d, drop = FALSE]
  hp$B <- hp$B[,1:d, drop = FALSE]

  ## calibrate
  if (is.null(calibrate)) {
    NULL ## nothing
  } else if (calibrate == "median+MAD") {
    hp$mad <- apply(hp$score, 2, mad)
    hp$median <- apply(hp$score, 2, median)
    hp$calibrated.score <- sweep(hp$score, 2, hp$median, FUN = "-")
    hp$calibrated.score <- sweep(hp$score, 2, hp$mad, FUN = "/")
    hp$calibrated.pvalue <- compute_pvalue_from_zscore(hp$calibrated.score)
  } else if (calibrate == "gif") {
    hp$gif <- compute_gif(hp$score)
    hp$calibrated.score2 <- sweep(hp$score ^ 2, 2, hp$gif, FUN = "/")
    hp$calibrated.pvalue <- compute_pvalue_from_zscore2(hp$calibrated.score2, df = 1)
  }

  hp
}
