##' LFMM $L_2$ regularized estimator
##'
##' This function compute the $L_2$ regularized least squares estimator of LFMM.
##' The algorithm optimize the loss function
##'
##' $$ Lridge(\U, \V, \B) =
##' \frac{1}{2} \norm{Y - U V^{T} - X B^T}_{F}^2 + \frac{\lambda}{2}
##' \norm{B}^{2}_{2} $$.
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
##' lfmm.res <- ridge_lfmm(Y = dat$Y, X = dat$X, K = 3, lambda = 1e-5)
##'
##' ## plot size effect matrix
##' id <- seq_along(lfmm.res$B)
##' cols <- c('red', 'green')[as.numeric(id %in% dat$outlier) + 1]
##' plot(id, lfmm.res$B, col = cols)
ridge_lfmm <- function(Y, X, K, lambda = 1e-5) {

  ## init
  m <- ridgeLFMM(K = K, lambda = lambda)
  dat <- LfmmDat(Y = Y, X = X)

  ## run
  m <- lfmm_fit(m, dat)

  ## return
  m
}

##' LFMM $L_1$ regularized estimator
##'
##' This function compute the $L_1$ regularized least squares estimator of LFMM.
##' The algorithm optimize the loss function
##'
##' $$ Llasso(\U, \V, \B) =
##' \frac{1}{2} \norm{Y - U V^{T} - X B^T}_{F}^2 + \frac{\lambda}{2}
##' \norm{B}^{2}_{2} $$.
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
##' lfmm.res <- lasso_lfmm(Y = dat$Y, X = dat$X, K = 3, nozero.prop= 0.2)
##'
##' ## plot size effect matrix
##' id <- seq_along(lfmm.res$B)
##' cols <- c('red', 'green')[as.numeric(id %in% dat$outlier) + 1]
##' plot(id, lfmm.res$B, col = cols)
lasso_lfmm <- function(Y, X, K,
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

##' This function compute the pvalue of the association test of each column of Y
##' with X. The hypothesis testing take into account latent variables computed
##' by lasso_lfmm or ridge_lfmm.
##'
##'
##' @title Hypothesis testing of association of Y with X with correction by latent factor.
##' @param Y Explained variables matrix. Each column is a variable.
##' @param X Explaining variables matrix. Each column is a variable.
##' @param lfmm Object returned by lasso_lfmm or ridge_lfmm
##' @param calibrate If TRUE pvalue are calibrated by computing mean and MAD of
##'   zscore.
##' @return A list with pvalue, zscore and effect.
##'
##' @export
##' @author cayek
##' @examples
hypothesis_test_lfmm <- function(Y, X, lfmm, calibrate = TRUE) {

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
  if (calibrate) {
    hp$mad <- mad(hp$score)
    hp$median <- mad(hp$score)
    hp$calibrated.score <- (hp$score - hp$median) / hp$mad
    hp$calibrated.pvalue <- compute_pvalue_from_zscore(hp$calibrated.score)
  }

  hp
}
