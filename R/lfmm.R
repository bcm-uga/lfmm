##' LFMM least-squares estimates with ridge penalty
##'
##' This function computes regularized least squares estimates 
##' for the parameters of latent factor mixed models using a ridge penalty.
##'
##' The algorithm minimizes the following penalized least-squares criterion
##' \deqn{ Lridge(U, V, B) = \frac{1}{2} ||Y - U V^{T} - X B^T||_{F}^2 
##' + \frac{\lambda}{2} norm{B}^{2}_{2} ,}
##' where Y is a response data matrix, X contains all explanatory variables, 
##' U denotes the score matrix, V is the loading matrix, and B is the effect size matrix.
##' 
##' @param Y a response variable matrix with n rows and p columns. 
##' Each column corresponds to a distinct response variable (e.g., SNP genotype, 
##' gene expression level, beta-normalized methylation profile, etc).
##' Response variables must be encoded as numeric.
##' @param X an explanatory variable matrix with n rows and d columns. 
##' Each column corresponds to a distinct explanatory variable (eg. phenotype).
##' Explanatory variables must be encoded as numeric.
##' @param K an integer for the number of latent factors in the regression model.
##' @param lambda a numeric value for the regularization parameter.
##' @return an object of class \code{lfmm} with the following attributes: 
##'  - U the latent variable score matrix with dimensions n x K.
##'  - V the latent variable axis matrix with dimensions p x K.
##'  - B the effect size matrix with dimensions p x d.
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
##' ## fit an LFMM with K = 3 latent factors
##' lfmm.res <- lfmm_ridge(Y = dat$Y, X = dat$X, K = 3, lambda = 1e-5)
##'
##' ## plot the effect size matrix
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

##'  LFMM least-squares estimates with lasso penalty
##'
##' This function computes regularized least squares estimates 
##' for the parameters of latent factor mixed models using a lasso penalty. 
##' 
##' The algorithm minimizes the following penalized least-squares criterion
##'
##' \deqn{ Llasso(U, V, B) =
##' frac{1}{2} ||Y - U V^{T} - X B^T||_{F}^2 + \frac{\lambda}{2}
##' \norm{B}^{2}_{2} , }
##' where Y is a response data matrix, X contains all explanatory variables, 
##' U denotes the score matrix, V is the loading matrix, and B is the effect 
##' size matrix.
##'
##' @param Y a response variable matrix with n rows and p columns. 
##' Each column is a response variable (e.g., SNP genotype, 
##' gene expression level, beta-normalized methylation profile, etc).
##' Response variables must be encoded as numeric.
##' @param X an explanatory variable matrix with n rows and d columns. 
##' Each column corresponds to a distinct explanatory variable (eg. phenotype).
##' Explanatory variables must be encoded as numeric.
##' @param K an integer for the number of latent factors in the regression model.
##' @param nozero.prop a numeric value for the expected proportion of non-zero effect sizes.
##' @param lambda.num a numeric value for the number of 'lambda' values (obscure).
##' @param lambda.min.ratio (obscure parameter) a numeric value for the smallest `lambda` value,
##'  A fraction of `lambda.max`, the data derived entry value (i.e. the smallest value for
##'   which all coefficients are zero).
##' @param lambda (obscure parameter) Smallest value of `lambda`. A fraction of 'lambda.max',
##'   the (data derived) entry value (i.e. the smallest value for which all
##'   coefficients are zero).
##' @param it.max an integer value for the number of iterations of the algorithm.
##' @param relative.err.epsilon a numeric value for a relative convergence error. Determine 
##' whether the algorithm converges or not.
##' @return an object of class \code{lfmm} with the following attributes: 
##'  - U the latent variable score matrix with dimensions n x K.
##'  - V the latent variable axes matrix with dimensions p x K.
##'  - B the effect size matrix with dimensions p x d.
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

##' Statistical tests of association between a response matrix and explanatory variables
##' with correction for unobserved confounders (latent factors).
##' This function returns significance values for the association between each column of the 
##' response matrix, Y, and the explanatory variables, X. The test is based on LFMM fitted with 
##' ridge or lasso penalties.
##'
##'
##' @param Y a response variable matrix with n rows and p columns. 
##' Each column is a response variable (numeric).
##' @param X an explanatory variable matrix with n rows and d columns. 
##' Each column corresponds to an explanatory variable (numeric).
##' @param lfmm an object of class \code{lfmm} returned by the \link{lfmm_lasso} 
##' or \link{lfmm_ridge} function
##' @param calibrate a character string, "gif" or "median+MAD". If the "gif" option is set (default), 
##' significance values are calibrated by using the genomic control method. Genomic control 
##' uses a robust estimate of the variance of z-scores called "genomic inflation factor". 
##' If the "median+MAD" option is set, the pvalues are calibrated by computing the median and MAD of the zscores. If \code{NULL}, the 
##' pvalues are not calibrated.
##' @return a list with the following attributes:
##'  - B the effect size matrix with dimensions p x d.
##'  - score a p x d matrix which contains z-scores for each explanatory variable (column of X)
##'  - pvalue a p x d matrix which contains p-values for each explanatory variable
##'  - calibrated.pvalue a p x d matrix which contains calibrated p-values for each explanatory variable
##'  - gif a numeric value for the genomic inflation factor.
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
