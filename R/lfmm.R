##' LFMM least-squares estimates with ridge penalty
##'
##' This function computes regularized least squares estimates 
##' for the parameters of latent factor mixed models using a ridge penalty.
##'
##' The algorithm minimizes the following penalized least-squares criterion
##' \deqn{ Lridge(U, V, B) = \frac{1}{2} ||Y - U V^{T} - X B^T||_{F}^2 
##' + \frac{\lambda}{2} ||B||^{2}_{2} ,}
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
##' data(example.data)
##' Y <- scale(example.data$genotype, scale = FALSE)
##' X <- scale(example.data$phenotype)
##' 
##' ## fits an lfmm model, i.e, computes B, U, V:
##' mod.lfmm <- lfmm_ridge(Y = Y, X = X, K = 6)
##' 
##' ## performs association testing using the fitted model:
##' pv <- lfmm_test(Y = Y, X = X, lfmm = mod.lfmm, calibrate = "gif")
##' 
##' ## Manhattan plot
##' plot(-log10(pv$calibrated.pvalue), pch = 19, cex = .2, col = "grey")
##' points(example.data$causal.set, 
##'       -log10(pv$calibrated.pvalue)[example.data$causal.set], 
##'        type = "h", col = "blue")
##'
##' 
##' ## Another example  
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

##' K-flod cross validation of LFMM least-squares estimates with ridge penalty
##'
##' This function split the dataset into a train set and test set to compute a
##' prediction error. The function \code{\link{lfmm_ridge}} is run with the
##' train set and the prediction error is evaluated with the test set.
##'
##'
##' @param Y a response variable matrix with n rows and p columns. 
##' Each column corresponds to a distinct response variable (e.g., SNP genotype, 
##' gene expression level, beta-normalized methylation profile, etc).
##' Response variables must be encoded as numeric.
##' @param X an explanatory variable matrix with n rows and d columns. 
##' Each column corresponds to a distinct explanatory variable (eg. phenotype).
##' Explanatory variables must be encoded as numeric.
##' @param Ks a list of integer for the number of latent factors in the regression model.
##' @param lambdas a list of numeric values for the regularization parameter.
##' @param n.fold.row number of folds along rows.
##' @param p.fold.col number of folds along columns.
##' @return a dataframe with prediction error for each value of lambda and K
##'
##' @export
##' @author cayek
##' @examples
##' library(ggplot2)
##' library(lfmm)
##'
##'  ## sample data
##'  K <- 3
##'  dat <- lfmm_sampler(n = 100, p = 1000, K = K,
##'                      outlier.prop = 0.1,
##'                      cs = c(0.8),
##'                      sigma = 0.2,
##'                      B.sd = 1.0,
##'                      U.sd = 1.0,
##'                      V.sd = 1.0)
##'
##'  ## run cross validation
##'  errs <- lfmm_ridge_CV(Y = dat$Y,
##'                          X = dat$X,
##'                          n.fold.row = 5,
##'                          n.fold.col = 5,
##'                          lambdas = c(1e-10, 1, 1e20),
##'                          Ks = c(1,2,3,4,5,6))
##'
##'  ## plot error
##'  ggplot(errs, aes(y = err, x = as.factor(K))) +
##'    geom_boxplot() +
##'    facet_grid(lambda ~ ., scale = "free")
##'
##'  ggplot(errs, aes(y = err, x = as.factor(lambda))) +
##'    geom_boxplot() +
##'    facet_grid(K ~ ., scales = "free")
##'
##' @seealso \code{\link{lfmm_ridge}}
lfmm_ridge_CV <- function(Y, X, n.fold.row, n.fold.col, lambdas , Ks) {

  ## init
  lfmm <- lfmm::ridgeLFMM(K = NULL,
                          lambda = NULL)
  dat <- LfmmDat(Y = Y, X = X)

  ## run and return
  return(lfmm::lfmm_CV(m  = lfmm, dat = dat,
                       n.fold.row = n.fold.row,
                       n.fold.col = n.fold.col,
                       Ks = Ks,
                       lambdas = lambdas))
}


##'  LFMM least-squares estimates with lasso penalty
##'
##' This function computes regularized least squares estimates 
##' for the parameters of latent factor mixed models using a lasso penalty. 
##' 
##' The algorithm minimizes the following penalized least-squares criterion
##'
##' \deqn{ Llasso(U, V, B) =
##' \frac{1}{2} ||Y - U V^{T} - X B^T||_{F}^2 + \frac{\lambda}{2}
##' ||B||^{2}_{2} , }
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

##' Statistical tests with latent factor mixed models
##' 
##' 
##' This function returns significance values for the association between each column of the 
##' response matrix, Y, and the explanatory variables, X, including correction for unobserved confounders 
##' (latent factors). The test is based on an LFMM fitted with a ridge or lasso penalty.
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
lfmm_test <- function(Y, X, lfmm, calibrate = "gif") {

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


##' Indirect effect sizes for latent factor models
##' 
##' 
##' This function returns 'indirect' effect sizes for the regression of X (of dimension 1) on the matrix Y,
##' as usually computed in genome-wide association studies.
##'
##' @param Y a response variable matrix with n rows and p columns. 
##' Each column is a response variable (numeric).
##' @param X an explanatory variable with n rows and d = 1 column (numeric). 
##' @param object an object of class \code{lfmm} returned by the \link{lfmm_lasso} 
##' or \link{lfmm_ridge} function.
##' @return a vector of length p containing all effect sizes for the regression 
##' of X on the matrix Y 
##' @export
##' @author cayek, francoio
##' @examples
##' library(lfmm)
##'
##' ## Simulation of 1000 phenotypes (x)
##' ## Only the 10 last variables have significant effect sizes
##' u <- matrix(rnorm(300, sd = 1), nrow = 100, ncol = 2)
##' v <- matrix(rnorm(3000, sd = 2), nrow = 2, ncol = 1000)
##' b <- matrix(c(rep(0, 990), rep(1000, 10)))
##' x <- 1000 + y%*%b + rnorm(100, sd = 100)
##' 
##' ## Compute direct effect sizes using lfmm_ridge
##' mod <- lfmm_ridge(Y = scale(y, scale = F), X = scale(x, scale = F), K = 2)
##' 
##' ## Compute indirect effect sizes using lfmm_ridge estimates
##' b.estimates <- effect_size(scale(y, scale = F), scale(x, scale = F), mod)
##' 
##' ## plot the 20 last effect sizes (true values are 0 and 1000)
##' plot(b.estimates[981:1000] )
effect_size <- function(Y, X, object){
  if (ncol(X) > 1) stop("Indirect effect sizes are computed for 
                        a single variable (d=1).")
  reg.lm <- function(i){
    dat <- data.frame(Y[,i], object$U)
    lm(X ~ ., data = dat)$coefficients[2]
  } 
  p <- ncol(Y)
  effect.sizes <- sapply(1:p, FUN = reg.lm)
  return(effect.sizes)
}

##' Predict polygenic risk scores from latent factor models
##' 
##' 
##' This function polygenic risk scores from latent factor models. It uses the 
##' 'indirect' effect sizes for the regression of X (a single phenotype) on the matrix Y,
##' for predicting phenotypes fromnew data.
##'
##' @param Y a response variable matrix with n rows and p columns. 
##' Each column is a response variable (numeric).
##' @param X an explanatory variable with n rows and d = 1 column (numeric). 
##' @param object an object of class \code{lfmm} returned by the \link{lfmm_lasso} 
##' or \link{lfmm_ridge} function.
##' @param fdr.level a numeric value for the FDR level in the lfmm test used to define
##' candidate variables for predicting new phenotypes.
##' @param newdata a matrix with n rows and p' columns, and similar to Y, on which 
##' predictions of X will be based. If NULL, Y is used as newdata. 
##' @return a list with the following attributes:
##'       - prediction: a vector of length n containing the predicted values for X. If 
##'       newdata = NULL, the fitted values are returned.
##'       - candidates: vector of candidate columns of Y on which the predictions are 
##'       built.
##' @export
##' @author cayek, francoio
##' @examples
##' library(lfmm)
##'
##' ## Simulation of 1000 phenotypes (x)
##' ## Only the 10 last variables have significant effect sizes
##' u <- matrix(rnorm(300, sd = 1), nrow = 100, ncol = 2)
##' v <- matrix(rnorm(3000, sd = 2), nrow = 2, ncol = 1000)
##' b <- matrix(c(rep(0, 990), rep(1000, 10)))
##' x <- 1000 + y%*%b + rnorm(100, sd = 100)
##' 
##' ## Compute direct effect sizes using lfmm_ridge
##' ## Note that centering is important (scale = F).
##' 
##' mod <- lfmm_ridge(Y = scale(y, scale = F), X = scale(x, scale = F), K = 2)
##' x.pred <- predict_lfmm(Y = scale(y, scale = F), scale(x, scale = F), mod)$pred
##' 
##' compare simulated and predicted/fitted phenotypes
##' plot(x, x.pred)
##' abline(0,1)
##' abline(lm(x.pred~x), col = 2)
predict_lfmm <- function(Y, X, object, fdr.level = 0.1, newdata = NULL){
  b.values <- effect_size(Y, X, object) 
  pvalues <- lfmm_test(Y, X, object,calibrate = "gif")$calibrated.pvalue
  p = length(pvalues)
  w = which(sort(pvalues) < fdr.level * (1:p)/p)
  candidates = order(pvalues)[w]
  if (is.null(newdata)) {newdata <- Y}
  x.pred <- newdata[,candidates] %*% matrix(b.values[candidates])
  return( list(prediction = x.pred, candidates = candidates) )
}
