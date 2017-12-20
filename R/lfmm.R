##' LFMM least-squares estimates with ridge penalty
##'
##' This function computes regularized least squares estimates 
##' for latent factor mixed models using a ridge penalty.
##'
##' The algorithm minimizes the following penalized least-squares criterion
##' \deqn{ L(U, V, B) = \frac{1}{2} ||Y - U V^{T} - X B^T||_{F}^2 
##' + \frac{\lambda}{2} ||B||^{2}_{2} ,}
##' where Y is a response data matrix, X contains all explanatory variables, 
##' U denotes the score matrix, V is the loading matrix, B is the (direct) effect 
##' size matrix, and lambda is a regularization parameter.
##'
##' @param Y a response variable matrix with n rows and p columns. 
##' Each column corresponds to a distinct response variable (e.g., SNP genotype, 
##' gene expression level, beta-normalized methylation profile, etc).
##' Response variables must be encoded as numeric.
##' @param X an explanatory variable matrix with n rows and d columns. 
##' Each column corresponds to a distinct explanatory variable (eg. phenotype).
##' Explanatory variables must be encoded as numeric variables.
##' @param K an integer for the number of latent factors in the regression model.
##' @param lambda a numeric value for the regularization parameter.
##' @param algorithm analytical solution algorithm or numerical alternated algorithm.
##' The analytical solution algorithm is based on global minimum of the loss function and
##' the computation is very quick. The numerical alternated algorithm converges toward a local
##' minimum of the loss function. The interest of the alternated algorithm is that it does not
##' imply the computation of a \eqn{n \times n} matrix. Whereas the analytical algorithm requires
##' the computation of a \eqn{n \times n} matrix, which can be a problem when n is very large.
##' @param it.max an integer value for the number of iterations of the algorithm. Only for the
##' alternated numerical algorithm
##' @param relative.err.epsilon a numeric value for a relative convergence error. Determine
##' whether the algorithm converges or not. Only for the
##' alternated numerical algorithm.
##' @return an object of class \code{lfmm} with the following attributes:
##'  - U the latent variable score matrix with dimensions n x K,
##'  - V the latent variable axis matrix with dimensions p x K,
##'  - B the effect size matrix with dimensions p x d.
##' @details The response variable matrix Y and the explanatory variable are centered.
##' @export
##' @author cayek, francoio
##' @examples
##' 
##' library(lfmm)
##' 
##' ## a GWAS example with Y = SNPs and X = phenotype
##' data(example.data)
##' Y <- example.data$genotype
##' X <- example.data$phenotype
##' 
##' ## Fit an LFMM with K = 6 factors
##' mod.lfmm <- lfmm_ridge(Y = Y, 
##'                        X = X, 
##'                        K = 6)
##' 
##' ## Perform association testing using the fitted model:
##' pv <- lfmm_test(Y = Y, 
##'                 X = X, 
##'                 lfmm = mod.lfmm, 
##'                 calibrate = "gif")
##' 
##' ## Manhattan plot with causal loci shown
##' 
##' pvalues <- pv$calibrated.pvalue
##' plot(-log10(pvalues), pch = 19, 
##'      cex = .2, col = "grey", xlab = "SNP")
##' points(example.data$causal.set, 
##'       -log10(pvalues)[example.data$causal.set], 
##'        type = "h", col = "blue")
##'
##' 
##' ## An EWAS example with Y = methylation data and X = exposure
##' Y <- scale(skin.exposure$beta.value)
##' X <- scale(as.numeric(skin.exposure$exposure))
##' 
##' ## Fit an LFMM with 2 latent factors
##' mod.lfmm <- lfmm_ridge(Y = Y,
##'                        X = X, 
##'                        K = 2)
##'                        
##' ## Perform association testing using the fitted model:
##' pv <- lfmm_test(Y = Y, 
##'                 X = X,
##'                 lfmm = mod.lfmm, 
##'                 calibrate = "gif")
##'                 
##' ## Manhattan plot with true associations shown
##' pvalues <- pv$calibrated.pvalue
##' plot(-log10(pvalues), 
##'      pch = 19, 
##'      cex = .3,
##'      xlab = "Probe",
##'      col = "grey")
##'      
##' causal.set <- seq(11, 1496, by = 80)
##' points(causal.set, 
##'       -log10(pvalues)[causal.set], 
##'        col = "blue")
lfmm_ridge <- function(Y, X, K, lambda = 1e-5, algorithm = c("analytical", "alternated"),
                       it.max = 100, relative.err.min = 1e-6) {

  ## init
  m <- ridgeLFMM(K = K, lambda = lambda)
  m$algorithm <- algorithm[1]
  dat <- LfmmDat(Y = scale(Y, scale = FALSE), X = scale(X, scale = FALSE))

  ## run
  m <- lfmm_fit(m, dat, it.max = it.max, relative.err.min = relative.err.min)

  ## return
  m
}

##' Cross validation of LFMM estimates with ridge penalty
##'
##' This function splits the data set into a train set and a test set, and returns 
##' a prediction error. The function \code{\link{lfmm_ridge}} is run with the
##' train set and the prediction error is evaluated from the test set.
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
##' @param n.fold.row number of cross-validation folds along rows.
##' @param p.fold.col number of cross-validation folds along columns.
##' @return a dataframe containing prediction errors for all values of lambda and K
##' @details The response variable matrix Y and the explanatory variable are centered.
##'
##' @export
##' @author cayek, francoio
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
lfmm_ridge_CV <- function(Y, X, n.fold.row, n.fold.col, lambdas, Ks) {

  ## init
  lfmm <- lfmm::ridgeLFMM(K = NULL,
                          lambda = NULL)
  dat <- LfmmDat(Y = scale(Y, scale = FALSE), X = scale(X, scale = FALSE))

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
##' for latent factor mixed models using a lasso penalty. 
##' 
##' The algorithm minimizes the following penalized least-squares criterion
##'
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
##'  - U the latent variable score matrix with dimensions n x K,
##'  - V the latent variable axes matrix with dimensions p x K,
##'  - B the effect size matrix with dimensions p x d.
##' @details The response variable matrix Y and the explanatory variable are centered.
##'
##' @export
##' @author cayek, francoio
##' @examples
##' 
##' library(lfmm)
##' 
##' ## a GWAS example with Y = SNPs and X = phenotype
##' data(example.data)
##' Y <- example.data$genotype
##' X <- example.data$phenotype
##' 
##' ## Fit an LFMM with 6 factors
##' mod.lfmm <- lfmm_lasso(Y = Y, 
##'                        X = X, 
##'                        K = 6,
##'                        nozero.prop = 0.01)
##'                        
##' ## Perform association testing using the fitted model:
##' pv <- lfmm_test(Y = Y, 
##'                 X = X, 
##'                 lfmm = mod.lfmm, 
##'                 calibrate = "gif")
##'                 
##' ## Manhattan plot with causal loci shown
##' pvalues <- pv$calibrated.pvalue
##' plot(-log10(pvalues), 
##'       pch = 19, cex = .2, 
##'       col = "grey", xlab = "SNP")
##'       
##' points(example.data$causal.set, 
##'        -log10(pvalues)[example.data$causal.set], 
##'        type = "h", col = "blue")
##'        
##' ## An EWAS example with Y = methylation data 
##' ## and X = exposure
##' Y <- scale(skin.exposure$beta.value)
##' X <- scale(as.numeric(skin.exposure$exposure))
##' 
##' ## Fit an LFMM with 2 latent factors
##' mod.lfmm <- lfmm_lasso(Y = Y,
##'                        X = X, 
##'                        K = 2,
##'                        nozero.prop = 0.01)
##'                        
##' ## Perform association testing using the fitted model:
##' pv <- lfmm_test(Y = Y, 
##'                 X = X,
##'                 lfmm = mod.lfmm, 
##'                 calibrate = "gif")
##'                 
##' ## Manhattan plot with true associations shown
##' pvalues <- pv$calibrated.pvalue
##' plot(-log10(pvalues), 
##'      pch = 19, 
##'      cex = .3,
##'      xlab = "Probe",
##'      col = "grey")
##'      
##' causal.set <- seq(11, 1496, by = 80)
##' points(causal.set, 
##'        -log10(pvalues)[causal.set], 
##'        col = "blue")
lfmm_lasso <- function(Y, X, K, 
                       nozero.prop = 0.01,
                       lambda.num = 100,
                       lambda.min.ratio = 0.01,
                       lambda = NULL,
                       it.max = 100, 
                       relative.err.epsilon = 1e-6) {
  ## init
  m <- lassoLFMM(K = K,
                 nozero.prop = nozero.prop,
                 lambda.num = lambda.num,
                 lambda.min.ratio = lambda.min.ratio,
                 lambda = lambda)
  dat <- LfmmDat(Y = scale(Y, scale = FALSE), X = scale(X, scale = FALSE))

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
##'  - score a p x d matrix which contains z-scores for each explanatory variable (columns of X),
##'  - pvalue a p x d matrix which contains p-values for each explanatory variable,
##'  - calibrated.pvalue a p x d matrix which contains calibrated p-values for each explanatory variable,
##'  - gif a numeric value for the genomic inflation factor.
##' @details The response variable matrix Y and the explanatory variable are centered.
##' @export
##' @author cayek, francoio
##' @examples
##' 
##' library(lfmm)
##' 
##' ## a GWAS example with Y = SNPs and X = phenotype
##' data(example.data)
##' Y <- example.data$genotype
##' X <- example.data$phenotype
##' 
##' ## Fit an LFMM with K = 6 factors
##' mod.lfmm <- lfmm_ridge(Y = Y, 
##'                        X = X, 
##'                        K = 6)
##' 
##' ## Perform association testing using the fitted model:
##' pv <- lfmm_test(Y = Y, 
##'                 X = X, 
##'                 lfmm = mod.lfmm, 
##'                 calibrate = "gif")
##' 
##' ## Manhattan plot with causal loci shown
##' 
##' pvalues <- pv$calibrated.pvalue
##' plot(-log10(pvalues), pch = 19, 
##'      cex = .2, col = "grey", xlab = "SNP")
##' points(example.data$causal.set, 
##'       -log10(pvalues)[example.data$causal.set], 
##'        type = "h", col = "blue")
##'
##' 
##' ## An EWAS example with Y = methylation data and X = exposure
##' Y <- scale(skin.exposure$beta.value)
##' X <- scale(as.numeric(skin.exposure$exposure))
##' 
##' ## Fit an LFMM with 2 latent factors
##' mod.lfmm <- lfmm_ridge(Y = Y,
##'                        X = X, 
##'                        K = 2)
##'                        
##' ## Perform association testing using the fitted model:
##' pv <- lfmm_test(Y = Y, 
##'                 X = X,
##'                 lfmm = mod.lfmm, 
##'                 calibrate = "gif")
##'                 
##' ## Manhattan plot with true associations shown
##' pvalues <- pv$calibrated.pvalue
##' plot(-log10(pvalues), 
##'      pch = 19, 
##'      cex = .3,
##'      xlab = "Probe",
##'      col = "grey")
##'      
##' causal.set <- seq(11, 1496, by = 80)
##' points(causal.set, 
##'       -log10(pvalues)[causal.set], 
##'        col = "blue")
lfmm_test <- function(Y, X, lfmm, calibrate = "gif") {

  ## init
  dat <- LfmmDat(Y = scale(Y, scale = FALSE), X = scale(X, scale = FALSE)) ## WARNING : scale here !!

  ## hp
  X <- cbind(dat$X, lfmm$U)
  d <- ncol(dat$X)
  if (class(lfmm) == "ridgeLFMM") {
    lfmm.la <- lfmm$lambda } else {
    lfmm.la <- 0.0   
    }
  hp <- hypothesis_testing_lm(dat, X, lfmm.la)
  hp$score <- hp$score[,1:d, drop = FALSE]
  hp$pvalue <- hp$pvalue[,1:d, drop = FALSE]
  hp$B <- hp$B[ ,1:d, drop = FALSE]

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


##' Direct effect sizes estimated from latent factor models
##' 
##' 
##' This function returns 'direct' effect sizes for the regression of X (of dimension 1) on the matrix Y,
##' as usually computed in genome-wide association studies.
##'
##' @param Y a response variable matrix with n rows and p columns. 
##' Each column is a response variable (numeric).
##' @param X an explanatory variable with n rows and d = 1 column (numeric). 
##' @param lfmm.object an object of class \code{lfmm} returned by the \link{lfmm_lasso} 
##' or \link{lfmm_ridge} function.
##' @return a vector of length p containing all effect sizes for the regression 
##' of X on the matrix Y 
##' @details The response variable matrix Y and the explanatory variable are centered.
##' @export
##' @author cayek, francoio
##' @examples
##' library(lfmm)
##'
##' ## Simulation of 1000 genotypes for 100 individuals (y)
##' u <- matrix(rnorm(300, sd = 1), nrow = 100, ncol = 2) 
##' v <- matrix(rnorm(3000, sd = 2), nrow = 2, ncol = 1000)
##' y <- matrix(rbinom(100000, size = 2, 
##'                   prob = 1/(1 + exp(-0.3*(u%*%v 
##'                   + rnorm(100000, sd = 2))))),
##'                   nrow = 100,
##'                   ncol = 1000)
##'
##' #PCA of genotypes, 3 main axes of variation (K = 2) 
##' plot(prcomp(y))
##'   
##' ## Simulation of 1000 phenotypes (x)
##' ## Only the last 10 genotypes have significant effect sizes (b)
##' b <- matrix(c(rep(0, 990), rep(6000, 10)))
##' x <- y%*%b + rnorm(100, sd = 100)
##' 
##' ## Compute effect sizes using lfmm_ridge
##' ## Note that centering is important (scale = F).
##' mod.lfmm <- lfmm_ridge(Y = y, 
##'                   X = x,
##'                   K = 2)
##'               
##' ## Compute direct effect sizes using lfmm_ridge estimates
##' b.estimates <- effect_size(y, x, mod.lfmm)
##' 
##' ## plot the last 30 effect sizes (true values are 0 and 6000)
##' plot(b.estimates[971:1000])
##' abline(0, 0)
##' abline(6000, 0, col = 2)
##' 
##' ## Prediction of phenotypes
##' candidates <- 991:1000 #set of causal loci 
##' x.pred <- scale(y[,candidates], scale = F) %*% matrix(b.estimates[candidates])
##' 
##' ## Check predictions
##' plot(x - mean(x), x.pred, 
##'      pch = 19, col = "grey", 
##'      xlab = "Observed phenotypes (centered)", 
##'      ylab = "Predicted from PRS")
##'      abline(0,1)
##'      abline(lm(x.pred ~ scale(x, scale = FALSE)), col = 2)
effect_size <- function(Y, X, lfmm.object){
  if (ncol(X) > 1) stop("Indirect effect sizes are computed for 
                        a single variable (d=1).")
  Y = scale(Y, scale = FALSE)
  X = scale(X, scale = FALSE)
  reg.lm <- function(i){
    datafr <- data.frame(Y[,i], lfmm.object$U)
    lm(X ~ ., data = datafr)$coefficients[2]
  } 
  p <- ncol(Y)
  effect.sizes <- sapply(1:p, FUN = reg.lm)
  return(effect.sizes)
}

##' Predict polygenic scores from latent factor models
##' 
##' 
##' This function computes polygenic risk scores from the estimates of latent factor models. 
##' It uses the indirect' effect sizes for the regression of X (a single phenotype) on the matrix Y,
##' for predicting phenotypic values for new genotype data.
##'
##' @param Y a response variable matrix with n rows and p columns, typically containing genotypes. 
##' Each column is a response variable (numeric).
##' @param X an explanatory variable with n rows and d = 1 column (numeric) representing a phenotype 
##' with zero mean across the sample. 
##' @param lfmm.object an object of class \code{lfmm} returned by the \link{lfmm_lasso} 
##' or \link{lfmm_ridge} function, computed for X and Y.
##' @param fdr.level a numeric value for the FDR level in the lfmm test used to define
##' candidate variables for predicting new phenotypes.
##' @param newdata a matrix with n rows and p' columns, and similar to Y, on which 
##' predictions of X will be based. If NULL, Y is used as new data. 
##' @return a list with the following attributes:
##' - prediction: a vector of length n containing the predicted values for X. If newdata 
##' = NULL, the fitted values are returned.
##' - candidates: a vector of candidate columns of Y on which the predictions are built.
##' @details The response variable matrix Y and the explanatory variable are centered.
##' @export
##' @author cayek, francoio
##' @examples
##' library(lfmm)
##'
##' ## Simulation of 1000 genotypes for 100 individuals (y)
##' u <- matrix(rnorm(300, sd = 1), nrow = 100, ncol = 2) 
##' v <- matrix(rnorm(3000, sd = 2), nrow = 2, ncol = 1000)
##' y <- matrix(rbinom(100000, size = 2, 
##'                   prob = 1/(1 + exp(-0.3*(u%*%v 
##'                   + rnorm(100000, sd = 2))))),
##'                   nrow = 100,
##'                   ncol = 1000)
##'
##' #PCA of genotypes, 2 main axes of variation (K = 2) 
##' plot(prcomp(y))
##'   
##' ## Simulation of 1000 phenotypes (x)
##' ## Only the last 10 genotypes have significant effect sizes (b)
##' b <- matrix(c(rep(0, 990), rep(6000, 10)))
##' x <- y%*%b + rnorm(100, sd = 100)
##' 
##' ## Compute effect sizes using lfmm_ridge
##' mod <- lfmm_ridge(Y = y, 
##'                   X = x,
##'                   K = 2)
##'               
##' x.pred <- predict_lfmm(Y = y, 
##'                        X = x,
##'                        fdr.level = 0.25, 
##'                        mod)
##'                     
##' x.pred$candidates
##' 
##' ##Compare simulated and predicted/fitted phenotypes
##' plot(x - mean(x), x.pred$pred, 
##'      pch = 19, col = "grey", 
##'      xlab = "Observed phenotypes (centered)", 
##'      ylab = "Predicted from PRS")
##' abline(0,1)
##' abline(lm(x.pred$pred ~ scale(x, scale = FALSE)), col = 2)
predict_lfmm <- function(Y, X, lfmm.object, fdr.level = 0.1, newdata = NULL){
  Y = scale(Y, scale = FALSE)
  X = scale(X, scale = FALSE)
  b.values <- effect_size(Y, X, lfmm.object) 
  pvalues <- lfmm_test(Y, X, lfmm.object, calibrate = "gif")$calibrated.pvalue
  pvalues[is.na(pvalues)] <- 1
  p = length(pvalues)
  w = which(sort(pvalues) < fdr.level * (1:p)/p)
  candidates = order(pvalues)[w]
  if (is.null(newdata)) {newdata <- Y}
  x.pred <- newdata[,candidates] %*% matrix(b.values[candidates])
  return( list(prediction = x.pred, candidates = candidates) )
}


##' Forward inclusion tests with latent factor mixed models
##' 
##' 
##' This function tests for association between each column of the response matrix, Y, 
##' and the explanatory variables, X, by recursively conditioning on the top hits in the set 
##' of explanatory variables. The conditional tests are based on LFMMs with ridge penalty.
##'
##'
##' @param Y a response variable matrix with n rows and p columns. 
##' Each column is a response variable (numeric).
##' @param X an explanatory variable matrix with n rows and d = 1 column (eg. phenotype). 
##' @param K an integer for the number of latent factors in the regression model.
##' @param lambda a numeric value for the regularization parameter.
##' @param niter an integer value for the number of forward inclusion tests.
##' @param scale a boolean value, \code{TRUE} if the explanatory variable, X, is scaled 
##' (recommended option). 
##' @param rev.confounder a boolean value. If \code{TRUE} confounders are revaluated in each 
##' conditional test. May take some time (default = \code{TRUE}. 
##' @param candidate.list a vector of integers corresponding to response variables (columns in Y), 
##' which are known candidates for association. If \code{NULL}, a list of candidates 
##' is built in during the algorithm run.
##' 
##' @return a list with the following attributes: 
##'  - candidates a vector of niter response variables (column labels in Y) detected as top hits in 
##'  each conditional association analysis. 
##'  - log.p a vector of uncorrected log p-values for checking that the algorithm behaves well (but not trustable for testing). 
##' @details The response variable matrix Y and the explanatory variable are centered.
##' 
##' @export
##' @author cayek, francoio
##' @examples
##' library(lfmm)
##' data("example.data")
##' Y <- example.data$genotype
##' X <- example.data$phenotype #scaled variable
##' 
##' ## fits an LFMM, i.e, computes B, U, V:
##' mod.lfmm <- lfmm_ridge(Y = Y,
##'                        X = X, 
##'                        K = 6)
##'                        
##' ## performs initial association testing using the fitted model:
##' pv <- lfmm_test(Y = Y, 
##'                 X = X,
##'                 lfmm = mod.lfmm,
##'                 calibrate = "gif")
##' ## Manhattan plot 
##' plot(-log10(pv$calibrated.pvalue), 
##'       pch = 19, 
##'       cex = .2,
##'       col = "grey")
##'       
##' ## Start forward tests (3 iterations)       
##' obj <- forward_test(Y, 
##'                     X, 
##'                     K = 6, 
##'                     niter = 3, 
##'                     scale = TRUE)
##' 
##' ## Record Log p.values for the 3 top hits
##' log.p <-  obj$log.p
##' log.p
##' 
##' ## Check perfect hits for each causal SNPs (labelled from 1 to 20)
##' obj$candidate %in% example.data$causal.set
##' 
##' ## Check for candidates at distance 20 SNPs (about 10kb)
##' theta <- 20
##' ## Number of hits for each causal SNPs (1-20)
##'  hit.3 <- as.numeric(
##'           apply(sapply(obj$candidate, 
##'           function(x) abs(x - example.data$causal.set) < theta), 
##'           2, 
##'           which))
##' ## Number of hits for each causal SNPs (1-20) 
##' table(hit.3)
##' 
##' 
##' ## Continue forward tests (2 additional iterations)       
##' obj <- forward_test(Y, 
##'                     X, 
##'                     K = 6, 
##'                     niter = 2,
##'                     candidate.list = obj$candidates,
##'                     scale = TRUE)
##' 
##' ## Record Log p.values for all 5 top hits
##' log.p <-  c(log.p, obj$log.p)
##' log.p
##' 
##' ## Check perfect hits for each causal SNPs (labelled from 1 to 20)
##' obj$candidate %in% example.data$causal.set
##' 
##' ## Check for candidates at distance 5 SNPs (about 2.5kb)
##' theta <- 5
##' ## Number of hits for each causal SNPs (1-20)
##'  hit.5 <- as.numeric(
##'           apply(sapply(obj$candidate, 
##'           function(x) abs(x - example.data$causal.set) < theta), 
##'           2, 
##'           which))
##' ## Number of hits for each causal SNPs (1-20)          
##' table(hit.5)
##' 
##' ## Plot log P
##' plot(log.p, xlab = "Conditional test iteration", ylab="Top hit log(p)")
forward_test <- function(Y, 
                         X, 
                         K, 
                         niter = 5, 
                         scale = FALSE,
                         candidate.list = NULL,
                         rev.confounder = FALSE,
                         lambda = 1e-5){
  
  if (ncol(X) > 1) stop("The function works with 
                          a single explanatory variable 
                          (d = 1).")
  
  if (!is.null(candidate.list)){
      if (!is.integer(candidate.list)) 
      {
         stop("Error in candidate.list: A list 
                of integers is required.")}
      if (max(candidate.list)>ncol(Y))
      {
        stop("Error in candidate.list: 
             max(candidate.list)>ncol(Y).")}
    }
  
  candidates <- candidate.list 
  log.pvalues <- NULL
  
  for (i in 1:niter){
    if (scale){
      X.m <- scale(cbind(X, Y[,candidates]))} 
    else {
        X.m <- cbind(X, Y[,candidates])  
      }
    
    if (i == 1 | rev.confounder == TRUE){ 
      mod.lfmm <- lfmm_ridge(Y = Y, 
                           X = X.m, 
                           K = K, 
                           lambda = lambda) 
      }
    pv <- lfmm_test(Y = Y, 
                    X = X.m, 
                    lfmm = mod.lfmm, 
                    calibrate = "gif")
    j <- 1
    while (order(-log10(pv$calibrated.pvalue[,1]), 
                 decreasing = TRUE)[j] %in% candidates) {
      j <- j + 1
    }
    
    candidates <- c(candidates,
                    order(-log10(pv$calibrated.pvalue[,1]),
                    decreasing = TRUE)[j])
    log.pvalues <- c(log.pvalues, 
                     sort(-log10(pv$calibrated.pvalue[,1]), 
                     decreasing = TRUE)[j])
    }
  return(list(candidates = candidates, log.p = log.pvalues))
}
