## ------------------------------------------------------------------------
#devtools::install_github("bcm-uga/lfmm")

## ------------------------------------------------------------------------
library(lfmm)

## ------------------------------------------------------------------------
## Simulated phenotypes for Arabidopsis thaliana SNP data
data("example.data")
## Simulated (and real) methylation levels for sun-exposed tissue sampled
data("skin.exposure")

## ------------------------------------------------------------------------
Y <- example.data$genotype
pc <- prcomp(Y)
plot(pc$sdev[1:20]^2, xlab = 'PC', ylab = "Variance explained")
points(6,pc$sdev[6]^2, type = "h", lwd = 3, col = "blue")

## ------------------------------------------------------------------------
Y <- skin.exposure$beta.value
pc <- prcomp(Y)
plot(pc$sdev[1:20]^2, xlab = 'PC', ylab = "Variance explained")
points(2,pc$sdev[2]^2, type = "h", lwd = 3, col = "blue")

## ------------------------------------------------------------------------
 Y <- example.data$genotype
 X <- example.data$phenotype #scaled phenotype

## ------------------------------------------------------------------------
## Fit an LFMM, i.e, compute B, U, V estimates
 mod.lfmm <- lfmm_ridge(Y = Y, 
                        X = X, 
                        K = 6)

## ------------------------------------------------------------------------
 ## performs association testing using the fitted model:
 pv <- lfmm_test(Y = Y, 
                 X = X, 
                 lfmm = mod.lfmm, 
                 calibrate = "gif")

## ------------------------------------------------------------------------
pvalues <- pv$calibrated.pvalue 
qqplot(rexp(length(pvalues), rate = log(10)),
       -log10(pvalues), xlab = "Expected quantile",
       pch = 19, cex = .4)
abline(0,1)

## ------------------------------------------------------------------------
 ## Manhattan plot
 plot(-log10(pvalues), 
      pch = 19, 
      cex = .2, 
      xlab = "SNP", ylab = "-Log P",
      col = "grey")
 points(example.data$causal.set, 
       -log10(pvalues)[example.data$causal.set], 
        type = "h", 
        col = "blue")

## ------------------------------------------------------------------------
 Y <- scale(skin.exposure$beta.value)
 X <- scale(as.numeric(skin.exposure$exposure))

## ------------------------------------------------------------------------
 ## Fit and LFMM, i.e, compute B, U, V estimates
 mod.lfmm <- lfmm_ridge(Y = Y, 
                        X = X, 
                        K = 2)
 
 ## Perform association testing using the fitted model:
 pv <- lfmm_test(Y = Y, 
                 X = X, 
                 lfmm = mod.lfmm, 
                 calibrate = "gif")

## ------------------------------------------------------------------------
 ## Manhattan plot
 plot(-log10(pv$calibrated.pvalue), 
      pch = 19, 
      cex = .3,
      xlab = "Probe", ylab = "-Log P",
      col = "grey")
 causal.set <- seq(11, 1496, by = 80)
 points(causal.set, 
       -log10(pv$calibrated.pvalue)[causal.set], 
        col = "blue")

## ------------------------------------------------------------------------
 Y <- example.data$genotype
 X <- example.data$phenotype #scaled phenotype

## ---- message = FALSE----------------------------------------------------
## Fit an LFMM, i.e, compute B, U, V estimates
 mod.lfmm <- lfmm_lasso(Y = Y, 
                        X = X, 
                        K = 6,
                        nozero.prop = 0.01)

## ------------------------------------------------------------------------
 ## performs association testing using the fitted model:
 pv <- lfmm_test(Y = Y, 
                 X = X, 
                 lfmm = mod.lfmm, 
                 calibrate = "gif")

## ------------------------------------------------------------------------
pvalues <- pv$calibrated.pvalue 
qqplot(rexp(length(pvalues), rate = log(10)),
       -log10(pvalues), xlab = "Expected quantile",
       pch = 19, cex = .4)
abline(0,1)

## ------------------------------------------------------------------------
 ## Manhattan plot
 plot(-log10(pvalues), 
      pch = 19, 
      cex = .2, 
      xlab = "SNP", ylab = "-Log P",
      col = "grey")
 points(example.data$causal.set, 
       -log10(pvalues)[example.data$causal.set], 
        type = "h", 
        col = "blue")

## ------------------------------------------------------------------------
 ## Simulation of 1000 genotypes for 100 individuals (y)
 u <- matrix(rnorm(300, sd = 1), nrow = 100, ncol = 2) 
 v <- matrix(rnorm(3000, sd = 2), nrow = 2, ncol = 1000)
 y <- matrix(rbinom(100000, size = 2, 
                   prob = 1/(1 + exp(-0.3*(u%*%v 
                   + rnorm(100000, sd = 2))))),
                   nrow = 100,
                   ncol = 1000)

 ## Simulation of 1000 phenotypes (x)
 ## Only the last 10 genotypes have significant effect sizes (b)
 b <- matrix(c(rep(0, 990), rep(6000, 10)))
 x <- y%*%b + rnorm(100, sd = 100)

## ------------------------------------------------------------------------
 mod <- lfmm_ridge(Y = y, 
                   X = x,
                   K = 2)

## ------------------------------------------------------------------------
  candidates <- 991:1000 #causal loci
  b.values <- effect_size(Y = y, X = x, lfmm.object = mod) 
  x.pred <- scale(y[,candidates], scale = F)%*% matrix(b.values[candidates])

## ------------------------------------------------------------------------
 ##Compare simulated and predicted/fitted phenotypes
 plot(x - mean(x), x.pred, 
      pch = 19, col = "grey", 
      xlab = "Observed phenotypes (centered)", 
      ylab = "Predicted from PRS")
 abline(0,1)
 abline(lm(x.pred ~ scale(x, scale = FALSE)), col = 2)

## ------------------------------------------------------------------------
 pred <- predict_lfmm(Y = y, 
                      X = x,
                      fdr.level = 0.25, 
                      mod)
 
 ##Compare simulated and predicted/fitted phenotypes
 plot(x - mean(x), pred$pred, 
      pch = 19, col = "grey", 
      xlab = "Observed phenotypes (centered)", 
      ylab = "Predicted from PRS")
 abline(0,1)
 abline(lm(pred$pred ~ scale(x, scale = FALSE)), col = 2)

