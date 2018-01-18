SimulatedLfmmDat.builder <- setRefClass("SimulatedLfmmDat", contains = "LfmmDat",
                                        fields = c("B", "Epsilon", "U", "V", "outlier"))



##' LFMM generative data sampler
##'
##' Simulate data from the latent factor model.
##'
##' `lfmm_sample()` sample a response matrix Y and a primary variable X such that
##'
##'          Y = U t(V) + X t(B) + Epsilon.
##'
##' U,V, B and Epsilon are simulated according to normal multivariate distributions.
##' Moreover U and X are such that `cor(U[,i], X) = cs[i]`.
##'
##' @return A list with simulated data.
##' @author kevin caye, olivier francois
##' @param n number of observations.
##' @param p number of response variables.
##' @param K number of latent variables (factors).
##' @param outlier.prop proportion of outlier.
##' @param cs correlation with between X and U.
##' @param sigma standard deviation of residual errors.
##' @param B.sd standard deviation for the effect size (B).
##' @param B.mean mean of B.
##' @param U.sd standard deviations for K factors.
##' @param V.sd standard deviations for loadings.
##' @export
##'
##' @examples
##'
##' dat <- lfmm_sampler(n = 100, 
##'                     p = 1000, 
##'                     K = 3,
##'                     outlier.prop = 0.1,
##'                     cs = c(0.8),
##'                     sigma = 0.2,
##'                     B.sd = 1.0, 
##'                     B.mean = 0.0,
##'                     U.sd = 1.0, 
##'                     V.sd = 1.0)
lfmm_sampler <- function(n, p, K,
                         outlier.prop,
                         cs,
                         sigma = 0.2,
                         B.sd = 1.0,
                         B.mean = 0.0,
                         U.sd = 1.0,
                         V.sd = 1.0)
{

  ## sample outlier
  outlier <- sample.int(p, outlier.prop * p)
  outlier.nb = length(outlier)

  ## test cs
  if (length(cs) < K) {
    message("length(cs) < K. Filling cs with zero")
    cs <- c(cs, rep(0, times = K - length(cs)))
  }

  ## simulate U and X
  theta2 <- sum((cs/U.sd)^2) + 0.01
    
  Sigma <- diag(x = U.sd^2, nrow = K, ncol = K)
  Sigma <- rbind(Sigma, matrix(cs, nrow = 1))
  Sigma <- cbind(Sigma, matrix(c(cs, theta2), ncol = 1))
  UX <- MASS::mvrnorm(n, mu = rep(0.0, K + 1), Sigma = Sigma)
  U <- UX[,1:K, drop = FALSE]
  X <- UX[,K + 1, drop = FALSE]

  ## simulate V
  V <- MASS::mvrnorm(p, mu = rep(0.0, K), Sigma = V.sd^2 * diag(K))


  ## simulate B
  B <- matrix(0, p, 1)
  B[outlier, 1] <- rnorm(outlier.nb, B.mean, B.sd)

  ## sample error
  Epsilon = MASS::mvrnorm(n, mu = rep(0.0, p), Sigma = sigma^2 * diag(p))

  ## syntheses
  Y = U %*% t(V) + X %*% t(B) + Epsilon

  SimulatedLfmmDat.builder(Y = Y,
                           X = X,
                           outlier = outlier,
                           U = U,
                           V = V,
                           B = B,
                           meta = list(),
                           missing.ind = c())
}
