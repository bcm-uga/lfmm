##' Hypothesis testing with lm
##'
##' linear model: 
##' Y = X B^T + E
##' 
##' @author cayek
##' @export
hypothesis_testing_lm <- function(dat, X) {

  d <- ncol(X)
  p <- ncol(dat$Y)
  effective.degree.freedom <- nrow(dat$Y) - ncol(X)

  res <- list()
  ## B
  Af <- function(x) {
    t(dat$productYt(x))
  }
  res$B <- compute_B_ridge(Af, X, 0.0)

  ## compute Var(E)
  res$epsilon.sigma2 <- dat$sigma2_lm(X, res$B, effective.degree.freedom)

  ## compute Var(B)
  aux <- solve(crossprod(X))
  res$B.sigma2 <- t(matrix(diag(aux), d, 1) %*% matrix(res$epsilon.sigma2, 1, p))

  ## compute zscore
  res$score <- res$B / sqrt(res$B.sigma2)

  ## compute pvalue
  res$pvalue <- compute_pvalue_from_tscore(res$score, df = effective.degree.freedom)

  res
}
