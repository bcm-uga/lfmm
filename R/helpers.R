##' @export
impute_median <- function(Y) {
  n <- nrow(Y)
  p <- ncol(Y)
  missing.index <- which(is.na(Y))
  meds <- apply(Y, 2, median, na.rm = TRUE)
  Y[missing.index] <- (matrix(1,n,1) %*% matrix(meds,1,p))[missing.index] ## :-(
  Y
}

compute_B_ridge <- function(A, X, lambda) {
  D <- diag(1, ncol(X), ncol(X))
  if (is.matrix(A)) {
    t(solve((crossprod(X,X) + lambda * D), crossprod(X, A)))
  } else if(is.function(A)) {
    D <- diag(1, ncol(X), ncol(X))
    t(solve((crossprod(X, X) + lambda * D), A(X)))
  }
}

compute_B_lasso <- function(A, X, lambda) {
  B_hat <- compute_B_ridge(A, X, 0.0)
  sign(B_hat) * sapply((abs(B_hat) - lambda), function(x) max(0,x))
}

#' score are assume to follow student distibution with df degre of freedom
compute_pvalue_from_tscore <- function(score, df) {
  apply(score, 1:2, function(z) 2 * pt(abs(z), df = df, lower.tail = FALSE))
}

#' score are assume to follow normal distibution
compute_pvalue_from_zscore <- function(score, mean = 0, sd = 1) {
  apply(score, 1:2, function(z) 2 * pnorm(abs(z), mean = mean, sd = sd, lower.tail = FALSE))
}
