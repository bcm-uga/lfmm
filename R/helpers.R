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
    solve((crossprod(X,X) + lambda * D), crossprod(X, A))
  } else if(is.function(A)) {
    D <- diag(1, ncol(dat$X), ncol(dat$X))
    solve((crossprod(dat$X, dat$X) + lambda * D), A(X))
  }
}
