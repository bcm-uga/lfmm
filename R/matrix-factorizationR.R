#' R package with matrix factorization algorithms
#'
#'
#' @docType package
#'
#' @name MatrixFactorizationR
#' @importFrom Rcpp evalCpp
#' @useDynLib MatrixFactorizationR
#' @import RcppEigen
NULL

#' Fit the model
#'
#' @export
MatrixFactorizationR_fit <- function(m, dat, ...) {
  UseMethod("MatrixFactorizationR_fit")
}

#' Fit the model when latent factor loadings are known
#'
#' @export
MatrixFactorizationR_fit_knowing_loadings <- function(m, dat, ...) {
  UseMethod("MatrixFactorizationR_fit_knowing_loadings")
}

#' Cross validation
#'
#' @export
MatrixFactorizationR_CV <- function(m, dat, n.fold.row, n.fold.col, ...) {
  UseMethod("MatrixFactorizationR_CV")
}

#' Impute Y with a fitted model.
#'
#' @export
MatrixFactorizationR_impute <- function(m, dat, ...) {
  UseMethod("MatrixFactorizationR_impute")
}

#' Compute the residual error
#'
#' @export
MatrixFactorizationR_residual_error2 <- function(m, dat, ...) {
  UseMethod("MatrixFactorizationR_residual_error2")
}
