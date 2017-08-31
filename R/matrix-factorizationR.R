#' R package with matrix factorization algorithms
#'
#'
#' @docType package
#'
#' @name lfmm
#' @importFrom Rcpp evalCpp
#' @importFrom foreach foreach %:% %do% %dopar%
#' @useDynLib lfmm
#' @import RcppEigen
NULL

#' Fit the model
#'
#' @export
lfmm_fit <- function(m, dat, ...) {
  UseMethod("lfmm_fit")
}

#' Fit the model when latent factor loadings are known
#'
#' @export
lfmm_fit_knowing_loadings <- function(m, dat, ...) {
  UseMethod("lfmm_fit_knowing_loadings")
}

#' Cross validation
#'
#' @export
lfmm_CV <- function(m, dat, n.fold.row, n.fold.col, ...) {
  UseMethod("lfmm_CV")
}

#' Impute Y with a fitted model.
#'
#' @export
lfmm_impute <- function(m, dat, ...) {
  UseMethod("lfmm_impute")
}

#' Compute the residual error
#'
#' @export
lfmm_residual_error2 <- function(m, dat, ...) {
  UseMethod("lfmm_residual_error2")
}
