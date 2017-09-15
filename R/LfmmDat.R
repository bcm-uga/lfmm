LfmmDat.builder <- setRefClass("LfmmDat", contains = "Dat",
                               fields = c("X"),
                               methods = list(
                                 predict_lfmm_knowing_loadings =
                                   function(V, B, unknown.j) {
                                     n <- nrow(.self$Y)
                                     ## Compute U
                                     U <- (.self$Y[,-unknown.j,drop = FALSE] -
                                           tcrossprod(.self$X, B[-unknown.j,, drop = FALSE])) %*%
                                       V[-unknown.j, ,drop = FALSE]
                                     ## predict Y
                                     tcrossprod(U, V[unknown.j,,drop = FALSE]) +
                                     tcrossprod(.self$X, B[unknown.j,,drop = FALSE])
                                   },
                                 impute_lfmm = function(U, V, B) {
                                   impute_lfmm_cpp(.self$Y, .self$X, U, V, B, .self$missing.ind)
                                 },
                                 err2_lfmm = function(U, V, B) {
                                   err2_lfmm_cpp(.self$Y, .self$X, U, V, B)
                                 },
                                 err2s_lfmm = function(U, V, B) {
                                   err2s_lfmm_cpp(.self$Y, .self$X, U, V, B)
                                 },
                                 sigma2_lm = function(X, B, nb.df) {
                                   if (is.matrix(.self$Y) && is.double(.self$Y)) {
                                     res <- sum2_lm_cpp(.self$Y, X, B) / nb.df
                                   } else {
                                     res <- 1:ncol(.self$Y)
                                     aux.f <- function(j) {
                                       aux <- .self$Y[,j] - tcrossprod(X , B[j,,drop = FALSE])
                                       sum(aux * aux)
                                     }
                                     res <- sapply(res,aux.f)
                                     res <- res / nb.df
                                     res
                                   }
                                 }
                               )
                               )

#' Class which store data
#'
#'
#' @export
LfmmDat <- function(Y, X, missing = TRUE) {
  dat <- LfmmDat.builder(Y = read_input(Y),
                         X = read_input(X),
                         missing.ind = NULL,
                         meta = list())
  if (missing) {
    dat$missing.ind <- which(is.na(dat$Y))
  }
  dat
}

#' Class which store data
#'
#'
#' @export
SimulatedLfmmDat <- function(Y, X, outlier, U, V, B) {
  dat <- SimulatedLfmmDat.builder(Y = read_input(Y),
                                  meta = list(),
                                  X = read_input(X),
                                  outlier = read_input(outlier),
                                  U = read_input(U),
                                  B = read_input(B),
                                  V = read_input(V))
  dat$missing.ind <- which(is.na(dat$Y))
  dat
}

