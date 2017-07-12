LfmmDat.builder <- setRefClass("LfmmDat", contains = "Dat",
                               fields = c("X"),
                               methods = list(
                                 impute_lfmm = function(U, V, B) {
                                   impute_lfmm_cpp(.self$Y, .self$X, U, V, B, .self$missing.ind)
                                 },
                                 err2_lfmm = function(U, V, B) {
                                   err2_lfmm_cpp(.self$Y, .self$X, U, V, B)
                                 },
                                 sigma2_lm = function(X, B, nb.df) {
                                   if (is.matrix(.self$Y) && is.double(.self$Y)) {
                                     res <- sum2_lm_cpp(.self$Y, X, B) / nb.df
                                   } else {
                                     res <- foreach(j = 1:ncol(Y), .combine = 'c') %dopar%
                                       {
                                         aux <- Y[,j] - tcrossprod(X , B[j,,drop = FALSE])
                                         sum(aux * aux)
                                       }
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

