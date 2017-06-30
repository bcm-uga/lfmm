LfmmDat.builder <- setRefClass("LfmmDat", contains = "Dat",
                               fields = c("X"),
                               methods = list(
                                 impute_lfmm = function(U, V, B) {
                                   impute_lfmm_cpp(.self$Y, X, U, V, B, .self$missing.ind)
                                 },
                                 err2_lfmm = function(U, V, B) {
                                   err2_lfmm_cpp(Y, X, U, V, B)
                                 }
                               )
                               )

#' Class which store data
#'
#'
#' @export
LfmmDat <- function(Y, X) {
  dat <- LfmmDat.builder(Y = read_input(Y),
                         X = read_input(X),
                         missing.ind = NULL,
                         meta = list())
  dat$missing.ind <- which(is.na(dat$Y))
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

