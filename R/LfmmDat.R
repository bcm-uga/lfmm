LfmmDat.builder <- setRefClass("LfmmDat", contains = "Dat",
                               fields = c("X"))

#' Class which store data
#'
#'
#' @export
LfmmDat <- function(Y, X) {
  LfmmDat.builder(Y = read_input(Y),
                  X = read_input(X))
}

#' Class which store data
#'
#'
#' @export
SimulatedLfmmDat <- function(Y, X, outlier, U, V, B) {
  SimulatedLfmmDat.builder(Y = read_input(Y),
                           X = read_input(X),
                           outlier = read_input(outlier),
                           U = read_input(U),
                           B = read_input(B),
                           V = read_input(V))
}

