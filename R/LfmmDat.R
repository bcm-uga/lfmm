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

