Dat.builder <- setRefClass("Dat", fields = c("Y", "meta"),
                           methods = list(
                             getY = function() {
                               return(Y)
                             },
                             productY = function(x) {
                               Y %*% x
                             },
                             productYt = function(x) {
                               crossprod(Y, x)
                             }
                           )
                           )

#' Class which store data
#'
#'
#' @export
Dat <- function(Y) {
  Dat.builder(Y = read_input(Y),
              meta = list())
}

