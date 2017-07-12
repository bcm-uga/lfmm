Dat.builder <- setRefClass("Dat", fields = c("Y", "meta", "missing.ind"),
                           methods = list(
                             getY = function() {
                               return(.self$Y)
                             },
                             productY = function(x) {
                               .self$Y %*% x
                             },
                             productYt = function(x) {
                               crossprod(.self$Y, x)
                             }
                           )
                           )

#' Class which store data
#'
#'
#' @export
Dat <- function(Y) {
  dat <- Dat.builder(Y = read_input(Y),
                     meta = list(),
                     missing.ind = NULL)
  dat$missing.ind <- which(is.na(dat$Y))
  dat
}

