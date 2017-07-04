Dat.builder <- setRefClass("Dat", fields = c("Y", "meta", "missing.ind"),
                           methods = list(
                             getY = function() {
                               return(Y)
                             },
                             productY = function(x) {
                               Y %*% x
                             },
                             productYt = function(x) {
                               crossprod(Y, x)
                             },
                             svd = function(k, nu, nv, opts = list(tol = 10e-10)) {
                               Af <- function(x, args) {
                                 .self$productY(x)
                               }
                               Atransf <- function(x, args) {
                                 .self$productYt(x)
                               }
                               RSpectra::svds(A = Af,
                                              Atrans = Atransf,
                                              k = k,
                                              nu = nu, nv = nv,
                                              opts = opts,
                                              dim = c(nrow(.self$Y), ncol(.self$Y)))
                             },
                             svd_soft= function(gamma, k, opts = list(tol = 10e-10)) {
                               Af <- function(x, args) {
                                 .self$productY(x)
                               }
                               Atransf <- function(x, args) {
                                 .self$productYt(x)
                               }
                               svd.res <- RSpectra::svds(A = Af,
                                                         Atrans = Atransf,
                                                         k = k + 2, ## mouais...
                                                         nu = k + 2, nv = k + 2,
                                                         opts = opts,
                                                         dim = c(nrow(.self$Y), ncol(.self$Y)))
                               svd.res$d <- svd.res$d - gamma
                               svd.res$d <- svd.res$d[svd.res$d > 0.0]
                               K <- length(svd.res$d)
                               if (K > k) {
                                 warning("K is increasing, now K = ", K)
                               }
                               if (K > (k + 2) || K <= 0.0) {
                                 stop("K too big or too small, OMG, call 911 !!")
                               }
                               svd.res$u <-svd.res$u[,1:K] 
                               svd.res$v <- svd.res$v[,1:K]
                               svd.res
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
  dat
}

