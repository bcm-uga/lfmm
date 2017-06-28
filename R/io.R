##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @param X input
##' @return output
##' @author cayek
##' @export
read_input <- function(X) {
  if (is.matrix(X)) {
    return(X)
  }

  if (is.data.frame(X)) {
    return(X)
  }

  if (is.numeric(X)) {
    return(X)
  }

  if (is.logical(X)) {
    return(X)
  }

  if (is.integer(X)) {
    return(X)
  }


  if (is.character(X)) {
    if (tools::file_ext(X) == "lfmm") {
      return(as.matrix(readr::read_delim(X, delim = " ",
                                         col_names = FALSE,
                                         col_types = readr::cols(.default = readr::col_integer()))
                       )
             )
    } else if (tools::file_ext(X) == "RData") {
      stop("TODO")
    } else if (tools::file_ext(X) == "rds") {
      return(readRDS(X))
    } else {
      stop("TODO")
    }
  }
  stop("X not handle")
}
