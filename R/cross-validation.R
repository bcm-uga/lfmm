#' return a list of train/test indices
#'
#' @export
left.out.kfold <- function(kfold, J) {
  if (kfold == 1) {
    cuts <- rep(factor("[1,m]"), J)
  } else {
    cuts <- cut(sample.int(J), breaks = kfold)
  }
  folds <- list()
  for (l in seq_along(levels(cuts))) {
    folds[[l]] <- which(cuts == levels(cuts)[l])
  }
  folds
}

CV <- function(m, dat, n.fold.row, n.fold.col, params) {

  n <- nrow(dat$Y)
  p <- ncol(dat$Y)

  param.names <- names(params)

  ## main loops
  res <- foreach(i = 1:nrow(params), .combine = 'rbind') %dopar%
    {
      errs <- tibble()

      param <- params[i, ]
      message("=== params")
      print.data.frame(param)

      m.train <- m

      ## copy dat object
      dat.train <- new(class(dat))
      dat.test <- new(class(dat))
      dat.predicted <- new(class(dat))

      ## row folds
      row.folds <- left.out.kfold(n.fold.row, n)
      for (row.fold in row.folds) {

        ## train/test
        dat.train$Y <- dat$Y[-row.fold,,drop = FALSE]
        dat.train$X <- dat$X[-row.fold,,drop = FALSE]
        dat.test$Y <- dat$Y[row.fold,,drop = FALSE]
        dat.test$X <- dat$X[row.fold,,drop = FALSE]

        ## method
        m.train[param.names] <- param

        ## fit method
        m.train <- MatrixFactorizationR_fit(m.train, dat.train)

        ## compute err
        col.folds <- left.out.kfold(n.fold.col, p)
        err <- data.frame()
        for (col.fold in col.folds) {
          m.train <- MatrixFactorizationR_fit_knowing_loadings(m = m.train,
                                                               dat = dat.test)
          dat.predicted$Y <- dat.test$Y
          dat.predicted$X <- dat.test$X
          dat.predicted$Y[,col.fold] <- NA
          dat.predicted <- MatrixFactorizationR_impute(m.train, dat.predicted)
          err <- rbind(err,
                       data.frame(err = mean((dat.predicted$Y[,col.fold] - dat.test$Y[,col.fold]) ^2),
                                  param))
        }
        errs <- rbind(errs, err)
      }
      errs
    }
  res
}
