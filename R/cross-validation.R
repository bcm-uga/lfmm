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

CV <- function(m, dat, kfold.row, kfold.col, params) {

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

      ## row folds
      row.folds <- left.out.kfold(kfold.row, n)
      for (row.fold in row.folds) {

        ## train/test
        dat.train <- dat
        dat.train$Y <- dat.train$Y[-row.fold,,drop = FALSE]
        dat.train$X <- dat.train$X[-row.fold,,drop = FALSE]
        dat.test <- dat
        dat.test$Y <- dat.test$Y[row.fold,,drop = FALSE]
        dat.test$X <- dat.test$X[row.fold,,drop = FALSE]

        ## method
        m.train[param.names] <- param

        ## fit method
        m.train <- MatrixFactorizationR_fit(m.train, dat.train)

        ## compute err
        col.folds <- left.out.kfold(kfold.col, p)
        err <- data.frame()
        for (col.fold in col.folds) {
          m.train <- MatrixFactorizationR_fit_knowing_loadings(m = m.train,
                                                               dat = dat.test)
          dat.predicted <- dat.test
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
