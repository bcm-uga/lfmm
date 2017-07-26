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

CV <- function(m, dat, n.fold.row, n.fold.col, params, col.prop = 1.0, ...) {

  n <- nrow(dat$Y)
  p <- ncol(dat$Y)

  param.names <- names(params)

  m.train <- m
  ## main loops
  res <- foreach(i = 1:nrow(params), .combine = 'rbind') %dopar%
    {
      errs <- tibble()

            ## param
      param <- params[i, , drop = FALSE]
      message("=== params")
      print.data.frame(param)

      ## copy dat object
      dat.train <- new(class(dat))
      dat.test <- new(class(dat))

      ## row folds
      row.folds <- left.out.kfold(n.fold.row, n)
      for (row.fold in row.folds) {

        ## init U
        m.train$U <- NULL

        ## train/test
        dat.train$Y <- dat$Y[-row.fold,,drop = FALSE]
        dat.train$X <- dat$X[-row.fold,,drop = FALSE]
        dat.test$Y <- dat$Y[row.fold,,drop = FALSE]
        dat.test$X <- dat$X[row.fold,,drop = FALSE]

        ## method
        m.train[param.names] <- param

        ## fit method
        m.train <- MatrixFactorizationR_fit(m.train, dat.train, ...)

        ## col with less error
        lfmm.err2s <- dat.train$err2s_lfmm(m.train$U, m.train$V, m.train$B)
        kept.col.ind <- order(lfmm.err2s)[1:(round(col.prop * p))]

        ## compute err
        col.folds <- left.out.kfold(n.fold.col, length(kept.col.ind))
        err <- data.frame()
        for (col.fold in col.folds) {
          out.col.id <- kept.col.ind[col.fold]

          ## predict
          predicted.Y <- dat.test$predict_lfmm_knowing_loadings(V = m.train$V,
                                                                B = m.train$B,
                                                                unknown.j = out.col.id)
          ## compute error
          err <- rbind(err,
                       data.frame(err = mean((predicted.Y -
                                              dat.test$Y[,out.col.id]) ^2),
                                  param,
                                  nozero.prop = mean(m.train$B != 0.0)
                                  ))
        }
        errs <- rbind(errs, err)
      }
      errs
    }
  res
}
