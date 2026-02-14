# Sum aligned summary statistics across studies

#' Sum Aligned HH, HX Matrix across studies
#'
#' @description This function takes a list of aligned summary statistics as input,
#' sum all HxHx Matrices and all HxX Matrices separately for the same x (either A or B).
#'
#' @noRd
#' @param meta.data A list of aligned summary statistics from different studies.
#' @param x.name A character. Specify A or B.
#'
#' @return A list with a meta object named `meta` that stores summed HxHx, summed HxX.
#'
#' @import Matrix
#' @importFrom methods as
#'
#' @keywords internal
.sumHxHx.Matrix <- function(meta.data, x.name) {

  if (x.name == "A") {
    HH <- "HAHA"
    HX <- "HAX"
  } else {
    HH <- "HBHB"
    HX <- "HBX"
  }

  meta_obj <- which(names(meta.data) == "meta")
  n_studies <- length(meta.data)
  idx <- c(1L:n_studies)[-meta_obj]

  sum_HH <- meta.data[[idx[1L]]][[HH]] -
    as(Matrix::tcrossprod(meta.data[[idx[1L]]][[HX]]), "symmetricMatrix")
  meta.data[[idx[1L]]][[HH]] <- NULL

  for (i in idx[-1L]) {
    sum_HH <- sum_HH + meta.data[[i]][[HH]] -
      as(Matrix::tcrossprod(meta.data[[i]][[HX]]), "symmetricMatrix")
    meta.data[[i]][[HH]] <- NULL
  }

  meta.data$meta[[HH]] <- sum_HH
  meta.data
}

#' Sum Aligned HAHB Matrix across studies
#'
#' @description This function takes a list of aligned summary statistics as input,
#' sum all HAHB Matrices across studies.
#'
#' @noRd
#' @param meta.data A list of aligned summary statistics from different studies.
#'
#' @return A list with a meta object named `meta` that stores summed HAHB.
#'
#' @import Matrix
#'
#' @keywords internal
.sumHxHy.Matrix <- function(meta.data) {

  meta_obj <- which(names(meta.data) == "meta")
  n_studies <- length(meta.data)
  idx <- c(1L:n_studies)[-meta_obj]

  sum_HH <- meta.data[[idx[1L]]]$HAHB -
    Matrix::tcrossprod(meta.data[[idx[1L]]]$HAX, meta.data[[idx[1L]]]$HBX)
  meta.data[[idx[1L]]]$HAHB <- NULL

  for (i in idx[-1L]) {
    sum_HH <- sum_HH + meta.data[[i]]$HAHB -
      Matrix::tcrossprod(meta.data[[i]]$HAX, meta.data[[i]]$HBX)
    meta.data[[i]]$HAHB <- NULL
  }
  meta.data$meta$HAHB <- sum_HH
  meta.data
}

#' Sum Aligned HH, HX FBMs across studies
#'
#' @description This function takes a list of aligned summary statistics as input,
#' sum all HxHx FBMs and all HxX FBMs separately for the same x (either A or B).
#'
#' @noRd
#' @param meta.data A list of aligned summary statistics (FBMs) from different studies.
#' @param x.name A character. Specify A or B.
#'
#' @return A list with a meta object named `meta` that stores summed HxHx, summed HxX as FBMs.
#'
#' @import bigstatsr
#'
#' @keywords internal
.sumHxHx.FBM <- function(meta.data, x.name) {

  if (x.name == "A") {
    HH <- "HAHA"
    HX <- "HAX"
  } else {
    HH <- "HBHB"
    HX <- "HBX"
  }

  meta_obj <- which(names(meta.data) == "meta")
  n_studies <- length(meta.data)
  idx <- c(1L:n_studies)[-meta_obj]

  xx <- bigstatsr::big_tcrossprodSelf(meta.data[[idx[1L]]][[HX]])

  sum_HH <- bigstatsr::FBM(nrow = nrow(xx), ncol = ncol(xx),
                           type = "double", init = 0.0)

  tmp <- bigstatsr::big_apply(
    X = meta.data[[idx[1L]]][[HH]],
    a.FUN = function(X, ind, xx, sum_HH) {
      sum_HH[, ind] <- X[, ind] - xx[, ind]
      NULL
    },
    a.combine = 'c',
    block.size = 100,
    xx = xx,
    sum_HH = sum_HH
  )

  meta.data[[idx[1L]]][[HH]] <- NULL

  for (i in idx[-1L]) {
    xx <- bigstatsr::big_tcrossprodSelf(meta.data[[i]][[HX]])

    tmp <- bigstatsr::big_apply(
      X = meta.data[[i]][[HH]],
      a.FUN = function(X, ind, sum_hh, xx) {
        sum_hh[, ind] <- sum_hh[, ind] + X[, ind] - xx[, ind]
        NULL
      },
      a.combine = 'c',
      block.size = 100,
      sum_hh = sum_HH,
      xx = xx
    )
    meta.data[[i]][[HH]] <- NULL
  }

  meta.data$meta[[HH]] <- sum_HH
  meta.data
}

#' Sum Aligned HAHB FBMs across studies
#'
#' @description This function takes a list of aligned summary statistics as input,
#' sum all HAHB FBMs.
#'
#' @noRd
#' @param meta.data A list of aligned summary statistics (FBMs) from different studies.
#'
#' @return A list with a meta object named `meta` that stores summed HAHB as FBMs.
#'
#' @import bigstatsr
#'
#' @keywords internal
.sumHxHy.FBM <- function(meta.data) {

  meta_obj <- which(names(meta.data) == "meta")
  n_studies <- length(meta.data)
  idx <- c(1L:n_studies)[-meta_obj]

  hbt <- bigstatsr::big_transpose(meta.data[[idx[1L]]]$HBX)

  sum_HH <- bigstatsr::FBM(nrow = nrow(meta.data[[idx[1L]]]$HAHB),
                           ncol = ncol(meta.data[[idx[1L]]]$HAHB),
                           type = "double", init = 0.0)

  tmp <- bigstatsr::big_apply(
    X = meta.data[[idx[1L]]]$HAHB,
    a.FUN = function(X, ind, hax, hbt, sum_HH) {
      sum_HH[, ind] <- X[, ind] - hax %*% hbt[, ind]
      NULL
    },
    a.combine = 'c',
    block.size = 100,
    hax = meta.data[[idx[1L]]]$HAX,
    hbt = hbt,
    sum_HH = sum_HH
  )
  meta.data[[idx[1L]]]$HAHB <- NULL

  for (i in idx[-1L]) {
    hbt <- bigstatsr::big_transpose(meta.data[[i]]$HBX)
    tmp <- bigstatsr::big_apply(
      X = meta.data[[i]]$HAHB,
      a.FUN = function(X, ind, sum_hh, hax, hbt) {
        sum_hh[, ind] <- sum_hh[, ind] + X[, ind] - hax %*% hbt[, ind]
        NULL
      },
      a.combine = 'c',
      block.size = 100,
      sum_hh = sum_HH,
      hax = meta.data[[i]]$HAX,
      hbt = hbt
    )
    meta.data[[i]]$HAHB <- NULL
  }

  meta.data$meta$HAHB <- sum_HH
  meta.data
}

# Switch for small HAHA/HAX Matrix or big memory HAHA/HAX FBM
.sumHAHA <- function(meta.data) {

  meta_obj <- which(names(meta.data) == "meta")
  n_studies <- length(meta.data)
  idx <- 1L:n_studies
  if (length(meta_obj) != 0) idx <- idx[-meta_obj]
  if (inherits(meta.data[[idx[1L]]]$HAHA, "FBM")) {
    .sumHxHx.FBM(meta.data, "A")
  } else {
    .sumHxHx.Matrix(meta.data, "A")
  }
}

# Switch for small HAHB Matrix or big memory HAHB FBM
.sumHAHB <- function(meta.data) {

  meta_obj <- which(names(meta.data) == "meta")
  n_studies <- length(meta.data)
  idx <- c(1L:n_studies)[-meta_obj]
  if (inherits(meta.data[[idx[1L]]]$HAHB, "FBM")) {
    .sumHxHy.FBM(meta.data)
  } else {
    .sumHxHy.Matrix(meta.data)
  }
}

# Switch for small HBHB/HBX Matrix or big memory HBHB/HBX FBM
.sumHBHB <- function(meta.data) {

  meta_obj <- which(names(meta.data) == "meta")
  n_studies <- length(meta.data)
  idx <- c(1L:n_studies)[-meta_obj]
  if (inherits(meta.data[[idx[1L]]]$HBHB, "FBM")) {
    .sumHxHx.FBM(meta.data, "B")
  } else {
    .sumHxHx.Matrix(meta.data, "B")
  }
}

# remove input study-specific SS
.rmStudies <- function(meta.data) {

  meta_obj <- which(names(meta.data) == "meta")
  n_studies <- length(meta.data)
  idx <- c(1L:n_studies)[-meta_obj]
  meta.data[idx] <- NULL
  meta.data
}

