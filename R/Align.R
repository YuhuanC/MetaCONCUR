
#' Shift Study-Specific Score Statistics to Aligned Score Statistic
#'
#' @description This function takes study-specific score statistics with length
#' equal to double the number of grids of a single study as input and returns
#' aligned score statistics with length equal to double the number of meta grids.
#'
#' @param SS A vector. The score statistic of a single study.
#' @param grid.study A vector. The grid boundaries of a single study.
#' @param grid.meta A vector. The grid boundaries of the combined studies.
#' @return A vector. The aligned score statistics.
#' @keywords internal

.dotMatch <- function(SS, grid.study, grid.meta) {

  stopifnot(
    "`SS` must be a numeric vector" =
      !missing(SS) && .isNumericVector(SS) && length(SS) > 1L,
    "`grid.study` must be a numeric vector" =
      !missing(grid.study) && .isNumericVector(grid.study) && length(grid.study) > 1L,
    "`grid.meta` must be a numeric vector" =
      !missing(grid.meta) && .isNumericVector(grid.meta) && length(grid.meta) > 1L,
    "`grid.meta` must contain `grid.study`" = all(grid.study %in% grid.meta))

  # n_bin_meta <- length(grid.meta) - 1L
  n_bin_study <- length(grid.study) - 1L

  index <- findInterval(grid.meta, grid.study)
  index <- index[-length(index)]

  # Taking care of zeros and out-of-boundary index
  store_zeros <- which(is.na(index) | {index == 0L} | {index > n_bin_study})
  index[store_zeros] <- NA

  SSm_dup <- SS[index]
  SSm_dup[is.na(SSm_dup)] <- 0.0

  SSm_del <- SS[index + n_bin_study]
  SSm_del[is.na(SSm_del)] <- 0.0

  c(SSm_dup, SSm_del)
}

# Switch for small matrix or big memory
.dotHH <- function(study.HH, ...) {
  UseMethod(".dotHH", study.HH)
}


#' Shift Study HH Matrix to Aligned HH Matrix
#'
#' @description This function takes study-specific summary statistic HxHx matrix
#' as input and returns an updated HxHx sparse matrix.
#'
#' @noRd
#' @param HH An HxHx Matrix for a single study.
#' @param grid.study A vector. The grid boundaries of a single study.
#' @param grid.meta A vector. The grid boundaries of the combined studies.
#'
#' @return A sparse Matrix. The updated HxHx sparse matrix.
#' @importFrom Matrix Matrix
#' @importFrom methods as
#'
#' @keywords internal
.dotHH.Matrix <- function(study.HH, grid.study, grid.meta) {

  stopifnot(
    "`grid.study` must be a numeric vector" =
      !missing(grid.study) && .isNumericVector(grid.study) && length(grid.study) > 1L,
    "`grid.meta` must be a numeric vector" =
      !missing(grid.meta) && .isNumericVector(grid.meta) && length(grid.meta) > 1L,
    "`grid.meta` must contain `grid.study`" = all(grid.study %in% grid.meta))

  # n_bin_meta <- length(grid.meta) - 1L
  n_bin_study <- length(grid.study) - 1L

  index <- findInterval(grid.meta, grid.study)
  index <- index[-length(index)]

  # Taking care of zeros and out-of-boundary index
  store_zeros <- which(is.na(index) | {index == 0L} | {index > n_bin_study})
  index[store_zeros] <- 1L

  # Update dup by dup
  HHm1 <- as(study.HH[index, index], "sparseMatrix")
  HHm1[store_zeros, ] <- Matrix(0, nrow = length(store_zeros), ncol = ncol(HHm1), sparse = TRUE)
  HHm1[, store_zeros] <- Matrix(0, nrow = nrow(HHm1), ncol = length(store_zeros), sparse = TRUE)

  # Update dup by del
  HHm2 <- as(study.HH[index, n_bin_study + index], "sparseMatrix")
  HHm2[store_zeros, ] <- Matrix(0, nrow = length(store_zeros), ncol = ncol(HHm2), sparse = TRUE)
  HHm2[, store_zeros] <- Matrix(0, nrow = nrow(HHm2), ncol = length(store_zeros), sparse = TRUE)

  # Update del by dup
  HHm3 <- as(study.HH[n_bin_study + index, index], "sparseMatrix")
  HHm3[store_zeros, ] <- Matrix(0, nrow = length(store_zeros), ncol = ncol(HHm3), sparse = TRUE)
  HHm3[, store_zeros] <- Matrix(0, nrow = nrow(HHm3), ncol = length(store_zeros), sparse = TRUE)

  # Update del by del
  HHm4 <- as(study.HH[index + n_bin_study, index + n_bin_study], "sparseMatrix")
  HHm4[store_zeros, ] <- Matrix(0, nrow = length(store_zeros), ncol = ncol(HHm4), sparse = TRUE)
  HHm4[, store_zeros] <- Matrix(0, nrow = nrow(HHm4), ncol = length(store_zeros), sparse = TRUE)

  as(rbind(cbind(HHm1, HHm2),
           cbind(HHm3, HHm4)), "sparseMatrix")

}

#' Shift Study HH Matrix to Meta Grid
#'
#' @noRd
#' @param HH An HxHx FBM object for a single study
#' @param grid.study grid of a single study
#' @param grid.meta grid of the combined studies
#'
#' @return An FBM object.
#' @import bigstatsr
#' @keywords internal
.dotHH.FBM <- function(study.HH, grid.study, grid.meta) {

  stopifnot(
    "`grid.study` must be a numeric vector" =
      !missing(grid.study) && .isNumericVector(grid.study) && length(grid.study) > 1L,
    "`grid.meta` must be a numeric vector" =
      !missing(grid.meta) && .isNumericVector(grid.meta) && length(grid.meta) > 1L,
    "`grid.meta` must contain `grid.study`" = all(grid.study %in% grid.meta))

  n_bin_meta <- length(grid.meta) - 1L
  n_bin_study <- length(grid.study) - 1L

  # vector of length n_bin_meta
  index <- findInterval(grid.meta, grid.study)
  index <- utils::head(index, -1L)

  # taking care of zeros and out-of-boundary index
  store_zeros <- which(is.na(index) | {index == 0L} | {index > n_bin_study})
  store_zeros_long <- c(store_zeros, store_zeros + n_bin_meta)
  index[store_zeros] <- 1L

  # print(paste0("Start creating matrix with dimension ", 2L*n_bin_meta, " by ", 2L*n_bin_meta))

  HHc <- bigstatsr::FBM(nrow = 2L * n_bin_meta,
                        ncol = 2L * n_bin_meta,
                        type = "double",
                        init = 0.0)

  #print("Matrix created")
  #gc()

  tmp <- bigstatsr::big_apply(
    X = HHc,
    a.FUN = function(X, ind, index, hh, hhc, store_zeros_long, n_bin_study, n_bin_meta) {
      # this fills the left side of the matrix
      if (any(ind <= n_bin_meta))
        hhc[, ind[ind <= n_bin_meta]] <- hh[c(index, index + n_bin_study), index[ind[ind <= n_bin_meta]]]

      # this fills the right side of the matrix
      if (any(ind > n_bin_meta))
        hhc[, ind[ind > n_bin_meta]] <- hh[c(index, index + n_bin_study),
                                           c(index, index + n_bin_study)[ind[ind > n_bin_meta]]]

      hhc[store_zeros_long, ind] <- 0.0
      if (any(ind %in% store_zeros_long))
        hhc[, ind[ind %in% store_zeros_long]] <- 0.0
      NULL
    },
    a.combine = 'c',
    block.size = 100,
    index = index,
    hh = study.HH,
    hhc = HHc,
    store_zeros_long = store_zeros_long,
    n_bin_study = n_bin_study,
    n_bin_meta = n_bin_meta
  )
  #gc()

  HHc

}

# Switch for small Matrix or big memory
.dotHX <- function(study.HX, ...) {
  UseMethod(".dotHX", study.HX)
}

#methods::setOldClass("Matrix")

#' Shift Study HX Matrix to Aligned HX Matrix
#'@description This function takes study-specific summary statistic HxX matrix
#' as input and returns an updated HxX sparse matrix.
#'
#' @noRd
#' @param HX An HxX Matrix for a single study
#' @param grid.study A vector. The grid boundaries of a single study.
#' @param grid.meta A vector. The grid boundaries of the combined studies.
#'
#' @return A Matrix.
#' @importFrom Matrix Matrix
#' @importFrom methods as
#' @keywords internal

.dotHX.Matrix <- function(study.HX, grid.study, grid.meta) {

  stopifnot(
    "`grid.study` must be a numeric vector" =
      !missing(grid.study) && .isNumericVector(grid.study) && length(grid.study) > 1L,
    "`grid.meta` must be a numeric vector" =
      !missing(grid.meta) && .isNumericVector(grid.meta) && length(grid.meta) > 1L,
    "`grid.meta` must contain `grid.study`" = all(grid.study %in% grid.meta))

  n_bin_meta <- length(grid.meta) - 1L
  n_bin_study <- length(grid.study) - 1L

  index <- findInterval(grid.meta, grid.study)
  index <- utils::head(index, -1L)

  # Taking care of zeros and out-of-boundary index
  store_zeros <- which(is.na(index) | {index == 0L} | {index > n_bin_study})
  index[store_zeros] <- 1L

  HXc <- as(study.HX[c(index, index + n_bin_study), ], "sparseMatrix")
  HXc[c(store_zeros, store_zeros + n_bin_meta), ] <-
    Matrix(0, nrow = length(c(store_zeros, store_zeros + n_bin_meta)),
           ncol = ncol(HXc), sparse = TRUE)
  HXc
}

#' Shift Study HX Matrix to Meta Grid
#' #'@description This function takes study-specific summary statistic HxX matrix
#' as input and returns an updated HxX sparse matrix.
#'
#' @noRd
#' @param HX An HxX FBM object for a single study
#' @param grid.study grid of a single study
#' @param grid.meta grid of the combined studies
#'
#' @return An FBM object.
#' @import bigstatsr
#' @keywords internal
.dotHX.FBM <- function(study.HX, grid.study, grid.meta) {

  stopifnot(
    "`grid.study` must be a numeric vector" =
      !missing(grid.study) && .isNumericVector(grid.study) && length(grid.study) > 1L,
    "`grid.meta` must be a numeric vector" =
      !missing(grid.meta) && .isNumericVector(grid.meta) && length(grid.meta) > 1L,
    "`grid.meta` must contain `grid.study`" = all(grid.study %in% grid.meta))

  n_bin_meta <- length(grid.meta) - 1L
  n_bin_study <- length(grid.study) - 1L

  index <- findInterval(grid.meta, grid.study)
  index <- utils::head(index, -1L)

  # Taking care of zeros and out-of-boundary index
  store_zeros <- which(is.na(index) | {index == 0L} | {index > n_bin_study})
  store_zeros_long <- c(store_zeros, store_zeros + n_bin_meta)
  index[store_zeros] <- 1L

  HXc <- bigstatsr::FBM(nrow = 2L*n_bin_meta, ncol = ncol(study.HX),
                        type = "double", init = 0.0)

  tmp <- bigstatsr::big_apply(
    X = study.HX,
    a.FUN = function(X, ind, index, store_zeros_long, HXc) {
      HXc[, ind] <- X[c(index, index + n_bin_study), ind]
      HXc[store_zeros_long, ind] <- 0.0
      NULL
    },
    a.combine = 'c',
    block.size = 100,
    index = index,
    store_zeros_long = store_zeros_long,
    HXc = HXc
  )

  HXc
}


#' Shift Study-Specific Summary Statistics to Aligned Summary Statistic
#'
#' @description This function takes a list of study-specific summary statistics,
#' including a study-specific grid vector, score statistics, covariance matrices,
#' as well as a meta grid vector as input and returns aligned summary statistics.
#'
#' @param study.data A list that contains `$grid` a vector of the grid boundaries;
#'  `$SSA` the score statistic for group A; `$SSB` the score statistic for group
#'  B; `$HAX` an 2Lxp matrix of H_A x; `$HBX` an 2Lxp matrix of H_B X;
#'  `$HAHA` an 2Lx2L matrix; `$HAHB` an 2Lx2L matrix; and `$HBHB` an 2Lx2L matrix.
#' @param meta.grid A vector. The grid boundaries of the combined studies.
#' @return A list. The aligned summary statistics.
#' @keywords internal
.alignStudyToMeta <- function(study.data, meta.grid) {

  stopifnot(
    "`study.data` must be a list of summary statistics in Worker-Output form" =
      !missing(study.data) && is.list(study.data) &&
      all(c("grid", "SSA", "SSB", "HAHA", "HAHB", "HBHB", "HAX", "HBX") %in% names(study.data)),
    "`meta.grid` must be a numeric vector" =
      !missing(meta.grid) && .isNumericVector(meta.grid) && length(meta.grid) > 1L,
    "`meta.grid` must contain `grid.study`" = all(study.data$grid %in% meta.grid))

  result <- list()

  study_grid <- study.data$grid

  result[["grid"]] <- study_grid
  # update score statistics
  result[["SSA"]] <-.dotMatch(study.data$SSA, study_grid, meta.grid)
  result[["SSB"]] <-.dotMatch(study.data$SSB, study_grid, meta.grid)

  # update variance-covariance matrix components
  # sparse matrices sum up to be numeric(0), can't get eigenvalue 0 x 0 matrix
  result[["HAHA"]] <-.dotHH(study.data$HAHA, study_grid, meta.grid)

  if (inherits (study.data$HAHA, "FBM")){
    #print(is(result[["HAHA"]]))
    #print(slotNames(result[["HAHA"]]))
    file.remove(study.data$HAHA$backingfile)
  }

  result[["HAHB"]] <-.dotHH(study.data$HAHB, study_grid, meta.grid)

  if (inherits (study.data$HAHB, "FBM")){

    file.remove(study.data$HAHB$backingfile)
  }

  result[["HBHB"]] <-.dotHH(study.data$HBHB, study_grid, meta.grid)

  if (inherits (study.data$HBHB, "FBM")){

    file.remove(study.data$HBHB$backingfile)
  }

  result[["HAX"]] <-.dotHX(study.data$HAX, study_grid, meta.grid)
  result[["HBX"]] <-.dotHX(study.data$HBX, study_grid, meta.grid)

  result
}

