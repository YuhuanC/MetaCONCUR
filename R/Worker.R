# Worker functions

#' Identify the records reflecting duplicates
#'
#' @noRd
#'
#' @param cnvData A data.frame object. Must contain at least column `TYPE`
#'   providing the copy number. All other columns are ignored.
#' @param cutoff A numeric object. Duplication is assumed if `TYPE is at or
#'   above the provided value.
#'
#' @return A logical vector indicating which records are duplications.
#'
#' @keywords internal
.idDuplicates <- function(cnvData, cutoff) {
  stopifnot(
    "`cnvData` must be a data.frame" =
      !missing(cnvData) && is.data.frame(cnvData),
    "`cnvData` must contain column `TYPE`" = "TYPE" %in% colnames(cnvData),
    "`cutoff` must be numeric" = !missing(cutoff) && is.numeric(cutoff))

  result <- cnvData$TYPE > {cutoff - 1e-8}
  if (all(result) || !any(result)) {
    warning(ifelse(all(result), "all ", "no "), "records code as duplication",
            call. = FALSE)
  }
  result
}

#' Identify the records reflecting deletions
#'
#' @noRd
#'
#' @param cnvData A data.frame object. Must contain at least column `TYPE`
#'   providing the copy number. All other columns are ignored.
#' @param cutoff A numeric object. Deletions is assumed if `TYPE is at or
#'   below the provided value.
#'
#' @return A logical vector indicating which records are deletions.
#'
#' @keywords internal
.idDeletions <- function(cnvData, cutoff) {
  stopifnot(
    "`cnvData` must be a data.frame"  =
      !missing(cnvData) && is.data.frame(cnvData),
    "`cnvData` must contain column `TYPE`" = "TYPE" %in% colnames(cnvData),
    "`cutoff` must be numeric" = !missing(cutoff) && is.numeric(cutoff))

  result <- cnvData$TYPE < {cutoff + 1e-8}
  if (all(result) || !any(result)) {
    warning(ifelse(all(result), "all ", "no "), "records code as deletion",
            call. = FALSE)
  }
  result
}


#' Identify the Coarsened Grid Boundaries
#'
#' Given start and end locations, identify the boundary values of the coarsened
#'   grid
#'
#' @noRd
#'
#' @param cnvData A data.frame object. Must contain column headers BP1 and BP2;
#'   all other columns are ignored.
#'
#' @return A vector of boundaries.
#'
#' @keywords internal
.defineGrid <- function(cnvData) {

  stopifnot(
    "`cnvData` must be a data.frame" = !missing(cnvData) && .isCNV(cnvData)
  )

  # grid is defined by the start/end points
  c(cnvData$BP1, cnvData$BP2) |> unique() |> sort()
}

#' Match the start and end locations to the coarsened grid units
#'
#' Match the start and end locations of each event to the units that define the
#'   coarsened grid.
#'
#' @noRd
#' @param cnvData A data.frame object. The CNV data. At a minimum this must
#'   contain column headers BP1 and BP2. All other columns are ignored.
#' @param grid A vector object. The boundaries that define the coarsened grid.
#'
#' @return A data.frame. The first column contains the ID of the event record.
#'   The second column (TYPE) contains the unaltered TYPE of the event record.
#'   The third column (lower) contains the index of the unit of the coarsened
#'   grid with lower boundary equal to the BP1 of the record. The fourth column
#'   (upper) contains the index of the unit of the coarsened grid with upper
#'   boundary equal to the BP2 of the record.
#'
#' @keywords internal
.matchGrid <- function(cnvData, grid) {
  stopifnot(
    "`cnvData` must be a data.frame" = !missing(cnvData) && .isCNV(cnvData),
    "`grid` must be a numeric vector" =
      !missing(grid) && .isNumericVector(grid) && length(grid) > 1L)

  data.frame("ID" = cnvData$ID,
             "TYPE" = cnvData$TYPE,
             "lower" = findInterval(cnvData$BP1, grid),
             "upper" = findInterval(cnvData$BP2, grid, left.open = TRUE))
}

#' Compute the Residual
#'
#' Calculates the residual of the provided linear model.
#'
#' @noRd
#' @param covData A data.frame object. The covariates.
#' @param phenoData A data.frame object. The phenotype data.
#' @param formula A formula object. The model. Must be of the form LHS ~ RHS.
#' @param phenoType A character object. Must be one of "bin" or "cont" indicating
#'   that that phenotype of interest is binary or continuous, respectively.
#'
#' @return A data.frame containing "ID", "Residual", "Fitted".
#'
#' @importFrom stats fitted.values glm lm residuals
#' @keywords internal
.getYstar <- function(covData, phenoData, formula, phenoType = c("bin", "cont")) {

  phenoType <- match.arg(phenoType)

  stopifnot(
    "`covData` must be a data.frame" =
      !missing(covData) && is.data.frame(covData),
    "`phenoData` must be a data.frame" =
      !missing(phenoData) && is.data.frame(phenoData),
    "`formula` must be a formula of the form LHS ~ RHS" =
      !missing(formula) && inherits(formula, "formula") && length(formula) == 3L,
    "`covData` and `phenoData` must contain the same IDs" =
      "ID" %in% colnames(covData) && "ID" %in% colnames(phenoData) &&
      all(covData$ID %in% phenoData$ID) && nrow(covData) == nrow(phenoData))

  id_order <- order(covData$ID)
  covData <- covData[id_order, ]
  phenoData <- phenoData[id_order, ]

  model <- tryCatch(switch(phenoType,
                           "cont" = stats::lm(formula, cbind(covData, phenoData)),
                           "bin" = stats::glm(formula, cbind(covData, phenoData),
                                              family = "binomial")),
                    error = function(e) {
                      stop("unable to fit phenotype\n\t", e$message,
                           call. = FALSE)
                    })

  res <- stats::residuals(model, type = "response")
  if (length(res) != nrow(covData)) {
    stop("a residual was not returned for all records", call. = FALSE)
  }

  data.frame("ID" = covData$ID,
             "Residuals" = res,
             "Fitted" = stats::fitted.values(model))
}

#' Compute Score Statistic on Coarsened Grid
#'
#' @noRd
#' @param Y A data.frame object. Contains the ID, Residuals, and Fitted values.
#' @param grid.ind A data.frame object. Contains the ID, TYPE, and the indices
#'   of the lower and upper grid units of each event record.
#' @param test.func A function. A function specifying the indicator function
#'   defining the H component under consideration. Note that this function
#'   takes as input the copy number vector, TYPE, and returns 0L/1L.
#' @param cutoff A numeric. The cutoff for the the indicator function specified
#'   by `test.func`.
#' @param px A numeric.The number of covariates plus one.
#' @param phenoType A character object. Must be one of "bin" or "cont" indicating
#'   that that phenotype of interest is binary or continuous, respectively.
#'
#' @return A vector. The score statistic on each grid unit. If a grid unit does
#'   note contain data, it will be zero.
#'
#' @keywords internal
.computeSS <- function(Y, grid.ind, test.func, cutoff, px,
                       phenoType = c("bin", "cont")) {

  phenoType <- match.arg(phenoType)

  stopifnot(
    "`Y` must be a data.frame with columns ID, Residuals, Fitted" = !missing(Y) &&
      is.data.frame(Y) && ncol(Y) == 3L &&
      all(c("ID", "Residuals", "Fitted") %in% colnames(Y)),
    "`grid.ind` must be a data.frame with columns lower, upper, TYPE, and ID" =
      !missing(grid.ind) && is.data.frame(grid.ind) &&
      all(c("lower", "upper", "ID", "TYPE") %in% colnames(grid.ind)),
    "`test.func` must be a function with inputs `x` and `cutoff`" = !missing(test.func) &&
      is.function(test.func) && length(names(formals(test.func))) == 2L,
    "`cutoff` must be a scalar" = !missing(cutoff) && .isNumericVector(cutoff, 1L),
    "IDs of `grid.ind` must be in IDs of `Y`" = all(grid.ind$ID %in% Y$ID))

  # calculate the indicator
  H <- do.call(test.func, list(grid.ind$TYPE, cutoff))

  ### Getting unbiased estimator of residual variance

  n <- nrow(Y)

  # normalization term
  phi <- switch(phenoType,
                "cont" = as.numeric(crossprod(Y$Residuals)/(n-px)) |> drop(),
                "bin" = 1.0)

  # shift Y to structure of the grid
  ids <- match(grid.ind$ID, Y$ID)
  Y_local <- Y$Residuals[ids]

  SS <- numeric(max(grid.ind$upper))
  for (i in seq_len(nrow(grid.ind))) {
    if (H[i] == 1L) {
      inds <- grid.ind$lower[i]:grid.ind$upper[i]
      SS[inds] <- SS[inds] +  Y_local[i] / phi
    }
  }
  SS
}

#' Calculate (X^T V^{-1} X)^{-1/2} Using Eigen Decomposition
#'
#' @noRd
#' @param X A data.frame object. The covariates. Must include column ID.
#' @param phi.hat A numeric vector. The variance.
#'
#' @return A matrix. (X^T V^{-1} X)^{-1/2}
#' @importFrom stats model.matrix update.formula
#' @keywords internal
.computeXTXinvsqrt <- function(X, phi.hat) {

  stopifnot(
    "`X` must be a matrix" = !missing(X) && is.matrix(X),
    "`phi.hat` must be a numeric vector" = !missing(phi.hat) &&
      .isNumericVector(phi.hat, nrow(X))
  )

  # calculate (X^T X)^{-1/2} using eigen decomposition
  xtvx_eigen <- tryCatch(eigen(crossprod(X, X / phi.hat), symmetric = TRUE),
                         error = function(e) {
                           stop("Unable to invert X^T V^{-1} X\n",
                                e$message, call. = FALSE)
                         })
  if (any(xtvx_eigen$values < 1e-16)) {

    d <- xtvx_eigen$values
    d <- {d + abs(d)} /2
    d2 <- 1.0 / sqrt(d)
    d2[abs(d) < 1e-16] <- 0.0

    xtx_inv_sqrt <- tryCatch(t(t(xtvx_eigen$vectors) * d2) %*%
                               solve(xtvx_eigen$vectors),
                             error = function(e) {
                               stop("Unable to calculate (X^T V^{-1} X)^{-1/2}\n",
                                    e$message, call. = FALSE)
                             })
  } else {
    # there should never be an error here, but JIC
    # p x p
    xtx_inv_sqrt <- tryCatch(t(t(xtvx_eigen$vectors) / sqrt(xtvx_eigen$values)) %*%
                               solve(xtvx_eigen$vectors),
                             error = function(e) {
                               stop("Unable to calculate (X^T V^{-1} X)^{-1/2}\n",
                                    e$message, call. = FALSE)
                             })
  }

  xtx_inv_sqrt
}

#' Compute Covariance Components on Coarsened Grid
#'
#' @noRd
#' @param y A data.frame object. Contains column ID, Fitted (the fitted values
#'   of the outcome regression), and Residual (the residuals of the outcome
#'   regression).
#' @param X A data.frame object. The covariates. Must include column ID.
#' @param xtx_inv_sqrt A matrix. (X^T V X)^{-1/2}
#' @param phi.hat A numeric vector. The variance.
#' @param grid.ind A data.frame object. Contains the ID, TYPE, and the indices
#'   of the lower and upper grid units of each event record.
#' @param HA A numeric vector. The indicator functions.
#' @param HB A numeric vector. The indicator functions.
#'
#' @return A list containing:
#'   \item{HAHA}{H_A^T V^{-1} H_A}
#'   \item{HAHB}{H_A^T V^{-1} H_B}
#'   \item{HBHB}{H_B^T V^{-1} H_B}
#'   \item{HAX}{H_A^T V^{-1} X (X^T V^{-1} X)^{-1/2}}
#'   \item{HBX}{H_B^T V^{-1} X (X^T V^{-1} X)^{-1/2}}
#'
#' @importFrom Matrix Matrix
#' @importFrom methods as
#' @importFrom stats model.matrix update.formula
#' @keywords internal
.computeCov.Matrix <- function(y, X, xtx_inv_sqrt, phi.hat, grid.ind, HA, HB) {

  stopifnot(
    "`y` must be a data.frame" = !missing(y) && is.data.frame(y) &&
      all(c("ID", "Fitted", "Residuals") %in% colnames(y)),
    "`X` must be a matrix" = !missing(X) && is.matrix(X),
    "`phi.hat` must be a numeric vector" = !missing(phi.hat) &&
      .isNumericVector(phi.hat, nrow(X)),
    "`grid.ind` must be a data.frame" = !missing(grid.ind) &&
      is.data.frame(grid.ind) &&
      all(c("ID", "TYPE", "lower", "upper") %in% colnames(grid.ind)),
    "`HA` must be an integer vector" = !missing(HA) &&
      .isIntegerVector(HA, 2L * nrow(grid.ind)),
    "`HB` must be a integer vector" = !missing(HB) &&
      .isIntegerVector(HB, 2L * nrow(grid.ind))
  )

  ids <- y$ID

  # Now need H^T V^{-1} X = 1/phi_hat H X and HT V^{-1} H

  # number of bins in coarse grid
  L <- max(grid.ind$upper)

  L2 <- 2L * L
  p <- ncol(X)

  # initialize covariance components matrices
  HAX <- Matrix::Matrix(0.0, nrow = L2, ncol = p, sparse = TRUE)
  HBX <- Matrix::Matrix(0.0, nrow = L2, ncol = p, sparse = TRUE)
  HAHA <- Matrix::Matrix(0.0, nrow = L2, ncol = L2, sparse = TRUE)
  HAHB <- Matrix::Matrix(0.0, nrow = L2, ncol = L2, sparse = TRUE)
  HBHB <- Matrix::Matrix(0.0, nrow = L2, ncol = L2, sparse = TRUE)

  # unique ID in the event data
  unique_ids <- unique(grid.ind$ID)
  grid.ind <- rbind(grid.ind,
                    data.frame("ID" = grid.ind$ID,
                               "TYPE" = grid.ind$TYPE,
                               "lower" = grid.ind$lower + L,
                               "upper" = grid.ind$upper + L))

  for (i in unique_ids) {

    # identify which rows of grid.ind pertain to the case under consideration
    # further limit to only those with H = 1
    cases_A <- grid.ind$ID == i & HA == 1L
    cases_B <- grid.ind$ID == i & HB == 1L

    if (any(cases_A)) {
      # get all grid units spanned by the events of this case
      idx <- sapply(which(cases_A),
                    function(j) {grid.ind$lower[j]:grid.ind$upper[j]}) |>
        unlist() |> c()

      use <- ids == i
      Xmat <- matrix(X[use, ]/ phi.hat[use], nrow = length(idx), ncol = p, byrow = TRUE)

      # update H^T H and H^T X
      HAHA[idx, idx] <- HAHA[idx, idx] + 1.0 / phi.hat[use]
      HAX[idx, ] <- HAX[idx, ] +  Xmat

      if (any(cases_B)) {
        idxB <- sapply(which(cases_B),
                       function(j) {grid.ind$lower[j]:grid.ind$upper[j]}) |>
          unlist() |> c()

        Xmat <- matrix(X[use, ]/ phi.hat[use], nrow = length(idxB), ncol = p, byrow = TRUE)

        HAHB[idx, idxB] <- HAHB[idx, idxB] + 1.0 / phi.hat[use]
        HBHB[idxB, idxB] <- HBHB[idxB, idxB] + 1.0 / phi.hat[use]
        HBX[idxB, ] <- HBX[idxB, ] + Xmat
      }

    } else if (any(cases_B)) {
      # get all grid units spanned by the events of this case
      idx <- sapply(which(cases_B),
                    function(j) {grid.ind$lower[j]:grid.ind$upper[j]}) |>
        unlist() |> c()
      use <- ids == i
      Xmat <- matrix(X[use, ]/ phi.hat[use], nrow = length(idx), ncol = p, byrow = TRUE)
      # update H^T H and H^T X
      HBHB[idx, idx] <- HBHB[idx, idx] + 1.0 / phi.hat[use]
      HBX[, idx] <- HBX[, idx] +  Xmat
    }
  }


  HAHA <- as(HAHA, "symmetricMatrix")
  HBHB <- as(HBHB, "symmetricMatrix")
  HAX2 <- HAX %*% xtx_inv_sqrt
  rm(HAX)
  HBX2 <- HBX %*% xtx_inv_sqrt
  rm(HBX)

  list("HAHA" = HAHA,
       "HAHB" = HAHB,
       "HBHB" = HBHB,
       "HAX" = HAX2,
       "HBX" = HBX2)
}

#' Compute Covariance Components on Coarsened Grid, Return FBM
#'
#' @noRd
#' @param y A data.frame object. Contains column ID, Fitted (the fitted values
#'   of the outcome regression), and Residual (the residuals of the outcome
#'   regression).
#' @param X A data.frame object. The covariates. Must include column ID.
#' @param xtx_inv_sqrt A matrix. (X^T V X)^{-1/2}
#' @param phi.hat A numeric vector. The variance.
#' @param grid.ind A data.frame object. Contains the ID, TYPE, and the indices
#'   of the lower and upper grid units of each event record.
#' @param HA A numeric vector. The indicator functions.
#' @param HB A numeric vector. The indicator functions.
#' @param save_dir A directory path where summary statistics will be saved.
#'
#' @return A list containing FBM objects:
#'   \item{HAHA}{H_A^T V^{-1} H_A}
#'   \item{HAHB}{H_A^T V^{-1} H_B}
#'   \item{HBHB}{H_B^T V^{-1} H_B}
#'   \item{HAX}{H_A^T V^{-1} X (X^T V^{-1} X)^{-1/2}}
#'   \item{HBX}{H_B^T V^{-1} X (X^T V^{-1} X)^{-1/2}}
#'
#' @import bigstatsr
#' @importFrom Matrix Matrix
#' @importFrom stats model.matrix update.formula
#' @keywords internal
.computeCov.FBM <- function(y, X, xtx_inv_sqrt, phi.hat, grid.ind, HA, HB, save_dir) {

  stopifnot(
    "`y` must be a data.frame" = !missing(y) && is.data.frame(y) &&
      all(c("ID", "Fitted", "Residuals") %in% colnames(y)),
    "`X` must be a matrix" = !missing(X) && is.matrix(X),
    "`phi.hat` must be a numeric vector" = !missing(phi.hat) &&
      .isNumericVector(phi.hat, nrow(X)),
    "`grid.ind` must be a data.frame" = !missing(grid.ind) &&
      is.data.frame(grid.ind) &&
      all(c("ID", "TYPE", "lower", "upper") %in% colnames(grid.ind)),
    "`HA` must be an integer vector" = !missing(HA) &&
      .isIntegerVector(HA, 2L * nrow(grid.ind)),
    "`HB` must be a integer vector" = !missing(HB) &&
      .isIntegerVector(HB, 2L * nrow(grid.ind))
  )

  ids <- y$ID

  # Now need H^T V^{-1} X = 1/phi_hat H X and HT V^{-1} H

  # number of bins in coarse grid
  L <- max(grid.ind$upper)

  L2 <- 2L * L
  p <- ncol(X)

  # initialize covariance components FBM
  HAX <- bigstatsr::FBM(nrow = L2, ncol = p,
                        type = "double", init = 0.0,
                        backingfile = file.path(save_dir, "HAX"))
  HBX <- bigstatsr::FBM(nrow = L2, ncol = p,
                        type = "double", init = 0.0,
                        backingfile = file.path(save_dir, "HBX"))
  HAHA <- bigstatsr::FBM(nrow = L2, ncol = L2,
                         type = "double", init = 0.0,
                         backingfile = file.path(save_dir, "HAHA"))
  HAHB <- bigstatsr::FBM(nrow = L2, ncol = L2,
                         type = "double", init = 0.0,
                         backingfile = file.path(save_dir, "HAHB"))
  HBHB <- bigstatsr::FBM(nrow = L2, ncol = L2,
                         type = "double", init = 0.0,
                         backingfile = file.path(save_dir, "HBHB"))

  # unique ID in the event data
  unique_ids <- unique(grid.ind$ID)
  grid.ind <- rbind(grid.ind,
                    data.frame("ID" = grid.ind$ID,
                               "TYPE" = grid.ind$TYPE,
                               "lower" = grid.ind$lower + L,
                               "upper" = grid.ind$upper + L))

  for (i in unique_ids) {

    # identify which rows of grid.ind pertain to the case under consideration
    # further limit to only those with H = 1
    cases_A <- grid.ind$ID == i & HA == 1L
    cases_B <- grid.ind$ID == i & HB == 1L

    if (any(cases_A)) {
      # get all grid units spanned by the events of this case
      idx <- sapply(which(cases_A),
                    function(j) {grid.ind$lower[j]:grid.ind$upper[j]}) |>
        unlist() |> c()

      use <- ids == i
      Xmat <- matrix(X[use, ]/ phi.hat[use], nrow = length(idx), ncol = p, byrow = TRUE)

      # update H^T H and H^T X
      HAHA[idx, idx] <- HAHA[idx, idx] + 1.0 / phi.hat[use]
      HAX[idx, ] <- HAX[idx, ] +  Xmat

      if (any(cases_B)) {
        idxB <- sapply(which(cases_B),
                       function(j) {grid.ind$lower[j]:grid.ind$upper[j]}) |>
          unlist() |> c()

        Xmat <- matrix(X[use, ]/ phi.hat[use], nrow = length(idxB), ncol = p, byrow = TRUE)

        HAHB[idx, idxB] <- HAHB[idx, idxB] + 1.0 / phi.hat[use]
        HBHB[idxB, idxB] <- HBHB[idxB, idxB] + 1.0 / phi.hat[use]
        HBX[idxB, ] <- HBX[idxB, ] + Xmat
      }

    } else if (any(cases_B)) {
      # get all grid units spanned by the events of this case
      idx <- sapply(which(cases_B),
                    function(j) {grid.ind$lower[j]:grid.ind$upper[j]}) |>
        unlist() |> c()
      use <- ids == i
      Xmat <- matrix(X[use, ]/ phi.hat[use], nrow = length(idx), ncol = p, byrow = TRUE)
      # update H^T H and H^T X
      HBHB[idx, idx] <- HBHB[idx, idx] + 1.0 / phi.hat[use]
      HBX[, idx] <- HBX[, idx] +  Xmat
    }
  }

  HAHA$save()
  HAHB$save()
  HBHB$save()

  # this is going to barf for large HAX/HBX
  HAX2 <- bigstatsr::FBM(nrow = nrow(HAX),
                         ncol = ncol(HAX),
                         type = "double", init = 0.0,
                         backingfile = file.path(save_dir, "HAX2"))

  tmp <- bigstatsr::big_apply(
    X = HAX,
    a.FUN = function(X, ind, xtx, hax) {
      hax[,] <- hax[,] + X[, ind] %*% xtx[ind, ]
      NULL
    },
    a.combine = 'c',
    block.size = 100,
    hax = HAX2,
    xtx = xtx_inv_sqrt
  )

  HAX2$save()

  HBX2 <- bigstatsr::FBM(nrow = nrow(HBX),
                         ncol = ncol(HBX),
                         type = "double", init = 0.0,
                         backingfile = file.path(save_dir, "HBX2"))
  tmp <- bigstatsr::big_apply(
    X = HBX,
    a.FUN = function(X, ind, xtx, hbx) {
      hbx[,] <- hbx[,] + X[, ind] %*% xtx[ind, ]
      NULL
    },
    a.combine = 'c',
    block.size = 100,
    hbx = HBX2,
    xtx = xtx_inv_sqrt
  )

  HBX2$save()

  list("HAHA" = HAHA,
       "HAHB" = HAHB,
       "HBHB" = HBHB,
       "HAX" = HAX2,
       "HBX" = HBX2)
}

#' Compute Covariance Components on Coarsened Grid
#'
#' @noRd
#' @param formula A formula object. The model. Must be of the form LHS ~ RHS.
#' @param X A data.frame object. The covariates. Must include column ID.
#' @param phi.hat A numeric vector. The variance.
#' @param grid.ind A data.frame object. Contains the ID, TYPE, and the indices
#'   of the lower and upper grid units of each event record.
#' @param HA A numeric vector. The indicator functions.
#' @param HB A numeric vector. The indicator functions.
#' @param use.big.memory A logical. If TRUE file-backed matrices are used to
#'   store returned matrices instead of Matrix objects
#' @param save_dir A directory path where summary statistics will be saved.
#'
#' @return A list containing:
#'   \item{HAHA}{H_A^T V^{-1} H_A}
#'   \item{HAHB}{H_A^T V^{-1} H_B}
#'   \item{HBHB}{H_B^T V^{-1} H_B}
#'   \item{HAX}{H_A^T V^{-1} X (X^T V^{-1} X)^{-1/2}}
#'   \item{HBX}{H_B^T V^{-1} X (X^T V^{-1} X)^{-1/2}}
#'
#' @import bigstatsr
#' @importFrom Matrix Matrix
#' @importFrom stats model.matrix update.formula
#' @keywords internal
.computeCov <- function(y, X, phi.hat, grid.ind, HA, HB, use.big.memory, save_dir) {

  stopifnot(
    "`y` must be a data.frame" = !missing(y) && is.data.frame(y) &&
      all(c("ID", "Fitted", "Residuals") %in% colnames(y)),
    "`X` must be a matrix" = !missing(X) && is.matrix(X),
    "`phi.hat` must be a numeric vector" = !missing(phi.hat) &&
      .isNumericVector(phi.hat, nrow(X)),
    "`grid.ind` must be a data.frame" = !missing(grid.ind) &&
      is.data.frame(grid.ind) &&
      all(c("ID", "TYPE", "lower", "upper") %in% colnames(grid.ind)),
    "`HA` must be an integer vector" = !missing(HA) &&
      .isIntegerVector(HA, 2L * nrow(grid.ind)),
    "`HB` must be an integer vector" = !missing(HB) &&
      .isIntegerVector(HB, 2L * nrow(grid.ind)),
    "`use.big.memory` must be a logical" = !missing(use.big.memory) &&
      .isLogicalVector(use.big.memory, 1L)
  )

  xtx_inv_sqrt <- .computeXTXinvsqrt(X, phi.hat)

  if (use.big.memory) {
    .computeCov.FBM(y, X, xtx_inv_sqrt, phi.hat, grid.ind, HA, HB, save_dir)
  } else{
    .computeCov.Matrix(y, X, xtx_inv_sqrt, phi.hat, grid.ind, HA, HB)
  }
}

#' Procedure for Estimating Score Statistic for Specific Event Type
#'
#' Given data, residuals, and the cutoff for the event type, estimates the
#'   score statistic.
#'
#' @noRd
#'
#' @param cnvData A data.frame object. The CNV event data.
#' @param cutoff.dup A numeric scalar. The cutoff at which TYPE values are deemed
#'   to be duplicates.
#' @param cutoff.del A numeric scalar. The cutoff at which TYPE values are deemed
#'   to be deletions
#' @param y A data.frame object. Contains column ID, Fitted (the fitted values
#'   of the outcome regression), and Residual (the residuals of the outcome
#'   regression).
#' @param X A matrix object. The design matrix of the outcome regression.
#' @param phenoType A character object. Must be one of "bin" or "cont" indicating
#'   that that phenotype of interest is binary or continuous, respectively.
#' @param use.big.memory A logical. If TRUE file-backed matrices are used to
#'   store H-based matrices instead of Matrix objects
#' @param save_dir A directory path where summary statistics will be saved.
#'
#' @return A list object containing `$grid` a vector of the grid boundaries;
#'  `$SSA` the score statistic for group A; `$SSB` the score statistic for group
#'  B; `$HAX` an Lxp matrix of H_A x; `$HBX` an Lxp matrix of H_B X;
#'  `$HAHA` an LxL matrix; `$HAHB` an LxL matrix; and `$HBHB` an LxL matrix.
#' @keywords internal
.scoreProcedure <- function(cnvData, cutoff.dup, cutoff.del,
                            y, X, phenoType = c("bin", "cont"),
                            use.big.memory, save_dir) {

  phenoType <- match.arg(phenoType)

  stopifnot(
    "`cnvData` must be a data.frame" = !missing(cnvData) && .isCNV(cnvData),
    "`cutoff.dup` must be a scalar" = !missing(cutoff.dup) && .isNumericVector(cutoff.dup, 1L),
    "`cutoff.del` must be a scalar" = !missing(cutoff.del) && .isNumericVector(cutoff.del, 1L),
    "`y` must be a data.frame" = !missing(y) && is.data.frame(y) &&
      all(c("ID", "Fitted", "Residuals") %in% colnames(y)),
    "`X` must be a matrix" = !missing(X) && is.matrix(X),
    "`use.big.memory` must be a logical" = !missing(use.big.memory) &&
      .isLogicalVector(use.big.memory, 1L)
  )

  # generate boundary points of coarsened grid
  L <- .defineGrid(cnvData)

  # match coarsened grid units to event records
  M <- .matchGrid(cnvData, L)

  # adding dup and del together
  func.dup <- function(x, cutoff.dup) { as.integer(x > {cutoff.dup - 1e-8}) }
  func.del <- function(x, cutoff.del) { as.integer(x < {cutoff.del + 1e-8}) }

  # calculate score
  px <- ncol(X)
  result <- list()
  result[["grid"]] <- L
  result[["SSA"]] <- c(.computeSS(y, M, func.dup, cutoff.dup, px, phenoType),
                       .computeSS(y, M, func.del, cutoff.del, px, phenoType))
  result[["SSB"]] <- c(.computeSS(y, M, func.dup, cutoff.dup + 1.0, px, phenoType),
                       .computeSS(y, M, func.del, cutoff.del - 1.0, px, phenoType))

  # Compute H outside the loop
  HA <- c(do.call(func.dup, list(M$TYPE, cutoff.dup)),
          do.call(func.del, list(M$TYPE, cutoff.del)))
  HB <- c(do.call(func.dup, list(M$TYPE, cutoff.dup + 1.0)),
          do.call(func.del, list(M$TYPE, cutoff.del - 1.0)))

  n <- nrow(y)

  phi.hat <- switch(phenoType,
                    "cont" = rep(as.numeric(crossprod(y$Residuals) / {n - px})|> drop(), nrow(y)),
                    "bin" = 1.0 / {y$Fitted * {1.0 - y$Fitted}})

  Convfun <- .computeCov(y, X, phi.hat, M, HA, HB, use.big.memory = use.big.memory,
                         save_dir)

  c(result, Convfun)

}

#' Compute Test statistics
#'
#' @param score Score statistic vector, SSA or SSB.
#' @param grid A vector of the grid boundaries.
#'
#' @return A scaler.
#' @keywords internal

.computeTS <- function(score, grid) {

  stopifnot(
    "`score` must be a numeric vector with double length of `grid`" =
      !missing(score) &&.isNumericVector(score, 2*(length(grid)-1)),
    "`grid` must be a numeric vector" =
      !missing(grid) && .isNumericVector(grid) && length(grid) > 1L)

  crossprod(score * diff(grid), score) |> drop()
}

#' Save summary statistics as portable FBM and RDS
#'
#' @param stats A list. Output from .scoreProcedure for specified chromosome.
#' @param chr_name A character. Name of specified chromosome, eg. chr1.
#' @param use.big.memory A logical. If TRUE, save "grid", "SSA", "SSB" as a rds
#'   file. If FALSE, specify save_dir to save stats as a rds file.
#' @param save_dir A directory path where summary statistics will be saved. If
#'   use.big.memory equal to TRUE, must provide a relative path, the folder will
#'   be created inside the current working directory (getwd()). Defaulted as NULL.
#'
#' @return Invisible NULL
#' @keywords internal
#'
.save_chr_stats <- function(stats, chr_name, use.big.memory = FALSE, save_dir = NULL) {

  fbm_names <- c("HAHA", "HAHB", "HBHB", "HAX", "HBX")
  rds_names <- c("grid", "SSA", "SSB")

  if (use.big.memory) {
    # Case 1: big.memory = TRUE, require save_dir
    stopifnot(!is.null(save_dir))
    chr_dir <- file.path(save_dir, chr_name)

    # FBMs are already persistent
    for (nm in fbm_names) {
      stopifnot(inherits(stats[[nm]], "FBM"))
    }
    # Save small objects together
    small_stats <- stats[rds_names]
    saveRDS(small_stats, file.path(chr_dir, "wo_FBM.rds"))

    message("Summary statistics (using big memory) saved to: ", normalizePath(chr_dir))

  } else {
    # Case 2: big.memory = FALSE
    if (!is.null(save_dir)) {
      chr_dir <- file.path(save_dir, chr_name)

      # Save the whole list in one .rds
      saveRDS(stats, file.path(chr_dir, "SS.rds"))

      message("Summary statistics saved to: ", normalizePath(chr_dir))
    } else {
      message("No save_dir provided, summary statistics not saved for ", chr_name)
    }
  }

  invisible(NULL)
}

#' Procedure for Estimating Score Statistic and Covariance Matrix Components
#'
#' Estimate the summary statistic and variance components for each chromosome.
#'
#' @param cnvData A data.frame object. The CNV event data.
#' @param phenoData A data.frame object. The phenotype data.
#' @param cutoff.dup A numeric scalar. The cutoff at which TYPE values are deemed
#'   to be duplicates.
#' @param cutoff.del A numeric scalar. The cutoff at which TYPE values are deemed
#'   to be deletions
#' @param formula A formula object. The model. Must be of the form LHS ~ RHS.
#' @param covData A data.frame object. The covariates. Must include column ID.
#' @param phenoType A character object. Must be one of "bin" or "cont" indicating
#'   that that phenotype of interest is binary or continuous, respectively.
#' @param chr A character of integer object. If character, must be "all" indicating
#'   all chromosomes present in `cnvData`. If integer, one or more chromosomes
#'   present in `cnvData`. The score statistic and variance components will be
#'   returned only for those chromosomes specified here.
#' @param use.big.memory A logical. If TRUE, large matrix objects are handled
#'   as file-backed matrices. Using the feature significantly impacts
#'   computation speed and should be used only when matrices are too large
#'   to store in RAM.
#' @param save_dir A directory path where summary statistics will be saved. If
#'   use.big.memory equal to TRUE, must provide a relative path, the folder will
#'   be created inside the current working directory (getwd()). For example,
#'   "study1" will create ~/study1 directory.
#'
#' @return A list object. Each element is itself a list corresponding to one or
#' more chromosomes specified in `chr` and contains:
#'   \describe{
#'     \item{grid}{A vector of the grid boundaries.}
#'     \item{SSA}{The score statistic for group A.}
#'     \item{SSB}{The score statistic for group B.}
#'     \item{HAHA}{An L x L matrix.}
#'     \item{HAHB}{An L x L matrix.}
#'     \item{HBHB}{An L x L matrix.}
#'     \item{HAX}{An L x p matrix of H_A x.}
#'     \item{HBX}{An L x p matrix of H_B x.}
#'   }
#' @export
#' @examples
#' # Load example data
#' data(cnvData)
#' data(phenoData)
#' data(covData)
#' # Run on example data, first 100 individuals
#' s1_cnv <- cnvData[as.numeric(sub("P", "", cnvData$ID)) >= 1 &
#'                    as.numeric(sub("P", "", cnvData$ID)) <= 100, ]
#' s1_pheno <- phenoData[as.numeric(sub("P", "", phenoData$ID)) >= 1 &
#'                        as.numeric(sub("P", "", phenoData$ID)) <= 100, ]
#' s1_cov <- covData[as.numeric(sub("P", "", covData$ID)) >= 1 &
#'                    as.numeric(sub("P", "", covData$ID)) <= 100, ]
#' cutoff.dup <- 2.3
#' cutoff.del <- 1.7
#' formula_bin <- PHEB ~ SEX
#' s1_SS <- CONCUR_worker(
#'   cnvData = s1_cnv,
#'   phenoData = s1_pheno,
#'   cutoff.dup = cutoff.dup,
#'   cutoff.del = cutoff.del,
#'   formula = formula_bin,
#'   covData = s1_cov,
#'   phenoType = "bin",
#'   chr = 1,
#'   use.big.memory = FALSE,
#'   save_dir = NULL
#' )
#' # Output is a list of summary statistics
#' names(s1_SS)
#' names(s1_SS$chr1)
CONCUR_worker <- function(cnvData, phenoData, cutoff.dup, cutoff.del,
                          formula, covData, phenoType = c("bin", "cont"),
                          chr = "all", use.big.memory = FALSE, save_dir = NULL) {
  phenoType <- match.arg(phenoType)

  stopifnot(
    "`cnvData` must be a data.frame with columns 'ID', 'CHR', 'BP1', 'BP2', and 'TYPE'" =
      !missing(cnvData) && .isCNV(cnvData),
    "`phenoData` must be a data.frame" = !missing(phenoData) &&
      is.data.frame(phenoData),
    "`cutoff.dup` must be a scalar" = !missing(cutoff.dup) && .isNumericVector(cutoff.dup, 1L),
    "`cutoff.del` must be a scalar" = !missing(cutoff.del) && .isNumericVector(cutoff.del, 1L),
    "`covData` must be a data.frame" =
      !missing(covData) && is.data.frame(covData),
    "`formula` must be a formula of the form LHS ~ RHS" =
      !missing(formula) && inherits(formula, "formula") && length(formula) == 3L,
    #"`covData` must contain the same IDs" = "ID" %in% colnames(covData),
    "`use.big.memory` must be a logical" = !missing(use.big.memory) &&
      .isLogicalVector(use.big.memory, 1L),
    "`save_dir` must be specified if using big memory" =
      !use.big.memory || {use.big.memory &&!is.null(save_dir)}
  )

  phenoData <- phenoData[order(phenoData$ID), ]
  covData <- covData[order(covData$ID), ]

  chr_in_data <- cnvData$CHR |> unique() |> sort()

  if (chr == "all") {
    chr <- chr_in_data
  } else {
    if (any(!(chr %in% chr_in_data))) {
      stop(chr[!chr %in% chr_in_data], " are not present in `cnvData`",
           call. = FALSE)
    }
    chr <- chr |> unique() |> sort()
  }

  # fit phenotype model
  y.star <- .getYstar(covData, phenoData, formula, phenoType)

  # removing response variable as model.matrix expects it and we don't need it
  formula <- stats::update.formula(formula, NULL ~ .)

  # keeping IDs to ensure that we can match across cov/cnv datasets
  ids <- covData[order(covData$ID),]$ID

  # covariates of model with rownames set as original ID values
  X <- stats::model.matrix(formula, covData[order(covData$ID),])

  # assign ID back to X matrix
  rownames(X) <- ids

  res <- list()

  for (ch in chr) {
    chr_name <- paste0("chr", ch)
    # generate a directory for each chromosome to store SS
    chr_dir <- if (!is.null(save_dir)) file.path(save_dir, chr_name) else NULL
    if (!is.null(chr_dir)){
      dir.create(chr_dir, recursive = TRUE, showWarnings = FALSE)
    }

    cnv_subset <- cnvData[cnvData$CHR == ch, ]
    cnv_subset <- cnv_subset[order(cnv_subset$ID, cnv_subset$BP1), ]

    chr_stats <- .scoreProcedure(cnv_subset,
                                 cutoff.dup,
                                 cutoff.del,
                                 y.star, X,
                                 phenoType,
                                 use.big.memory = use.big.memory,
                                 save_dir = chr_dir)
    # save SS for later share
    .save_chr_stats(chr_stats, chr_name, use.big.memory, save_dir)

    res[[chr_name]] <- chr_stats
  }

  res

}
