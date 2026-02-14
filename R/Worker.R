# Worker functions

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
    "`cnvData` is not properly specified" = !missing(cnvData) && 
      inherits(cnvData, "cnvData")
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
#' @include helpful_tests.R
#' @keywords internal
.matchGrid <- function(cnvData, grid) {
  stopifnot(
    "`cnvData` is not properly specified" = !missing(cnvData) && 
      inherits(cnvData, "cnvData"),
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

  model <- tryCatch(
    switch(phenoType,
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
#' @include helpful_tests.R
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

  # calculate the indicator - note that this is now a logical
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

  SS <- numeric(max(grid.ind$upper) + 1L)
  for (i in seq_len(nrow(grid.ind))) {
    if (H[i]) {
      val <- Y_local[i] / phi
      SS[grid.ind$lower[i]] <- SS[grid.ind$lower[i]] + val
      SS[grid.ind$upper[i] + 1L] <- SS[grid.ind$upper[i] + 1L] - val
    }
  }
  SS <- cumsum(SS)
  head(SS, -1L)
}

#' Calculate (X^T V^{-1} X)^{-1/2} Using Eigen Decomposition
#'
#' @noRd
#' @param X A data.frame object. The covariates. Must include column ID.
#' @param phi.hat A numeric vector. The variance.
#'
#' @return A matrix. (X^T V^{-1} X)^{-1/2}
#' 
#' @importFrom stats model.matrix update.formula
#' @include helpful_tests.R
#' @keywords internal
.computeXTXinvsqrt <- function(X, phi.hat) {

  stopifnot(
    "`X` must be a matrix" = !missing(X) && is.matrix(X),
    "`phi.hat` must be a numeric vector" = !missing(phi.hat) &&
      .isNumericVector(phi.hat, nrow(X))
  )

  # calculate (X^T X)^{-1/2} using eigen decomposition
  Xw <- X / sqrt(phi.hat)
  A  <- crossprod(Xw)
  xtvx_eigen <- eigen(A, symmetric = TRUE)
  
  Q <- xtvx_eigen$vectors
  d <- xtvx_eigen$values
  
  if (any(xtvx_eigen$values < 1e-16)) {
    d <- {d + abs(d)} / 2.0
    d2 <- 1.0 / sqrt(d)
    d2[abs(d) < 1e-16] <- 0.0
  } else {
    d2 <- 1.0 / sqrt(d)
  }

  Qscaled <- sweep(Q, 2L, d2, `*`)
  tcrossprod(Qscaled, Q)
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
    "`HA` must be a logical vector" = !missing(HA) &&
      .isLogicalVector(HA, nrow(grid.ind)),
    "`HB` must be a logical vector" = !missing(HB) &&
      .isLogicalVector(HB, nrow(grid.ind))
  )

  ids <- y$ID

  # Now need H^T V^{-1} X = 1/phi_hat H X and HT V^{-1} H

  # number of bins in coarse grid
  L <- max(grid.ind$upper)

  p <- ncol(X)

  # initialize covariance components matrices
  HAX <- Matrix::Matrix(0.0, nrow = L, ncol = p, sparse = TRUE)
  HBX <- Matrix::Matrix(0.0, nrow = L, ncol = p, sparse = TRUE)
  HAHA <- Matrix::Matrix(0.0, nrow = L, ncol = L, sparse = TRUE)
  HAHB <- Matrix::Matrix(0.0, nrow = L, ncol = L, sparse = TRUE)
  HBHB <- Matrix::Matrix(0.0, nrow = L, ncol = L, sparse = TRUE)

  # unique ID in the event data
  unique_ids <- unique(grid.ind$ID)

  # group the row indices of grid.ind that pertain to each ID
  rows_by_id <- split(seq_len(nrow(grid.ind)), grid.ind$ID)

  for (i in unique_ids) {

    # identify which rows of grid.ind pertain to the case under consideration
    # further limit to only those with H = 1
    rows <- rows_by_id[[i]]
    cases_A <- rows[HA[rows]]
    cases_B <- rows[HB[rows]]
    
    if (length(cases_A) + length(cases_B) > 0L) {
      use <- match(i, ids)
      if (is.na(use)) next
      v <- 1.0 / phi.hat[use]
      Xvec <- X[use, ]
    } else {
      next
    }
    
    if (length(cases_A) > 0L) {
      # get all grid units spanned by the events of this case
      idx <- unlist(Map(seq, grid.ind$lower[cases_A], grid.ind$upper[cases_A]), 
                    use.names = FALSE)
      uA <- Matrix::sparseVector(x = 1.0, i = idx, length = L)
      
      HAHA <- HAHA + v * tcrossprod(uA)

      HAX <- HAX + tcrossprod(uA, Xvec * v)

      if (length(cases_B) > 0L) {
        idxB <- unlist(Map(seq, grid.ind$lower[cases_B], grid.ind$upper[cases_B]), 
                       use.names = FALSE)
        uB <- Matrix::sparseVector(x = 1.0, i = as.numeric(idxB), length = L)

        HAHB <- HAHB + v * tcrossprod(uA, uB)
        HBHB <- HBHB + v * tcrossprod(uB)
        HBX <- HBX + tcrossprod(uB, Xvec * v)
      }

    } else if (length(cases_B) > 0L) {
      # get all grid units spanned by the events of this case
      idx <- unlist(Map(seq, grid.ind$lower[cases_B], grid.ind$upper[cases_B]), 
                    use.names = FALSE)
      uB <- Matrix::sparseVector(x = 1.0, i = idx, length = L)

      HBHB <- HBHB + v * tcrossprod(uB)
      HBX <- HBX + tcrossprod(uB, Xvec)
    }
  }


  HAHA <- Matrix::forceSymmetric(HAHA, uplo = "U")
  HBHB <- Matrix::forceSymmetric(HBHB, uplo = "U")
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
#' @param HA A logical vector. The indicator functions.
#' @param HB A logical vector. The indicator functions.
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
      .isLogicalVector(HA, nrow(grid.ind)),
    "`HB` must be a integer vector" = !missing(HB) &&
      .isLogicalVector(HB, nrow(grid.ind))
  )

  ids <- y$ID

  # Now need H^T V^{-1} X = 1/phi_hat H X and HT V^{-1} H

  # number of bins in coarse grid
  L <- max(grid.ind$upper)

  p <- ncol(X)

  # initialize covariance components FBM
  HAX <- bigstatsr::FBM(nrow = L, ncol = p,
                        type = "double", init = 0.0,
                        backingfile = file.path(save_dir, "HAX"))
  HBX <- bigstatsr::FBM(nrow = L, ncol = p,
                        type = "double", init = 0.0,
                        backingfile = file.path(save_dir, "HBX"))
  HAHA <- bigstatsr::FBM(nrow = L, ncol = L,
                         type = "double", init = 0.0,
                         backingfile = file.path(save_dir, "HAHA"))
  HAHB <- bigstatsr::FBM(nrow = L, ncol = L,
                         type = "double", init = 0.0,
                         backingfile = file.path(save_dir, "HAHB"))
  HBHB <- bigstatsr::FBM(nrow = L, ncol = L,
                         type = "double", init = 0.0,
                         backingfile = file.path(save_dir, "HBHB"))

  # unique ID in the event data
  unique_ids <- unique(grid.ind$ID)

  # group the row indices of grid.ind that pertain to each ID
  rows_by_id <- split(seq_len(nrow(grid.ind)), grid.ind$ID)
  
  for (i in unique_ids) {
    
    # identify which rows of grid.ind pertain to the case under consideration
    # further limit to only those with H = 1
    rows <- rows_by_id[[i]]
    cases_A <- rows[HA[rows]]
    cases_B <- rows[HB[rows]]
    
    if (length(cases_A) + length(cases_B) > 0L) {
      use <- match(i, ids)
      if (is.na(use)) next
      v <- 1.0 / phi.hat[use]
      Xvec <- X[use, , drop = FALSE]
    } else {
      next
    }
    
    if (length(cases_A) > 0L) {
      # get all grid units spanned by the events of this case
      idx <- unlist(Map(seq, grid.ind$lower[cases_A], grid.ind$upper[cases_A]), 
                    use.names = FALSE)

      # update H^T H and H^T X
      HAHA[idx, idx] <- HAHA[idx, idx] + v
      HAX[idx, ] <- HAX[idx, ] +  Xvec[rep(1L, length(idx)), ] * v

      if (length(cases_B) > 0L) {
        idxB <- unlist(Map(seq, grid.ind$lower[cases_B], grid.ind$upper[cases_B]), 
                       use.names = FALSE)

        HAHB[idx, idxB] <- HAHB[idx, idxB] + v
        HBHB[idxB, idxB] <- HBHB[idxB, idxB] + v
        HBX[idxB, ] <- HBX[idxB, ] +  Xvec[rep(1L, length(idxB)), ] * v
      }

    } else if (length(cases_B) > 0L) {
      # get all grid units spanned by the events of this case
      idx <- unlist(Map(seq, grid.ind$lower[cases_B], grid.ind$upper[cases_B]), 
                    use.names = FALSE)
      
      # update H^T H and H^T X
      HBHB[idx, idx] <- HBHB[idx, idx] + v
      HBX[idx,] <- HBX[idx, ] +   Xvec[rep(1L, length(idx)), ]
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
  
  for (i in seq(1L, L, by = 100)) {
    rows <- i:min(i + 100 - 1L, L)
    HAX2[rows, ] <- HAX[rows, ] %*% xtx_inv_sqrt
  }
  
  HAX2$save()

  HBX2 <- bigstatsr::FBM(nrow = nrow(HBX),
                         ncol = ncol(HBX),
                         type = "double", init = 0.0,
                         backingfile = file.path(save_dir, "HBX2"))
  
  for (i in seq(1L, L, by = 100)) {
    rows <- i:min(i + 100 - 1L, L)
    HBX2[rows, ] <- HBX[rows, ] %*% xtx_inv_sqrt
  }
  
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
    "`HA` must be an logical vector" = !missing(HA) &&
      .isLogicalVector(HA, nrow(grid.ind)),
    "`HB` must be an logical vector" = !missing(HB) &&
      .isLogicalVector(HB, nrow(grid.ind)),
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

#' Steps for obtaining smaller dimension components for DUP or DEL
#'
#' @noRd
#'
#' @param L A vector. Grid boundaries.
#' @param M A data frame. Contains ID, type and index for lower and upper grid
#'   boundaries.
#' @param cutoff A numeric scalar. The cutoff at which TYPE values are deemed
#'   to be duplicates or deletions
#' @param y A data.frame object. Contains column ID, Fitted (the fitted values
#'   of the outcome regression), and Residual (the residuals of the outcome
#'   regression).
#' @param X A matrix object. The design matrix of the outcome regression.
#' @param phenoType A character object. Must be one of "bin" or "cont" indicating
#'   that that phenotype of interest is binary or continuous, respectively.
#' @param event A character object. Must be one of "dup" or "del" indicating that
#'   the CNV event of interest is duplication or deletion, respectively.
#' @param use.big.memory A logical. If TRUE file-backed matrices are used to
#'   store H-based matrices instead of Matrix objects.
#' @param save_dir A directory path where summary statistics will be saved.
#'
#' @return A list object containing `$grid` a vector of the grid boundaries;
#'  `$SSA` the score statistic for group A; `$SSB` the score statistic for group
#'  B; `$HAX` an Lxp matrix of H_A x; `$HBX` an Lxp matrix of H_B X;
#'  `$HAHA` an LxL matrix; `$HAHB` an LxL matrix; and `$HBHB` an LxL matrix.
#' @keywords internal
.subscoreProcedure <- function(L, M, cutoff,
                               y, X, phenoType = c("bin", "cont"),
                               event = c("dup","del"),
                               use.big.memory, save_dir) {
  
  phenoType <- match.arg(phenoType)
  event <- match.arg(event)
  
  stopifnot(
    #"`cnvData` must be a data.frame" = !missing(cnvData) && .isCNV(cnvData),
    "`cutoff` must be a scalar" = !missing(cutoff) && .isNumericVector(cutoff, 1L),
    "`y` must be a data.frame" = !missing(y) && is.data.frame(y) &&
      all(c("ID", "Fitted", "Residuals") %in% colnames(y)),
    "`X` must be a matrix" = !missing(X) && is.matrix(X),
    "`use.big.memory` must be a logical" = !missing(use.big.memory) &&
      .isLogicalVector(use.big.memory, 1L)
  )
  
  # switch function and cutoff for SSB based on dup or del
  func <- switch(event,
                 "dup" = function(x, cutoff) { x > {cutoff - 1e-8} },
                 "del" = function(x, cutoff) { x < {cutoff + 1e-8} })
  
  cutoff2 <- switch(event,
                    "dup" = cutoff + 1.0,
                    "del" = cutoff - 1.0)
  
  # calculate score
  px <- ncol(X)
  result <- list()
  result[["grid"]] <- L
  result[["SSA"]] <- .computeSS(y, M, func, cutoff, px, phenoType)
  result[["SSB"]] <- .computeSS(y, M, func, cutoff2, px, phenoType)
  
  # Compute H outside the loop -- note that HA and HB are now Logical
  HA <- do.call(func, list(M$TYPE, cutoff))
  HB <- do.call(func, list(M$TYPE, cutoff2))

  n <- nrow(y)
  
  phi.hat <- switch(phenoType,
                    "cont" = rep(as.numeric(crossprod(y$Residuals) / {n - px})|> drop(), nrow(y)),
                    "bin" = 1.0 / {y$Fitted * {1.0 - y$Fitted}})
  
  Convfun <- .computeCov(y, X, phi.hat, M, HA, HB,
                         use.big.memory = use.big.memory, save_dir)
  
  c(result, Convfun)
  
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
#' @param save_dir.dup A directory path where summary statistics of "DUP" will be saved.
#' @param save_dir.del A directory path where summary statistics of "DEL" will be saved.
#'
#' @return A list object containing two lists: "DUP" and "DEL", each list contains
#' the output of .subscoreProcedure.
#' @keywords internal
.scoreProcedure <- function(cnvData, cutoff.dup, cutoff.del,
                            y, X, phenoType = c("bin", "cont"),
                            use.big.memory, save_dir.dup, save_dir.del) {
  
  # generate boundary points of coarsened grid
  L <- .defineGrid(cnvData)
  
  # match coarsened grid units to event records
  M <- .matchGrid(cnvData, L)
  
  res_DUP <- .subscoreProcedure(L, M, cutoff.dup,
                                y, X, phenoType, event = "dup",
                                use.big.memory, save_dir.dup)
  
  res_DEL <- .subscoreProcedure(L, M, cutoff.del,
                                y, X, phenoType, event = "del",
                                use.big.memory, save_dir.del)
  
  list("DUP" = res_DUP,
       "DEL" = res_DEL)
  
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
    "`score` must be a numeric vector with length of `grid`" =
      !missing(score) &&.isNumericVector(score, (length(grid)-1)),
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
#'     \item{HAHA}{An 2L x 2L matrix.}
#'     \item{HAHB}{An 2L x 2L matrix.}
#'     \item{HBHB}{An 2L x 2L matrix.}
#'     \item{HAX}{An 2L x p matrix of H_A x.}
#'     \item{HBX}{An 2L x p matrix of H_B x.}
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
  
  class(cnvData) <- c(class(cnvData), "cnvData")

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
    
    chr_dir.dup <- if (!is.null(chr_dir)) file.path(chr_dir, "DUP") else NULL
    
    chr_dir.del <- if (!is.null(chr_dir)) file.path(chr_dir, "DEL") else NULL
    
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
                                 save_dir.dup = chr_dir.dup,
                                 save_dir.del = chr_dir.del)
    # save SS for later share
    .save_chr_stats(chr_stats$DUP, chr_name, use.big.memory, chr_dir.dup)
    .save_chr_stats(chr_stats$DEL, chr_name, use.big.memory, chr_dir.del)
    
    res[[chr_name]] <- chr_stats
  }
  
  res

}
