# Functions specific for mega

# one iteration of worker, use.big.memory = FALSE
.worker <- function(cnvData, cutoff.dup, cutoff.del, Y.star, X, phenoType, use.big.memory) {

  SS <- .scoreProcedure(cnvData, cutoff.dup, cutoff.del, Y.star, X, phenoType, use.big.memory)

  grid <- SS$grid

  TS <- .computeTS(SS$SSA, grid) + .computeTS(SS$SSB, grid)

  # Generate Phi back
  SS <- .combineCovMatrices(SS)

  # Multiply grid length
  Phi <- .combineBlocks(SS, grid)

  ee_values <- .eigen(Phi)
  rm("Phi")

  .getPValues(q = TS, lambda = ee_values, method = "worker")
}

# use.big.memory = TRUE
.mega <- function(meta.data) {

  #meta_grid <- meta.data[[1]]$grid
  meta_grid <- meta.data$grid

  meta.data$meta <- list()

  meta.data$meta$TS <- .computeTS(meta.data$SSA, meta_grid) + .computeTS(meta.data$SSB, meta_grid)
  meta.data <- .mgHxHx.FBM(meta.data, "A")
  meta.data <- .mgHxHy.FBM(meta.data)
  meta.data <- .mgHxHx.FBM(meta.data, "B")
  # meta.data[-which(names(meta.data) == "meta")] <- NULL

  Phi <- .combineBlocks(meta.data$meta, meta_grid)

  ee_meta <- .eigen(Phi)
  rm("Phi")

  .getPValues(q = meta.data$meta$TS, lambda = ee_meta, method = "meta")
}

.mgHxHx.FBM <- function(meta.data, x.name) {

  if (x.name == "A") {
    HH <- "HAHA"
    HX <- "HAX"
  } else {
    HH <- "HBHB"
    HX <- "HBX"
  }
  xx <- bigstatsr::big_tcrossprodSelf(meta.data[[HX]])

  sum_HH <- bigstatsr::FBM(nrow = nrow(xx), ncol = ncol(xx),
                           type = "double", init = 0.0)

  tmp <- bigstatsr::big_apply(
    X = meta.data[[HH]],
    a.FUN = function(X, ind, xx, sum_HH) {
      sum_HH[, ind] <- X[, ind] - xx[, ind]
      NULL
    },
    a.combine = 'c',
    block.size = 100,
    xx = xx,
    sum_HH = sum_HH
  )

  meta.data[[HH]] <- NULL
  meta.data$meta[[HH]] <- sum_HH
  meta.data
}

.mgHxHy.FBM <- function(meta.data) {

  hbt <- bigstatsr::big_transpose(meta.data$HBX)

  sum_HH <- bigstatsr::FBM(nrow = nrow(meta.data$HAHB),
                           ncol = ncol(meta.data$HAHB),
                           type = "double", init = 0.0)

  tmp <- bigstatsr::big_apply(
    X = meta.data$HAHB,
    a.FUN = function(X, ind, hax, hbt, sum_HH) {
      sum_HH[, ind] <- X[, ind] - hax %*% hbt[, ind]
      NULL
    },
    a.combine = 'c',
    block.size = 100,
    hax = meta.data$HAX,
    hbt = hbt,
    sum_HH = sum_HH
  )
  meta.data$HAHB <- NULL

  meta.data$meta$HAHB <- sum_HH
  meta.data
}

#' Internal convenience function to create full covariance matrix
#'
#' @noRd
#' @param cov.result A list object. The result returned by `.scoreProcedure()`
#'
#' @returns The input cov.result with HAHA, HBHB, and HAHB augmented by the
#'   appropriate HXHX. The HAX and HBX terms are removed.
#' @keywords internal
.combineCovMatrices <- function(cov.result) {

  if (inherits(cov.result$HAHA, "FBM")) {

    HXHX <- bigstatsr::big_tcrossprodSelf(cov.result$HAX)
    tmp <- bigstatsr::big_apply(
      X = cov.result$HAHA,
      a.FUN = function(X, ind, hxhx, hh) {
        hh[, ind] <- X[, ind] - hxhx[, ind]
        NULL  # Return NULL to save memory
      },
      a.combine = 'c',
      block.size = 100,
      hxhx = HXHX,
      hh = cov.result$HAHA
    )
    file.remove(HXHX$backingfile)

    HBXt <- bigstatsr::big_transpose(cov.result$HBX)
    tmp <- bigstatsr::big_apply(
      X = cov.result$HAHB,
      a.FUN = function(X, ind, hax, hbxt, hh) {
        hh[, ind] <- X[, ind] - hax %*% hbxt[, ind]
        NULL
      },
      a.combine = 'c',
      block.size = 100,
      hax = cov.result$HAX,
      hbxt = HBXt,
      hh = cov.result$HAHB
    )
    file.remove(HBXt$backingfile)

    HXHX <- bigstatsr::big_tcrossprodSelf(cov.result$HBX)

    tmp <- bigstatsr::big_apply(
      X = cov.result$HBHB,
      a.FUN = function(X, ind, hxhx, hh) {
        hh[, ind] <- X[, ind] - hxhx[, ind]
        NULL
      },
      a.combine = 'c',
      block.size = 100,
      hxhx = HXHX,
      hh = cov.result$HBHB
    )
    file.remove(HXHX$backingfile)

  } else {
    if (!inherits(cov.result$HAHA, "symmetricMatrix"))
      cov.result$HAHA <- as(cov.result$HAHA, "symmetricMatrix")
    if (!inherits(cov.result$HBHB, "symmetricMatrix"))
      cov.result$HBHB <- as(cov.result$HBHB, "symmetricMatrix")
    if (!inherits(cov.result$HAX, "sparseMatrix"))
      cov.result$HAX <- as(cov.result$HAX, "sparseMatrix")
    if (!inherits(cov.result$HBX, "sparseMatrix"))
      cov.result$HBX <- as(cov.result$HBX, "sparseMatrix")

    cov.result$HAHA <- cov.result$HAHA -
      as(tcrossprod(cov.result$HAX), "symmetricMatrix")
    cov.result$HAHB <- cov.result$HAHB - tcrossprod(cov.result$HAX, cov.result$HBX)
    cov.result$HBHB <- cov.result$HBHB -
      as(tcrossprod(cov.result$HBX), "symmetricMatrix")
  }

  cov.result$HAX <- NULL
  cov.result$HBX <- NULL

  cov.result
}

