# Multiply grid length and combine blocks to get Phi for eigen-decomp

#' Use grid length as weight to update covariance Matrix Phi
#'
#' @param blocks A `meta` object from `meta.data` list.
#' @param grid A vector. The grid boundaries of the combined studies.
#'
#' @return A sparse Matrix. The combined Phi variance-covariance Matrix.
#' @import Matrix
#'
#' @keywords internal

.combineBlocksSmall <- function(blocks, grid) {

  sG_H <- rep(diff(grid), times = 2L) |> sqrt()

  colN <- diff(blocks$HAHA@p)
  weight <- sG_H[blocks$HAHA@i + 1L] * sG_H[rep(seq_along(colN), colN)]
  blocks$HAHA@x <- blocks$HAHA@x * weight

  colN <- diff(blocks$HAHB@p)
  weight <- sG_H[blocks$HAHB@i + 1L] * sG_H[rep(seq_along(colN), colN)]
  blocks$HAHB@x <- blocks$HAHB@x * weight

  colN <- diff(blocks$HBHB@p)
  weight <- sG_H[blocks$HBHB@i + 1L] * sG_H[rep(seq_along(colN), colN)]
  blocks$HBHB@x <- blocks$HBHB@x * weight

  nonzeros <- list()
  nonzeros$HAHA <- apply(blocks$HAHA, 2, function(x) any(x > 1e-15 | x < -1e-15))
  nonzeros$HAHB <- apply(blocks$HAHB, 2, function(x) any(x > 1e-15 | x < -1e-15))
  nonzeros$HBHA <- apply(blocks$HAHB, 1, function(x) any(x > 1e-15 | x < -1e-15))
  nonzeros$HBHB <- apply(blocks$HBHB, 2, function(x) any(x > 1e-15 | x < -1e-15))

  # handle the columns first:
  cols_non_zero_left <- nonzeros$HAHA | nonzeros$HBHA #(left side of phi)
  cols_non_zero_right <- nonzeros$HAHB | nonzeros$HBHB #(right side of phi)

  phi_top <- cbind(blocks$HAHA[cols_non_zero_left, cols_non_zero_left, drop=FALSE],
                   blocks$HAHB[cols_non_zero_left, cols_non_zero_right, drop=FALSE])
  phi_bottom <- cbind(Matrix::t(blocks$HAHB)[cols_non_zero_right, cols_non_zero_left, drop=FALSE],
                      blocks$HBHB[cols_non_zero_right, cols_non_zero_right, drop=FALSE])
  phi <- rbind(phi_top, phi_bottom)

  phi
}

#' Use grid length as weight to update FBM Phi
#'
#' @param blocks A `meta` object from `meta.data` list.
#' @param grid A vector. The grid boundaries of the combined studies.
#'
#' @return A FBM object. The combined Phi FBM.
#' @import bigstatsr
#'
#' @keywords internal
#'
.combineBlocksBig <- function(blocks, grid) {
  # multiply all blocks by the grid unit size

  sG_H <- rep(diff(grid), times = 2L) |> sqrt()

  for (block in c("HAHA", "HAHB", "HBHB")) {
    tmp <- bigstatsr::big_apply(
      X = blocks[[block]],
      a.FUN = function(X, ind, gg, hh) {
        hh[, ind] <- X[, ind] %*% diag(gg[ind]) * gg
      },
      a.combine = 'c',
      block.size = 100,
      gg = sG_H,
      hh = blocks[[block]]
    )
  }

  blocks$HBHA <- bigstatsr::big_transpose(blocks$HAHB)

  # identify columns with no non-zero values
  nonzeros <- list()
  for (block in c("HAHA", "HAHB", "HBHA", "HBHB")) {
    nonzeros[[block]] <- bigstatsr::big_apply(
      X = blocks[[block]],
      a.FUN = function(X, ind) {
        colSums(abs(X[, ind])) > 1e-15
      },
      a.combine = 'c',
      block.size = 100
    )
  }

  # handle the columns first:
  cols_non_zero_left <- {nonzeros$HAHA} | {nonzeros$HBHA} #(left side of phi)
  cols_non_zero_right <- {nonzeros$HAHB} | {nonzeros$HBHB} #(right side of phi)

  HAHA_small <- bigstatsr::big_copy(blocks$HAHA, which(cols_non_zero_left), which(cols_non_zero_left))
  file.remove(blocks$HAHA$backingfile)
  HAHB_small <- bigstatsr::big_copy(blocks$HAHB, which(cols_non_zero_left), which(cols_non_zero_right))
  file.remove(blocks$HAHB$backingfile)
  HBHB_small <- bigstatsr::big_copy(blocks$HBHB, which(cols_non_zero_right), which(cols_non_zero_right))
  file.remove(blocks$HBHB$backingfile)
  HBHA_small <- bigstatsr::big_copy(blocks$HBHA, which(cols_non_zero_right), which(cols_non_zero_left))
  file.remove(blocks$HBHA$backingfile)

  phi <- bigstatsr::FBM(nrow = nrow(HAHA_small) + nrow(HBHA_small),
                        ncol = ncol(HAHA_small) + ncol(HAHB_small),
                        type = "double", init = 0.0)

  offset <- nrow(HAHA_small)
  offsetc <- ncol(HAHA_small)

  tmp <- bigstatsr::big_apply(
    X = HAHA_small,
    a.FUN = function(X, ind, phi) {
      phi[seq_len(nrow(X)), ind] <- X[, ind]
      NULL
    },
    a.combine = 'c',
    block.size = 100,
    phi = phi
  )
  file.remove(HAHA_small$backingfile)

  tmp <- bigstatsr::big_apply(
    X = HAHB_small,
    a.FUN = function(X, ind, phi, offsetc) {
      phi[seq_len(nrow(X)), ind + offsetc] <- X[, ind]
      NULL
    },
    a.combine = 'c',
    block.size = 100,
    phi = phi,
    offsetc = offsetc
  )
  file.remove(HAHB_small$backingfile)

  tmp <- bigstatsr::big_apply(
    X = HBHA_small,
    a.FUN = function(X, ind, phi, offset) {
      phi[seq_len(nrow(X)) + offset, ind] <- X[, ind]
      NULL
    },
    a.combine = 'c',
    block.size = 100,
    phi = phi,
    offset = offset
  )
  file.remove(HBHA_small$backingfile)

  tmp <- bigstatsr::big_apply(
    X = HBHB_small,
    a.FUN = function(X, ind, phi, offsetc) {
      phi[seq_len(nrow(X)) + offset, ind + offsetc] <- X[, ind]
      NULL
    },
    a.combine = 'c',
    block.size = 100,
    phi = phi,
    offsetc = offsetc
  )
  file.remove(HBHB_small$backingfile)

  phi
}

# Switch for small Matrix Phi or big memory Phi FBM
.combineBlocks <- function(blocks, grid) {

  stopifnot(
    "`blocks` must be a list of lists of covariance matrix components" =
      !missing(blocks) && is.list(blocks) &&
      all(c("HAHA", "HAHB", "HBHB") %in% names(blocks)),
    "`grid` must be a numeric vector" =
      !missing(grid) && .isNumericVector(grid) && length(grid) > 1L
  )

  if (inherits(blocks$HAHA, "FBM")) {
    .combineBlocksBig(blocks, grid)
  } else {
    .combineBlocksSmall(blocks, grid)
  }

}
