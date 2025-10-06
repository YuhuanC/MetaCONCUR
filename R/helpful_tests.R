#' Test if provided object is a numeric vector of specified length
#'
#' @noRd
#' @param x An R object.
#' @param length An integer or NULL. If NULL, length is tested to be non-zero;
#'   if positive, length of x must match.
#'
#' @returns A logical. TRUE if x is a numeric vector of appropriate length
#' @keywords internal
.isNumericVector <- function(x, length = NULL) {
  tst <- is.numeric(x) && is.vector(x, mode = "numeric") && length(x) > 0L
  if (!is.null(length)) tst <- tst && length(x) == length
  tst
}

#' Test if provided object is a character vector of specified length
#'
#' @noRd
#' @param x An R object.
#' @param length An integer or NULL. If NULL, length is tested to be non-zero;
#'   if positive, length of x must match.
#'
#' @returns A logical. TRUE if x is a character vector of appropriate length.
#' @keywords internal
.isCharacterVector <- function(x, length = NULL) {
  if (is.null(length)) {
    is.character(x) && is.vector(x, mode = "character") && length(x) > 0L
  } else {
    is.character(x) && is.vector(x, mode = "character") && length(x) == length
  }
}

#' Test if provided object is an integer vector of specified length
#'
#' @noRd
#' @param x An R object.
#' @param length An integer or NULL. If NULL, length is tested to be non-zero;
#'   if positive, length of x must match.
#'
#' @returns A logical. TRUE if x is an integer vector of appropriate length.
#' @keywords internal
.isIntegerVector <- function(x, length = NULL) {
  if (is.null(length)) {
    is.integer(x) && is.vector(x, mode = "integer") && length(x) > 0L
  } else {
    is.integer(x) && is.vector(x, mode = "integer") && length(x) == length
  }
}

#' Test if provided object is a logical vector of specified length
#'
#' @noRd
#' @param x An R object.
#' @param length An integer or NULL. If NULL, length is tested to be non-zero;
#'   if positive, length of x must match.
#'
#' @returns A logical. TRUE if x is a logical vector of appropriate length.
#' @keywords internal
.isLogicalVector <- function(x, length = NULL) {
  if (is.null(length)) {
    is.logical(x) && is.vector(x, mode = "logical") && length(x) > 0L
  } else {
    is.logical(x) && is.vector(x, mode = "logical") && length(x) == length
  }
}

#' Test if provided object is a list containing at a minimum elements nms
#'
#' @noRd
#' @param x An R object.
#' @param nms A character vector of element names
#'
#' @returns A logical. TRUE if x is a list containing the nms elements
#' @keywords internal
.isNamedList <- function(x, nms) {
  tst <- is.vector(x, mode = "list") && length(x) >= length(nms)
  tst <- tst && all(nms %in% names(x))
  tst
}

.isCNV <- function(CNV) {
  
  if (missing(CNV)) stop("`CNV` must be provided", call. = FALSE)
  
  if (!is.data.frame(CNV) || !all(c("ID", "CHR", "BP1", "BP2", "TYPE") %in% colnames(CNV))) {
    stop("`CNV` data.frame in PLINK format with columns: `ID`, `CHR`, `BP1`, `BP2`,`TYPE`",
         call. = FALSE)
  }
  
#  if (!is.vector(CNV$CHR, mode = "integer") ||
#      length(unique(CNV$CHR)) != 1L || !all(as.character(CNV$CHR) %in% as.character(c(1L:22L)))) {
#    stop("`CNV$CHR` must take a single value in [1, 22]", call. = FALSE)
#  }
  
  if (!is.vector(CNV$BP1, mode = "numeric") || !all(CNV$BP1 >= 0L) ||
      !isTRUE(all.equal(CNV$BP1, round(CNV$BP1)))) {
    stop("`CNV$BP1 must contain non-negative integers", call. = FALSE)
  }
  
  if (!is.vector(CNV$BP2, mode = "numeric") || !all(CNV$BP2 >= 0L) ||
      !isTRUE(all.equal(CNV$BP2, round(CNV$BP2)))) {
    stop("`CNV$BP2 must contain non-negative integers", call. = FALSE)
  }
  
  if (any(CNV$BP1 > CNV$BP2)) {
    stop("`CNV$BP2 must be great than CNV$BP1 for all records", call. = FALSE)
  }
  
  if (!is.vector(CNV$TYPE, mode = "integer") ||
      !all(CNV$TYPE %in% c(0L, 1L, 3L:9999L)) ||
      sum(CNV$TYPE == 2L) != 0L) {
    stop("`CNV$Type` must be integer in [0, 1, 3, ..., 9999]", call. = FALSE)
  }
  
  TRUE
  
}

