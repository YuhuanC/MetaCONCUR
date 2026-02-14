# Meta step internal functions

#' Compute Test Statistic TS for Meta
#'
#' @param meta.data A list of lists of summary statistics that contains `$grid`
#' a vector of the grid boundaries; `$SSA` the aligned score statistic for group A;
#' `$SSB` the aligned score statistic for group B; `$HAX` the aligned matrix of HAX;
#' `$HBX` the aligned matrix of HBX.
#' @param meta.grid A vector. The grid boundaries of the combined studies.
#'
#' @return A list that contains `meta` object with scalar `TS`.
#' @include Worker.R
#' @keywords internal

.computeTSMeta <- function(meta.data, meta.grid) {

  meta_obj <- which(names(meta.data) == "meta")
  meta_hold <- meta.data$meta
  meta.data$meta <- NULL
  
  sum_SSA <- lapply(meta.data, "[[", "SSA") |> do.call(what = cbind)
  sum_SSA <- rowSums(sum_SSA)

  sum_SSB <- lapply(meta.data, "[[", "SSB") |> do.call(what = cbind)
  sum_SSB <- rowSums(sum_SSB)
  
  meta.data <- lapply(meta.data, function(x) {x$SSA <- NULL; x$SSB <- NULL; x})
  
  meta_hold$TS <- .computeTS(sum_SSA, meta.grid) + .computeTS(sum_SSB, meta.grid)

  meta.data$meta <- meta_hold
  meta.data

}

# Switch for small Matrix or big memory
.eigen <- function(phi, ...) {
  UseMethod(".eigen", phi)
}

#' Obtain eigenvalues for Matrix Phi
#'
#' @param phi A Matrix. The assembled var-covariance Matrix.
#' @param .force_method Internal: force "base" or "rspectra" for testing only.
#'
#' @return A vector. Eigenvalues
#'
#' @importFrom rlang %||%
#' @import RSpectra
#' @keywords internal
.eigen.Matrix <- function(phi, .force_method = NULL) {

  # Decide which method to use
  method <- .force_method %||% if (nrow(phi) < 9000L) "base" else "rspectra"

  if (method == "base") {
    ee_values <- tryCatch(
      eigen(phi, only.values = TRUE)$values,
      error = function(e) {
        stop("unable to obtain eigenvalues\n\t", e$message, call. = FALSE)
      }
    )
  } else if (method == "rspectra") {
    ee_values <- Inf
    k <- min(500L, nrow(phi) - 1L)

    while (all(abs(ee_values) > 1e-8)) {
      ee_values <- tryCatch(
        RSpectra::eigs_sym(phi, k=k, opts = list(retvec = FALSE))$values,
        error = function(e) {
          stop("unable to obtain eigenvalues using RSpectra\n\t", e$message, call. = FALSE)
        }
      )

#      if (k >= nrow(phi)) break
      k <- min(nrow(phi), k + 500L)

      if(k == nrow(phi)){
        warning("Phi is expected to be sparse, the input Phi is not sparse")

        ee_values <- tryCatch(
          eigen(phi, only.values = TRUE)$values,
          error = function(e) {
            stop("unable to obtain eigenvalues\n\t", e$message, call. = FALSE)
            })
        break
      }

    }
  } else {
    stop("Unknown method: ", method)
  }

  ee_values
}


#' Obtain eigenvalues for FBM Object Phi
#'
#' @param phi An phi FBM object.
#'
#' @return A vector. Eigenvalues
#'
#' @import bigstatsr
#'
#' @keywords internal

.eigen.FBM <- function(phi) {

  phi <- phi$bm()

  n <- nrow(phi)

  wr <- matrix(0.0, nrow = n, ncol = 1L)
  wi <- matrix(0.0, nrow = n, ncol = 1L)
  vl <- bigmemory::big.matrix(1L, n, type = "double", init = 0.0)
  work <- bigmemory::big.matrix(1L, n, type = "double", init = 0.0)
  tryCatch(
    bigalgebra::dgeev(JOBVL = "N", JOBVR = "N",
                      N = n,
                      A = phi,
                      WR = wr, WI = wi,
                      VL = vl, LDVL = 1, LDVR = 1,
                      WORK = work, LWORK = -1),
    error = function(e) {
      stop("unable to determine size of work matrix for eigen decomposition\n\t",
           e$message, call. = FALSE)
    })
  lwork <- work[1L, 1L]

  work <- bigmemory::big.matrix(1L, lwork, type = "double", init = 0.0)
  tryCatch(
    bigalgebra::dgeev(JOBVL = "N", JOBVR = "N",
                      N = n,
                      A = phi,
                      WR = wr, WI = wi,
                      VL = vl, LDVL = 1, LDVR = 1,
                      WORK = work, LWORK = lwork),
    error = function(e) {
      stop("unable to complete eigen decomposition\n\t",
           e$message, call. = FALSE)
    })
  drop(wr)
}

#' Calculated Davies and Liu P-values
#'
#' @noRd
#' @param q A numeric vector. The value at which the distribution function is
#'   to be evaluated. Equivalent of `q` input in CompQuadForm::davies() and
#'   CompQuadForm::liu().
#' @param lambda A numeric vector. The eigenvalues.
#' @param method A character. Specifying the name of method.
#'
#' @return A named numeric vector. Element "davies" contains the probability
#'   as estimated by CompQuadForm::davies(); element "liu" contains the
#'   probability as estimated by CompQuadForm::liu().
#'
#' @importFrom CompQuadForm davies liu
#' @keywords internal

.getPValues <- function(q, lambda, method) {

  stopifnot(
    "`q` must be a numeric vector" = !missing(q) && .isNumericVector(q),
    "`lambda` must be a numeric vector" = !missing(lambda) &&
      .isNumericVector(lambda)
  )

  lambda <- lambda[abs(lambda) >= 1e-8]

  if (length(lambda) == 0L) stop("all eigenvalues are 0", call. = FALSE)

  pv_davies <- tryCatch(CompQuadForm::davies(q, lambda = lambda)$Qq,
                        error = function(e) {
                          stop("unable to complete davies in ", method, "\n\t",
                               e$message, call. = FALSE)
                        })

  pv_liu <- tryCatch(CompQuadForm::liu(q, lambda = lambda),
                     error = function(e) {
                       stop("unable to complete liu in ", method, "\n\t",
                            e$message, call. = FALSE)
                     })

  c("davies" = pv_davies, "liu" = pv_liu)

}

# Worker -> Align -> Meta
#' Input a list of study-specific summary statistics in worker output format,
#' output p-values.
#'
#' @param ... One or more summary statistics lists as generated by CONCUR_Worker
#' each for different studies.
#'
#' @return A named numeric vector. Element "davies" contains the probability
#'   as estimated by CompQuadForm::davies(); element "liu" contains the
#'   probability as estimated by CompQuadForm::liu().
#' @name Meta
#' @include Sum.R Phi.R
#' @examples
#' # Examples for meta with two studies
#' # Load example data
#' data(cnvData)
#' data(phenoData)
#' data(covData)
#' # Run on example data, first 100 individuals for study 1
#' s1_cnv <- cnvData[as.numeric(sub("P", "", cnvData$ID)) >= 1 &
#'                    as.numeric(sub("P", "", cnvData$ID)) <= 100, ]
#' s1_pheno <- phenoData[as.numeric(sub("P", "", phenoData$ID)) >= 1 &
#'                        as.numeric(sub("P", "", phenoData$ID)) <= 100, ]
#' s1_cov <- covData[as.numeric(sub("P", "", covData$ID)) >= 1 &
#'                    as.numeric(sub("P", "", covData$ID)) <= 100, ]
#' # Run on example data, 101-200 individuals for study 2
#' s2_cnv <- cnvData[as.numeric(sub("P", "", cnvData$ID)) >= 101 &
#'                    as.numeric(sub("P", "", cnvData$ID)) <= 200, ]
#' s2_pheno <- phenoData[as.numeric(sub("P", "", phenoData$ID)) >= 101 &
#'                        as.numeric(sub("P", "", phenoData$ID)) <= 200, ]
#' s2_cov <- covData[as.numeric(sub("P", "", covData$ID)) >= 101 &
#'                    as.numeric(sub("P", "", covData$ID)) <= 200, ]
#' cutoff.dup <- 2.3
#' cutoff.del <- 1.7
#' formula_bin <- PHEB ~ SEX
#' # Getting summary statistics for study 1,2
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
#' s2_SS <- CONCUR_worker(
#'   cnvData = s2_cnv,
#'   phenoData = s2_pheno,
#'   cutoff.dup = cutoff.dup,
#'   cutoff.del = cutoff.del,
#'   formula = formula_bin,
#'   covData = s2_cov,
#'   phenoType = "bin",
#'   chr = 1,
#'   use.big.memory = FALSE,
#'   save_dir = NULL
#' )
#' # Running Meta with two studies' SS chr1
#' res_p <- Meta(s1_SS$chr1, s2_SS$chr1)
#'
#' print(res_p)
#' @export

Meta <- function(...) {

  meta.data <- list(...)
  stopifnot(
    "`meta.data` must be a list of lists of summary statistics in Worker-Output form" =
      is.list(meta.data) && 
      all(
        vapply(meta.data, function(x) {
          all(c("DUP", "DEL") %in% names(x)) &&
          is.list(x$DUP) &&
            all(c("grid", "SSA", "SSB", "HAHA", "HAHB", "HBHB", "HAX", "HBX") %in% names(x$DUP)) &&
            is.list(x$DEL) &&
            all(c("grid", "SSA", "SSB", "HAHA", "HAHB", "HBHB", "HAX", "HBX") %in% names(x$DEL))
        }, logical(1))
      ),
    "`meta.data` must contain summary statistics of more than one study" =
      length(meta.data) > 1L
  )
  
  # split the data by DEL/DUP
  meta <- list()
  meta$DEL <- lapply(meta.data, "[[", "DEL")
  meta$DUP <- lapply(meta.data, "[[", "DUP")
  rm(meta.data)

  meta_grid <- lapply(meta$DEL, "[[", "grid") |> unlist() |> unique() |> sort()
  
  for (i in 1L:2L) {

    meta[[i]] <- lapply(meta[[i]], .alignStudyToMeta, meta.grid = meta_grid)

    meta[[i]]$meta <- list()

    meta[[i]] <- .computeTSMeta(meta[[i]], meta_grid)
    meta[[i]] <- .sumHAHA(meta[[i]])
    meta[[i]] <- .sumHAHB(meta[[i]])
    meta[[i]] <- .sumHBHB(meta[[i]])
    meta[[i]] <- .rmStudies(meta[[i]])
  }
  TSmeta <- meta$DEL$meta$TS + meta$DUP$meta$TS
  meta <- .combineDELDUP(meta$DEL, meta$DUP)
  
  Phi <- .combineBlocks(meta$meta, meta_grid)

  ee_meta <- .eigen(Phi)
  rm("Phi")

  .getPValues(q = TSmeta, lambda = ee_meta, method = "meta")
}

#' Load summary statistics (FBM + RDS or all-in-RDS)
#'
#' @param chr_name A character. Chromosome name, e.g. "chr1".
#' @param save_dir A directory path where summary statistics were saved.
#' @param use.big.memory Logical. If TRUE, load FBMs plus small rds;
#'   if FALSE, load full list from SS.rds.
#'
#' @return A list identical in structure to the output of Worker
#' @describeIn Meta A convenient function to retrieve Worker result
#' @export
load_chr_stats <- function(chr_name, save_dir, use.big.memory = FALSE) {
  stopifnot(
    "`save_dir` must be provided" =
      !is.null(save_dir),
    "`chr_name` must be specified" =
      !missing(chr_name)
    )
  chr_dir <- file.path(save_dir, chr_name)

  if (!dir.exists(chr_dir)) {
    stop(paste0("Chromosome directory not found: ", chr_dir))
  }

  # get file directory dir function test if wo_FBM.rds exists, message

  fbm_names <- c("HAHA", "HAHB", "HBHB", "HAX", "HBX")

  # List all files in the directory
  files_in_dir <- dir(chr_dir)

  # Append .rds extension
  fbm_files <- paste0(fbm_names, ".rds")

  # Find which are present
  found <- fbm_files[fbm_files %in% files_in_dir]
  missing <- setdiff(fbm_files, found)

  # if false, found FBMs > 0, message... sure about FBM false?
  # if true, print messages, specific names of FBM files, wo_FBM..

  if (isTRUE(use.big.memory)) {
    # Case 1: big.memory = TRUE

    # check if SS small file exists
    file_name <- "wo_FBM.rds"
    if (file_name %in% files_in_dir) {
      message("Found Small Statistics Bundle: ", file_name, " in ", chr_dir)
    } else {
      stop(paste0("Missing Small Statistics Bundle: ", file_name, " in ", chr_dir))
    }

    # Load the small objects bundle
    small_stats <- readRDS(file.path(chr_dir, file_name))

    # check if FBMs exist
    if (length(missing) > 0) {
      stop(paste0("Missing FBM files in ", chr_dir, ":\n  - ",
                  paste(missing, collapse = "\n  - ")))
    }else{
      message("Found FBM Files: ", paste(found, collapse = ", "))
    }

    # Reattach FBMs from their .rds metadata
    fbm_stats <- lapply(fbm_names, function(nm) {
      bigstatsr::big_attach(file.path(chr_dir, paste0(nm, ".rds")))
    })
    names(fbm_stats) <- fbm_names

    stats <- c(small_stats, fbm_stats)

  } else {
    # Case 2: big.memory = FALSE

    # remind if they want to use big memory when FBMs exist
    if(length(found) > 0){
      message("FBM Files are present but will be ignored (use.big.memory = FALSE).")
    }

    # check if SS exist
    file_name2 <- "SS.rds"
    if (file_name2 %in% files_in_dir) {
      message("Found Summary Statistics (Matrices Version): ", file_name2, " in ", chr_dir)
    } else {
      stop(paste0("Missing Summary Statistics (Matrices Version): ", file_name2, " in ", chr_dir))
    }
    stats <- readRDS(file.path(chr_dir, file_name2))
  }

  stats
}
