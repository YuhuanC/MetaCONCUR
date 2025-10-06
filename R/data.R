#' Example cnvData
#'
#' Simulated CNV data in PLINK CNV data format for 600 individuals.
#'
#' @format A data frame containing 803 observations with 5 columns:
#' ID character patient identifier.
#' CHR CNV chromosome.
#' BP1 starting location in base pairs.
#' BP2 ending location in base pairs.
#' TYPE copy number (0,1,3,or 4).
"cnvData"

#' Example phenotype phenoData
#'
#' Simulated phenotype data, include a binary phenotype for 600 individuals.
#'
#' @format A data frame containing 600 observations with 2 columns:
#' ID character patient identifier.
#' PHEB binary phenotype.
"phenoData"

#' Example covariates covData
#'
#' Simulated covariate data, include SEX covariate for 600 individuals.
#'
#' @format A data frame containing 600 observations with 2 columns:
#' ID character patient identifier.
#' SEX binary indicator of M/F.
"covData"
