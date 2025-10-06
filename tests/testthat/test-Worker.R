test_that(".idDuplicates works", {
  test_cnv <- data.frame("TYPE" = c(4,1,1))
  test_all <- within(test_cnv, TYPE <- c(4, 4, 3))
  test_none <- within(test_cnv, TYPE <- c(1, 1, 1))
  test_one <- data.frame(TYPE = 3)

  test_df <- data.frame(ID = 1:5)

  expect_identical(.idDuplicates(test_cnv,2), c(TRUE,FALSE,FALSE))
  expect_warning(.idDuplicates(test_all, 2), "all records code as duplication")
  expect_warning(.idDuplicates(test_none, 2), "no records code as duplication")
  expect_warning(
    expect_identical(.idDuplicates(test_one,2), TRUE),
    "all records code as duplication"
  )

  expect_error(.idDuplicates(NULL,2), "`cnvData` must be a data.frame")
  expect_error(.idDuplicates(test_df,2), "`cnvData` must contain column `TYPE`")
  expect_error(.idDuplicates(test_cnv,NULL), "`cutoff` must be numeric")

})

test_that(".idDeletions works", {
  test_cnv <- data.frame("TYPE" = c(4,1,1))
  test_all <- within(test_cnv, TYPE <- c(1, 0, 1))
  test_none <- within(test_cnv, TYPE <- c(4, 3, 3))
  test_one <- data.frame(TYPE = 1)

  test_df <- data.frame(ID = 1:5)

  expect_identical(.idDeletions(test_cnv,2), c(FALSE,TRUE,TRUE))
  expect_warning(.idDeletions(test_all, 2), "all records code as deletion")
  expect_warning(.idDeletions(test_none, 2), "no records code as deletion")
  expect_warning(
    expect_identical(.idDeletions(test_one,2), TRUE),
    "all records code as deletion"
  )

  expect_error(.idDeletions(NULL,2), "`cnvData` must be a data.frame")
  expect_error(.idDeletions(test_df,2), "`cnvData` must contain column `TYPE`")
  expect_error(.idDeletions(test_cnv,NULL), "`cutoff` must be numeric")

})

test_that(".defineGrid works",{
  test_cnv <- data.frame("ID" = c("P1","P2","P3"),
                         "CHR" = c(1,1,1),
                         "BP1" = c(150,200,200),
                         "BP2" = c(200,300,300),
                         "TYPE" = c(4L,1L,1L))

  test_cnv2 <- data.frame("TYPE" = c(4,1,1))

  expect_identical(
    .defineGrid(test_cnv),
    sort(unique(c(test_cnv$BP1, test_cnv$BP2)))
  )

  expect_error(.defineGrid())
  expect_error(.defineGrid(test_cnv2))
})

test_that(".matchGrid works",{
  test_cnv <- data.frame("ID" = c("P1","P2","P3"),
                         "CHR" = c(1,1,1),
                         "BP1" = c(150,200,200),
                         "BP2" = c(200,300,300),
                         "TYPE" = c(4L,1L,1L))
  grid <- .defineGrid(test_cnv)

  expected <- as.data.frame(list("ID" = test_cnv$ID,
                            "TYPE" = test_cnv$TYPE,
                            "lower" = as.integer(findInterval(test_cnv$BP1, grid)),
                            "upper" = as.integer(findInterval(test_cnv$BP2, grid,
                                                              left.open = TRUE))))

  expect_equal(.matchGrid(test_cnv,grid),expected)

  expect_error(.matchGrid(,grid))
  expect_error(.matchGrid(NULL,grid))
  expect_error(.matchGrid(test_cnv,NULL),"`grid` must be a numeric vector")
})

test_that(".getYstar works",{

  withr::with_seed(12345,
                   {n <- 10
                   covData <- data.frame(
                     ID = paste0("P", 1:n),
                     X1 = runif(n, 0, 1),
                     X2 = sample(0:1, n, replace = TRUE)
                   )

                   z <- 0.5 * covData$X1 + 0.3 * covData$X2 + rnorm(n)

                   # Simulate phenotype data
                   Y <- data.frame(
                     ID = paste0("P", 1:n),
                     Yc = z,
                     Yb = rbinom(n, size = 1, prob = 1 / (1 + exp(-z)))
                   )})

  Y2 <- Y[-1,]

  formula_c <- Yc ~ X1 + X2
  formula_b <- Yb ~ X1 + X2

  id_order <- order(covData$ID)
  covD <- covData[id_order, ]
  phenoD <- Y[id_order, ]

  model_cont <- stats::lm(formula_c, cbind(covD, phenoD))

  model_bin <- stats::glm(formula_b, cbind(covD, phenoD), family = "binomial")

  expected_cont <- as.data.frame(list("ID" = covD$ID,
                                      "Residuals" = stats::residuals(model_cont,
                                                                    type = "response"),
                                      "Fitted" = stats::fitted.values(model_cont)))

  expected_bin <- as.data.frame(list("ID" = covD$ID,
                                     "Residuals" = stats::residuals(model_bin,
                                                                    type = "response"),
                                     "Fitted" = stats::fitted.values(model_bin)))

  expect_equal(.getYstar(covData, Y, formula_c, "cont"), expected = expected_cont)

  expect_equal(.getYstar(covData, Y, formula_b, "bin"), expected = expected_bin)

  expect_error(.getYstar(NULL, Y, formula_b, "bin"), "`covData` must be a data.frame")
  expect_error(.getYstar(covData, NULL, formula_b, "bin"), "`phenoData` must be a data.frame")
  expect_error(.getYstar(covData, Y, NULL, "bin"),
               "`formula` must be a formula of the form LHS ~ RHS")
  expect_error(.getYstar(covData, Y2, formula_b, "bin"),
               "`covData` and `phenoData` must contain the same IDs")
})

test_that(".computeSS works", {

  withr::with_seed(12345,
                   {n <- 10
                   covData <- data.frame(
                     ID = paste0("P", 1:n),
                     X1 = runif(n, 0, 1),
                     X2 = sample(0:1, n, replace = TRUE)
                   )

                   z <- 0.5 * covData$X1 + 0.3 * covData$X2 + rnorm(n)

                   # Simulate phenotype data
                   Y <- data.frame(
                     ID = paste0("P", 1:n),
                     Yc = z,
                     Yb = rbinom(n, size = 1, prob = 1 / (1 + exp(-z)))
                   )})

  formula_c <- Yc ~ X1 + X2
  formula_b <- Yb ~ X1 + X2

  y_cont <- .getYstar(covData, Y, formula_c, "cont")
  y_bin <- .getYstar(covData, Y, formula_b, "bin")

  test_cnv <- data.frame("ID" = c("P1","P2","P3"),
                         "CHR" = c(1,1,1),
                         "BP1" = c(150,200,200),
                         "BP2" = c(200,300,300),
                         "TYPE" = c(4L,1L,1L))

  func.dup <- function(x, cutoff.dup) { as.integer(x > {cutoff.dup - 1e-8}) }
  func.del <- function(x, cutoff.del) { as.integer(x < {cutoff.del + 1e-8}) }

  grid <- .defineGrid(test_cnv)
  M <- .matchGrid(test_cnv, grid)
  M2 <- M
  M2$ID <- c("P10","P20","P30")

  X <- stats::model.matrix(stats::update.formula(formula_b, NULL ~ .),
                           covData[order(covData$ID),])
  px <- ncol(X)
  n_bin <- nrow(y_bin)
  n_cont <- nrow(y_cont)

  phi_cont <-drop(crossprod(y_cont$Residuals)/(n_cont-px))
  phi_bin <- 1.0

  SS_cont <- rep(0, max(M$upper))
  SS_bin <- rep(0, max(M$upper))

  # Expected value for continuous, deletions
  for(i in seq_len(n_cont)){
    M_i <- M[M$ID==y_cont$ID[i],]
    SS_i <- rep(0, max(M$upper))
    if(nrow(M_i) == 0){
      SS_i <- rep(0, max(M$upper))
    } else{
      for (j in seq_len(nrow(M_i))) {
        if(func.del(M_i$TYPE[j],1.7)){
          indx <- M_i$lower[j]:M_i$upper[j]
          SS_i[indx] <- SS_i[indx] + y_cont$Residuals[i]/phi_cont
        }
      }
    }
    SS_cont <- SS_cont + SS_i
  }

  # Expected value for continuous, duplications
  SS_cont_dup <- rep(0, max(M$upper))
  for(i in seq_len(n_cont)){
    M_i <- M[M$ID==y_cont$ID[i],]
    SS_i <- rep(0, max(M$upper))
    if(nrow(M_i)==0){
      SS_i <- rep(0, max(M$upper))
    } else{
      for (j in seq_len(nrow(M_i))) {
        if(func.dup(M_i$TYPE[j],2.3)){
          indx <- M_i$lower[j]:M_i$upper[j]
          SS_i[indx] <- SS_i[indx] + y_cont$Residuals[i]/phi_cont
        }
      }
    }
    SS_cont_dup <- SS_cont_dup + SS_i
  }

  # Expected value for binary, deletions
  for(i in seq_len(n_bin)){
    M_i <- M[M$ID == y_bin$ID[i], ]
    SS_i <- rep(0, max(M$upper))
    if(nrow(M_i)==0){
      SS_i <- rep(0, max(M$upper))
    } else{
      for (j in seq_len(nrow(M_i))) {
        if(func.del(M_i$TYPE[j],1.7)){
          indx <- M_i$lower[j]:M_i$upper[j]
          SS_i[indx] <- SS_i[indx] + y_bin$Residuals[i]/phi_bin
        }
      }
    }
    SS_bin <- SS_bin + SS_i
  }

  # Expected value for binary, duplications
  SS_bin_dup <- rep(0, max(M$upper))
  for(i in seq_len(n_bin)){
    M_i <- M[M$ID == y_bin$ID[i], ]
    SS_i <- rep(0, max(M$upper))
    if(nrow(M_i)==0){
      SS_i <- rep(0, max(M$upper))
    } else{
      for (j in seq_len(nrow(M_i))) {
        if(func.dup(M_i$TYPE[j],2.3)){
          indx <- M_i$lower[j]:M_i$upper[j]
          SS_i[indx] <- SS_i[indx] + y_bin$Residuals[i]/phi_bin
        }
      }
    }
    SS_bin_dup <- SS_bin_dup + SS_i
  }

  expect_equal(.computeSS(y_cont,M,func.del,1.7,px,"cont"),SS_cont)
  expect_equal(.computeSS(y_bin,M,func.del,1.7,px,"bin"),SS_bin)

  expect_equal(.computeSS(y_cont,M,func.dup,2.3,px,"cont"),SS_cont_dup)
  expect_equal(.computeSS(y_bin,M,func.dup,2.3,px,"bin"),SS_bin_dup)

  expect_error(.computeSS(NULL,M,func.del,1.7,px,"cont"),
               "`Y` must be a data.frame with columns ID, Residuals, Fitted" )
  expect_error(.computeSS(y_cont,NULL,func.del,1.7,px,"cont"),
               "`grid.ind` must be a data.frame with columns lower, upper, TYPE, and ID")
  expect_error(.computeSS(y_cont,M,NULL,1.7,px,"cont"),
               "`test.func` must be a function with inputs `x` and `cutoff`")
  expect_error(.computeSS(y_cont,M,func.del,NULL,px,"cont"),
               "`cutoff` must be a scalar")
  expect_error(.computeSS(y_cont,M2,func.del,1.7,px,"cont"),
               "IDs of `grid.ind` must be in IDs of `Y`")

})

test_that(".computeXTXinvsqrt works", {

  withr::with_seed(12345,
                   {n <- 10
                   covData <- data.frame(
                     ID = paste0("P", 1:n),
                     X1 = runif(n, 0, 1),
                     X2 = sample(0:1, n, replace = TRUE)
                   )

                   z <- 0.5 * covData$X1 + 0.3 * covData$X2 + rnorm(n)

                   # Simulate phenotype data
                   Y <- data.frame(
                     ID = paste0("P", 1:n),
                     Yc = z,
                     Yb = rbinom(n, size = 1, prob = 1 / (1 + exp(-z)))
                   )})

  formula_c <- Yc ~ X1 + X2
  formula_b <- Yb ~ X1 + X2

  y_cont <- .getYstar(covData, Y, formula_c, "cont")
  y_bin <- .getYstar(covData, Y, formula_b, "bin")

  X <- stats::model.matrix(stats::update.formula(formula_b, NULL ~ .),
                           covData[order(covData$ID),])
  px <- ncol(X)

  phi.cont <- drop(crossprod(y_cont$Residuals)/(nrow(y_cont)-px))
  phi.hat.cont <- rep(phi.cont, nrow(y_cont))
  phi.hat.bin <- 1.0 / {y_bin$Fitted * {1.0 - y_bin$Fitted}}

  expected_cont <- expm::sqrtm(solve(crossprod(X, X / phi.hat.cont)))
  expected_bin <- expm::sqrtm(solve(crossprod(X, X / phi.hat.bin)))

  expect_equal(.computeXTXinvsqrt(X,phi.hat.cont), expected_cont)
  expect_equal(.computeXTXinvsqrt(X,phi.hat.bin), expected_bin)

  expect_error(.computeXTXinvsqrt(NULL, phi.hat.cont), "`X` must be a matrix")
  expect_error(.computeXTXinvsqrt(X, NULL), "`phi.hat` must be a numeric vector")
})

test_that(".computeCov.Matrix works", {
  library(Matrix)

  withr::with_seed(12345,
                   {n <- 10
                   covData <- data.frame(
                     ID = paste0("P", 1:n),
                     X1 = runif(n, 0, 1),
                     X2 = sample(0:1, n, replace = TRUE)
                   )

                   z <- 0.5 * covData$X1 + 0.3 * covData$X2 + rnorm(n)

                   # Simulate phenotype data
                   Y <- data.frame(
                     ID = paste0("P", 1:n),
                     Yc = z,
                     Yb = rbinom(n, size = 1, prob = 1 / (1 + exp(-z)))
                   )})

  formula_c <- Yc ~ X1 + X2
  formula_b <- Yb ~ X1 + X2

  y_cont <- .getYstar(covData, Y, formula_c, "cont")
  y_bin <- .getYstar(covData, Y, formula_b, "bin")

  X <- stats::model.matrix(stats::update.formula(formula_b, NULL ~ .),
                           covData[order(covData$ID),])
  px <- ncol(X)

  phi.cont <- drop(crossprod(y_cont$Residuals)/(nrow(y_cont)-px))
  phi.hat.cont <- rep(phi.cont, nrow(y_cont))
  phi.hat.bin <- 1.0 / {y_bin$Fitted * {1.0 - y_bin$Fitted}}

  xtx_cont <- .computeXTXinvsqrt(X,phi.hat.cont)
  xtx_bin <- .computeXTXinvsqrt(X,phi.hat.bin)

  ids_cont <- y_cont$ID
  ids_bin <- y_bin$ID

  test_cnv <- data.frame("ID" = c("P1","P2","P3"),
                         "CHR" = c(1,1,1),
                         "BP1" = c(150,200,200),
                         "BP2" = c(200,300,300),
                         "TYPE" = c(4L,1L,1L))

  func.dup <- function(x, cutoff.dup) { as.integer(x > {cutoff.dup - 1e-8}) }
  func.del <- function(x, cutoff.del) { as.integer(x < {cutoff.del + 1e-8}) }

  grid <- .defineGrid(test_cnv)
  M <- .matchGrid(test_cnv, grid)

  HA <- c(do.call(func.dup, list(M$TYPE, 2.3)),
          do.call(func.del, list(M$TYPE, 1.7)))
  HB <- c(do.call(func.dup, list(M$TYPE, 2.3 + 1.0)),
          do.call(func.del, list(M$TYPE, 1.7 - 1.0)))

  L <- max(M$upper)

  L2 <- 2L * L
  p <- ncol(X)

  # initialize covariance components matrices
  HAX <- Matrix::Matrix(0.0, nrow = L2, ncol = p, sparse = TRUE)
  HBX <- Matrix::Matrix(0.0, nrow = L2, ncol = p, sparse = TRUE)
  HAHA <- Matrix::Matrix(0.0, nrow = L2, ncol = L2, sparse = TRUE)
  HAHB <- Matrix::Matrix(0.0, nrow = L2, ncol = L2, sparse = TRUE)
  HBHB <- Matrix::Matrix(0.0, nrow = L2, ncol = L2, sparse = TRUE)

  # unique ID in the event data
  unique_ids <- unique(M$ID)
  grid.ind <- rbind(M,
                    data.frame("ID" = M$ID,
                               "TYPE" = M$TYPE,
                               "lower" = M$lower + L,
                               "upper" = M$upper + L))

  for (i in unique_ids) {
    # rows corresponding to this ID
    use <- ids_cont == i

    # find A and B cases for this ID
    cases_A <- which(grid.ind$ID == i & HA == 1L)
    cases_B <- which(grid.ind$ID == i & HB == 1L)

    if (length(cases_A) > 0) {
      for (j in cases_A) {
        for (idx in grid.ind$lower[j]:grid.ind$upper[j]) {
          # HAHA diagonal block
          HAHA[idx, idx] <- HAHA[idx, idx] + 1.0 / phi.hat.cont[use]

          # HAX row
          for (pp in 1:p) {
            HAX[idx, pp] <- HAX[idx, pp] + X[use, pp] / phi.hat.cont[use]
          }
        }
      }

      if (length(cases_B) > 0) {
        for (jA in cases_A) {
          for (idxA in grid.ind$lower[jA]:grid.ind$upper[jA]) {
            for (jB in cases_B) {
              for (idxB in grid.ind$lower[jB]:grid.ind$upper[jB]) {
                # HAHB cross-block
                HAHB[idxA, idxB] <- HAHB[idxA, idxB] + 1.0 / phi.hat.cont[use]

                # HBHB diagonal block
                HBHB[idxB, idxB] <- HBHB[idxB, idxB] + 1.0 / phi.hat.cont[use]

                # HBX row
                for (pp in 1:p) {
                  HBX[idxB, pp] <- HBX[idxB, pp] + X[use, pp] / phi.hat.cont[use]
                }
              }
            }
          }
        }
      }
    } else if (length(cases_B) > 0) {
      for (jB in cases_B) {
        for (idxB in grid.ind$lower[jB]:grid.ind$upper[jB]) {
          # HBHB diagonal
          HBHB[idxB, idxB] <- HBHB[idxB, idxB] + 1.0 / phi.hat.cont[use]

          # HBX row
          for (pp in 1:p) {
            HBX[idxB, pp] <- HBX[idxB, pp] + X[use, pp] / phi.hat.cont[use]
          }
        }
      }
    }
  }


  expected_cont <- list("HAHA" = as(HAHA, "symmetricMatrix"),
                        "HAHB" = HAHB,
                        "HBHB" = as(HBHB, "symmetricMatrix"),
                        "HAX" = HAX %*% xtx_cont,
                        "HBX" = HBX %*% xtx_cont)

  expect_equal(.computeCov.Matrix(y_cont, X, xtx_cont, phi.hat.cont, M, HA, HB),
               expected_cont)

  # reset for binary
  # dense matrices, initialized to zero
  HAX <- Matrix::Matrix(0.0, nrow = L2, ncol = p, sparse = TRUE)
  HBX <- Matrix::Matrix(0.0, nrow = L2, ncol = p, sparse = TRUE)
  HAHA <- Matrix::Matrix(0.0, nrow = L2, ncol = L2, sparse = TRUE)
  HAHB <- Matrix::Matrix(0.0, nrow = L2, ncol = L2, sparse = TRUE)
  HBHB <- Matrix::Matrix(0.0, nrow = L2, ncol = L2, sparse = TRUE)

  for (i in unique_ids) {
    # rows corresponding to this ID
    use <- ids_bin == i

    # find A and B cases for this ID
    cases_A <- which(grid.ind$ID == i & HA == 1L)
    cases_B <- which(grid.ind$ID == i & HB == 1L)

    if (length(cases_A) > 0) {
      for (j in cases_A) {
        for (idx in grid.ind$lower[j]:grid.ind$upper[j]) {
          # HAHA diagonal block
          HAHA[idx, idx] <- HAHA[idx, idx] + 1.0 / phi.hat.bin[use]

          # HAX row
          for (pp in 1:p) {
            HAX[idx, pp] <- HAX[idx, pp] + X[use, pp] / phi.hat.bin[use]
          }
        }
      }

      if (length(cases_B) > 0) {
        for (jA in cases_A) {
          for (idxA in grid.ind$lower[jA]:grid.ind$upper[jA]) {
            for (jB in cases_B) {
              for (idxB in grid.ind$lower[jB]:grid.ind$upper[jB]) {
                # HAHB cross-block
                HAHB[idxA, idxB] <- HAHB[idxA, idxB] + 1.0 / phi.hat.bin[use]

                # HBHB diagonal block
                HBHB[idxB, idxB] <- HBHB[idxB, idxB] + 1.0 / phi.hat.bin[use]

                # HBX row
                for (pp in 1:p) {
                  HBX[idxB, pp] <- HBX[idxB, pp] + X[use, pp] / phi.hat.bin[use]
                }
              }
            }
          }
        }
      }
    } else if (length(cases_B) > 0) {
      for (jB in cases_B) {
        for (idxB in grid.ind$lower[jB]:grid.ind$upper[jB]) {
          # HBHB diagonal
          HBHB[idxB, idxB] <- HBHB[idxB, idxB] + 1.0 / phi.hat.bin[use]

          # HBX row
          for (pp in 1:p) {
            HBX[idxB, pp] <- HBX[idxB, pp] + X[use, pp] / phi.hat.bin[use]
          }
        }
      }
    }
  }

  expected_bin <- list("HAHA" = as(HAHA, "symmetricMatrix"),
                        "HAHB" = Matrix(HAHB, sparse = TRUE),
                        "HBHB" = as(HBHB, "symmetricMatrix"),
                        "HAX" = HAX %*% xtx_bin,
                        "HBX" = HBX %*% xtx_bin)

  expect_equal(.computeCov.Matrix(y_bin, X, xtx_bin, phi.hat.bin, M, HA, HB),
               expected_bin)
})


test_that(".computeCov.FBM works",{
  library(Matrix)
  library(bigstatsr)

  withr::with_seed(12345,
                   {n <- 10
                   covData <- data.frame(
                     ID = paste0("P", 1:n),
                     X1 = runif(n, 0, 1),
                     X2 = sample(0:1, n, replace = TRUE)
                   )

                   z <- 0.5 * covData$X1 + 0.3 * covData$X2 + rnorm(n)

                   # Simulate phenotype data
                   Y <- data.frame(
                     ID = paste0("P", 1:n),
                     Yc = z,
                     Yb = rbinom(n, size = 1, prob = 1 / (1 + exp(-z)))
                   )})

  formula_c <- Yc ~ X1 + X2
  formula_b <- Yb ~ X1 + X2

  y_cont <- .getYstar(covData, Y, formula_c, "cont")
  y_bin <- .getYstar(covData, Y, formula_b, "bin")

  X <- stats::model.matrix(stats::update.formula(formula_b, NULL ~ .),
                           covData[order(covData$ID),])
  px <- ncol(X)

  phi.cont <- drop(crossprod(y_cont$Residuals)/(nrow(y_cont)-px))
  phi.hat.cont <- rep(phi.cont, nrow(y_cont))
  phi.hat.bin <- 1.0 / {y_bin$Fitted * {1.0 - y_bin$Fitted}}

  xtx_cont <- .computeXTXinvsqrt(X,phi.hat.cont)
  xtx_bin <- .computeXTXinvsqrt(X,phi.hat.bin)

  test_cnv <- data.frame("ID" = c("P1","P2","P3"),
                         "CHR" = c(1,1,1),
                         "BP1" = c(150,200,200),
                         "BP2" = c(200,300,300),
                         "TYPE" = c(4L,1L,1L))

  func.dup <- function(x, cutoff.dup) { as.integer(x > {cutoff.dup - 1e-8}) }
  func.del <- function(x, cutoff.del) { as.integer(x < {cutoff.del + 1e-8}) }

  grid <- .defineGrid(test_cnv)
  M <- .matchGrid(test_cnv, grid)

  HA <- c(do.call(func.dup, list(M$TYPE, 2.3)),
          do.call(func.del, list(M$TYPE, 1.7)))
  HB <- c(do.call(func.dup, list(M$TYPE, 2.3 + 1.0)),
          do.call(func.del, list(M$TYPE, 1.7 - 1.0)))

  Matrix.cont <- .computeCov.Matrix(y_cont, X, xtx_cont, phi.hat.cont, M, HA, HB)
  Matrix.bin <- .computeCov.Matrix(y_bin, X, xtx_bin, phi.hat.bin, M, HA, HB)

  mat.cont <- lapply(Matrix.cont, function(mat) {
    mat_dense <- as.matrix(mat)
    mat_dense
  })

  mat.bin <- lapply(Matrix.bin, function(mat) {
    mat_dense <- as.matrix(mat)
    mat_dense
  })

  #save_dir <- "~/Documents/PackageDev/MetaCONCUR/tests/testthat/tmp_share"
  tmp <- withr::local_tempdir()
  save_dir1 <- file.path(tmp, "study1")

  act.cont <- .computeCov.FBM(y_cont, X, xtx_cont, phi.hat.cont, M, HA, HB,
                              save_dir1)
  act.mat.cont <- lapply(act.cont, function(fbm) {
    mat_dense <- as.matrix(fbm[, ])
    mat_dense
  })

  save_dir2 <- file.path(tmp, "study2")
  act.bin <- .computeCov.FBM(y_bin, X, xtx_bin, phi.hat.bin, M, HA, HB,
                             save_dir2)

  # 1. Return a list
  expect_type(act.cont, "list")
  # 2. Check all elements are FBM objects
  expect_true(all(vapply(act.cont, inherits, logical(1), "FBM")))
  # 1. Return a list
  expect_type(act.bin, "list")
  # 2. Check all elements are FBM objects
  expect_true(all(vapply(act.bin, inherits, logical(1), "FBM")))

  act.mat.bin <- lapply(act.bin, function(fbm) {
    mat_dense <- as.matrix(fbm[, ])
    mat_dense
  })

  expect_equal(act.mat.cont, mat.cont)

  expect_equal(act.mat.bin, mat.bin)

  # Clean up
  #unlink(save_dir1, recursive = TRUE, force = TRUE)
  #unlink(save_dir2, recursive = TRUE, force = TRUE)
})

test_that(".computeCov works",{

  library(Matrix)
  library(bigstatsr)

  withr::with_seed(12345,
                   {n <- 10
                   covData <- data.frame(
                     ID = paste0("P", 1:n),
                     X1 = runif(n, 0, 1),
                     X2 = sample(0:1, n, replace = TRUE)
                   )

                   z <- 0.5 * covData$X1 + 0.3 * covData$X2 + rnorm(n)

                   # Simulate phenotype data
                   Y <- data.frame(
                     ID = paste0("P", 1:n),
                     Yc = z,
                     Yb = rbinom(n, size = 1, prob = 1 / (1 + exp(-z)))
                   )})

  formula_c <- Yc ~ X1 + X2
  formula_b <- Yb ~ X1 + X2

  y_cont <- .getYstar(covData, Y, formula_c, "cont")
  y_bin <- .getYstar(covData, Y, formula_b, "bin")

  X <- stats::model.matrix(stats::update.formula(formula_b, NULL ~ .),
                           covData[order(covData$ID),])
  px <- ncol(X)

  phi.cont <- drop(crossprod(y_cont$Residuals)/(nrow(y_cont)-px))
  phi.hat.cont <- rep(phi.cont, nrow(y_cont))
  phi.hat.bin <- 1.0 / {y_bin$Fitted * {1.0 - y_bin$Fitted}}

  xtx_cont <- .computeXTXinvsqrt(X,phi.hat.cont)
  xtx_bin <- .computeXTXinvsqrt(X,phi.hat.bin)

  test_cnv <- data.frame("ID" = c("P1","P2","P3"),
                         "CHR" = c(1,1,1),
                         "BP1" = c(150,200,200),
                         "BP2" = c(200,300,300),
                         "TYPE" = c(4L,1L,1L))

  func.dup <- function(x, cutoff.dup) { as.integer(x > {cutoff.dup - 1e-8}) }
  func.del <- function(x, cutoff.del) { as.integer(x < {cutoff.del + 1e-8}) }

  grid <- .defineGrid(test_cnv)
  M <- .matchGrid(test_cnv, grid)

  HA <- c(do.call(func.dup, list(M$TYPE, 2.3)),
          do.call(func.del, list(M$TYPE, 1.7)))
  HB <- c(do.call(func.dup, list(M$TYPE, 2.3 + 1.0)),
          do.call(func.del, list(M$TYPE, 1.7 - 1.0)))

  #save_dir <- "~/Documents/PackageDev/MetaCONCUR/tests/testthat/tmp_share"
  tmp <- withr::local_tempdir()
  save_dir1 <- file.path(tmp, "study1")

  expect_equal(.computeCov(y_bin, X, phi.hat.bin, M, HA, HB, FALSE, NULL),
               .computeCov.Matrix(y_bin, X, xtx_bin, phi.hat.bin, M, HA, HB))

  act.cont <- .computeCov(y_cont, X, phi.hat.cont, M, HA, HB, TRUE, save_dir1)
  act.mat.cont <- lapply(act.cont, function(fbm) {
    mat_dense <- as.matrix(fbm[, ])
    mat_dense
  })

  save_dir2 <- file.path(tmp, "study2")

  exp.cont <- .computeCov.FBM(y_cont, X, xtx_cont, phi.hat.cont, M, HA, HB,
                              save_dir2)
  exp.mat.cont <- lapply(exp.cont, function(fbm) {
    mat_dense <- as.matrix(fbm[, ])
    mat_dense
  })

  expect_equal(act.mat.cont,exp.mat.cont)

  expect_true(all(vapply(act.cont, inherits, logical(1), "FBM")))

  #unlink(save_dir1, recursive = TRUE, force = TRUE)
  #unlink(save_dir2, recursive = TRUE, force = TRUE)

  expect_error(.computeCov(NULL, X, phi.hat.bin, M, HA, HB, FALSE, NULL),
               "`y` must be a data.frame")
  expect_error(.computeCov(y_bin, NULL, phi.hat.bin, M, HA, HB, FALSE, NULL),
               "`X` must be a matrix")
  expect_error(.computeCov(y_bin, X, NULL, M, HA, HB, FALSE, NULL),
               "`phi.hat` must be a numeric vector")
  expect_error(.computeCov(y_bin, X, phi.hat.bin, NULL, HA, HB, FALSE, NULL),
               "`grid.ind` must be a data.frame")
  expect_error(.computeCov(y_bin, X, phi.hat.bin, M, NULL, HB, FALSE, NULL),
               "`HA` must be an integer vector")
  expect_error(.computeCov(y_bin, X, phi.hat.bin, M, HA, NULL, FALSE, NULL),
               "`HB` must be an integer vector")
  expect_error(.computeCov(y_bin, X, phi.hat.bin, M, HA, HB, NULL, NULL),
               "`use.big.memory` must be a logical")

})


test_that(".scoreProcedure works",{

  library(Matrix)
  library(bigstatsr)

  withr::with_seed(12345,
                   {n <- 10
                   covData <- data.frame(
                     ID = paste0("P", 1:n),
                     X1 = runif(n, 0, 1),
                     X2 = sample(0:1, n, replace = TRUE)
                   )

                   z <- 0.5 * covData$X1 + 0.3 * covData$X2 + rnorm(n)

                   # Simulate phenotype data
                   Y <- data.frame(
                     ID = paste0("P", 1:n),
                     Yc = z,
                     Yb = rbinom(n, size = 1, prob = 1 / (1 + exp(-z)))
                   )})

  formula_c <- Yc ~ X1 + X2
  formula_b <- Yb ~ X1 + X2

  y_cont <- .getYstar(covData, Y, formula_c, "cont")
  y_bin <- .getYstar(covData, Y, formula_b, "bin")

  X <- stats::model.matrix(stats::update.formula(formula_b, NULL ~ .),
                           covData[order(covData$ID),])
  px <- ncol(X)

  phi.cont <- drop(crossprod(y_cont$Residuals)/(nrow(y_cont)-px))
  phi.hat.cont <- rep(phi.cont, nrow(y_cont))
  phi.hat.bin <- 1.0 / {y_bin$Fitted * {1.0 - y_bin$Fitted}}

  xtx_cont <- .computeXTXinvsqrt(X,phi.hat.cont)
  xtx_bin <- .computeXTXinvsqrt(X,phi.hat.bin)

  test_cnv <- data.frame("ID" = c("P1","P2","P3"),
                         "CHR" = c(1,1,1),
                         "BP1" = c(150,200,200),
                         "BP2" = c(200,300,300),
                         "TYPE" = c(4L,1L,1L))

  #save_dir <- "~/Documents/PackageDev/MetaCONCUR/tests/testthat/tmp_share"

  func.dup <- function(x, cutoff.dup) { as.integer(x > {cutoff.dup - 1e-8}) }
  func.del <- function(x, cutoff.del) { as.integer(x < {cutoff.del + 1e-8}) }

  grid <- .defineGrid(test_cnv)
  M <- .matchGrid(test_cnv, grid)

  HA <- c(do.call(func.dup, list(M$TYPE, 2.3)),
          do.call(func.del, list(M$TYPE, 1.7)))
  HB <- c(do.call(func.dup, list(M$TYPE, 2.3 + 1.0)),
          do.call(func.del, list(M$TYPE, 1.7 - 1.0)))

  res.cont <- list(grid = grid,
                       SSA = c(.computeSS(y_cont, M, func.dup, 2.3, px, "cont"),
                               .computeSS(y_cont, M, func.del, 1.7, px, "cont")),
                       SSB = c(.computeSS(y_cont, M, func.dup, 2.3 + 1.0, px, "cont"),
                               .computeSS(y_cont, M, func.del, 1.7 - 1.0, px, "cont")))

  exp.cont <- c(res.cont,
                .computeCov(y_cont, X, phi.hat.cont, M, HA, HB, FALSE, NULL))

  expect_equal(.scoreProcedure(test_cnv, 2.3, 1.7,y_cont, X, "cont", FALSE, NULL),
               exp.cont)

  res.bin <- list(grid = grid,
                  SSA = c(.computeSS(y_bin, M, func.dup, 2.3, px, "bin"),
                          .computeSS(y_bin, M, func.del, 1.7, px, "bin")),
                  SSB = c(.computeSS(y_bin, M, func.dup, 2.3 + 1.0, px, "bin"),
                          .computeSS(y_bin, M, func.del, 1.7 - 1.0, px, "bin")))

  tmp <- withr::local_tempdir()
  save_dir1 <- file.path(tmp, "study1")
  save_dir2 <- file.path(tmp, "study2")

  exp.bin <- .computeCov(y_bin, X, phi.hat.bin, M, HA, HB, TRUE, save_dir1)

  exp.mat.bin <- lapply(exp.bin, function(fbm) {
    mat_dense <- as.matrix(fbm[, ])
    mat_dense
  })

  exp.res.bin <- c(res.bin, exp.mat.bin)

  act.bin <- .scoreProcedure(test_cnv, 2.3, 1.7, y_bin, X, "bin", TRUE, save_dir2)

  act.mat.bin <-c(
    act.bin[1:3],
    lapply(act.bin[-(1:3)], function(fbm) {
      mat_dense <- as.matrix(fbm[, ])
      mat_dense
    })
  )

  expect_equal(act.mat.bin, exp.res.bin)

  #unlink(save_dir1, recursive = TRUE, force = TRUE)
  #unlink(save_dir2, recursive = TRUE, force = TRUE)

  expect_error(.scoreProcedure(NULL, 2.3, 1.7, y_bin, X, "bin", FALSE, NULL))
  expect_error(.scoreProcedure(test_cnv, 2.3, 1.7, NULL, X, "bin", FALSE, NULL),
               "`y` must be a data.frame")
  expect_error(.scoreProcedure(test_cnv, 2.3, 1.7, y_bin, NULL, "bin", FALSE, NULL),
               "`X` must be a matrix")
  expect_error(.scoreProcedure(test_cnv, NULL, 1.7, y_bin, X, "bin", FALSE, NULL),
              "`cutoff.dup` must be a scalar")
  expect_error(.scoreProcedure(test_cnv, 2.3, NULL, y_bin, X, "bin", FALSE, NULL),
              "`cutoff.del` must be a scalar")
  expect_error(.scoreProcedure(test_cnv, 2.3, 1.7, y_bin, X, "bin", NULL, NULL),
               "`use.big.memory` must be a logical")
})

test_that(".computeTS works",{
  withr::with_seed(12345,{
    grid <- integer(11)
    grid[1] <- sample(1:100, 1)

    # each next element is at least 25 greater
    for (i in 2:11) {
      grid[i] <- grid[i - 1] + sample(25:50, 1)
    }
    SSA <- runif(20)})

  exp.TS <- sum(SSA * diff(grid) * SSA)

  expect_equal(.computeTS(SSA, grid), exp.TS)
  expect_error(.computeTS(NULL, grid),
               "`score` must be a numeric vector with double length of `grid`")
  expect_error(.computeTS(SSA, NULL))
})

test_that(".save_chr_stats works in-memory mode", {
  tmp <- withr::local_tempdir()  # auto-cleanup after test
  chr_name <- "chr1"

  chr1_s <- list(grid = c(150,200,300),
                 SSA = c(1,0,0,2),
                 SSB = c(1,0,0,0),
                 HAHA = rbind(cbind(Matrix(c(1,0,0,0),nrow = 2,ncol = 2),
                                    Matrix(rep(0,4),nrow = 2,ncol = 2)),
                              cbind(Matrix(rep(0,4),nrow = 2,ncol = 2),
                                    Matrix(c(0,0,0,2),nrow = 2,ncol = 2))),
                 HAHB = forceSymmetric(rbind(cbind(Matrix(c(1,0,0,0),nrow = 2,ncol = 2),
                                                   Matrix(rep(0,4),nrow = 2,ncol = 2)),
                                             cbind(Matrix(rep(0,4),nrow = 2,ncol = 2),
                                                   Matrix(rep(0,4),nrow = 2,ncol = 2)))),
                 HBHB = forceSymmetric(rbind(cbind(Matrix(c(1,0,0,0),nrow = 2,ncol = 2),
                                                   Matrix(rep(0,4),nrow = 2,ncol = 2)),
                                             cbind(Matrix(rep(0,4),nrow = 2,ncol = 2),
                                                   Matrix(rep(0,4),nrow = 2,ncol = 2)))),
                 HAX = Matrix(c(1,0,0,1,0,0,0,1),nrow = 4,ncol = 2),
                 HBX = Matrix(c(1,0,0,0,0,0,0,0),nrow = 4,ncol = 2)
  )

  chr_dir <- file.path(tmp, chr_name)
  dir.create(chr_dir)

  .save_chr_stats(chr1_s, chr_name, use.big.memory = FALSE, save_dir = tmp)

  # check file exists
  expect_true(file.exists(file.path(chr_dir, "SS.rds")))

  reloaded <- readRDS(file.path(chr_dir, "SS.rds"))
  expect_equal(reloaded$SSA, chr1_s$SSA)
  expect_equal(reloaded$SSB, chr1_s$SSB)
})

test_that(".save_chr_stats works in big memory mode", {
  tmp <- withr::local_tempdir()
  chr_name <- "chr1"
  chr_dir <- file.path(tmp, chr_name)
  dir.create(chr_dir)

  # create fake FBMs in chr_dir
  HAHA <- bigstatsr::FBM(4, 4, backingfile = file.path(chr_dir, "HAHA"))
  HAHB <- bigstatsr::FBM(4, 4, backingfile = file.path(chr_dir, "HAHB"))
  HBHB <- bigstatsr::FBM(4, 4, backingfile = file.path(chr_dir, "HBHB"))
  HAX  <- bigstatsr::FBM(4, 4, backingfile = file.path(chr_dir, "HAX"))
  HBX  <- bigstatsr::FBM(4, 4, backingfile = file.path(chr_dir, "HBX"))

  stats <- list(
    grid = c(150,200,300),
    SSA = c(1,0,0,2),
    SSB = c(1,0,0,0),
    HAHA = HAHA, HAHB = HAHB, HBHB = HBHB, HAX = HAX, HBX = HBX
  )

  .save_chr_stats(stats, chr_name, use.big.memory = TRUE, save_dir = tmp)

  expect_true(file.exists(file.path(chr_dir, "wo_FBM.rds")))

  reloaded_small <- readRDS(file.path(chr_dir, "wo_FBM.rds"))
  expect_equal(reloaded_small$SSA, stats$SSA)
  expect_equal(reloaded_small$SSB, stats$SSB)
})

test_that("CONCUR_worker works",{
  library(Matrix)
  library(bigstatsr)

  withr::with_seed(12345,
                   {n <- 10
                   covData <- data.frame(
                     ID = paste0("P", 1:n),
                     X1 = runif(n, 0, 1),
                     X2 = sample(0:1, n, replace = TRUE)
                   )

                   z <- 0.5 * covData$X1 + 0.3 * covData$X2 + rnorm(n)

                   # Simulate phenotype data
                   Y <- data.frame(
                     ID = paste0("P", 1:n),
                     Yc = z,
                     Yb = rbinom(n, size = 1, prob = 1 / (1 + exp(-z)))
                   )})

  formula_c <- Yc ~ X1 + X2
  formula_b <- Yb ~ X1 + X2

  y_cont <- .getYstar(covData, Y, formula_c, "cont")
  y_bin <- .getYstar(covData, Y, formula_b, "bin")

  X <- stats::model.matrix(stats::update.formula(formula_b, NULL ~ .),
                           covData[order(covData$ID),])

  test_cnv <- data.frame("ID" = c("P1","P1","P2","P3","P4"),
                         "CHR" = c(1L,2L,1L,1L,2L),
                         "BP1" = c(150,400,200,200,450),
                         "BP2" = c(200,500,300,300,500),
                         "TYPE" = c(4L,3L,1L,1L,3L))

  #save_dir <- "~/Documents/PackageDev/MetaCONCUR/tests/testthat/tmp_share"

  res.cont <- list(chr1 = .scoreProcedure(test_cnv[test_cnv$CHR == 1L, ],
                                          2.3, 1.7, y_cont, X, "cont", FALSE,
                                          NULL),
                   chr2 = .scoreProcedure(test_cnv[test_cnv$CHR == 2L, ],
                                          2.3, 1.7, y_cont, X, "cont", FALSE,
                                          NULL))

  res.bin <- list(chr1 = .scoreProcedure(test_cnv[test_cnv$CHR == 1L, ],
                                         2.3, 1.7, y_bin, X, "bin", FALSE,
                                         NULL))

  expect_equal(CONCUR_worker(test_cnv, Y, 2.3, 1.7, formula_c, covData, "cont",
                             "all", FALSE, NULL), res.cont)

  expect_equal(CONCUR_worker(test_cnv, Y, 2.3, 1.7, formula_b, covData, "bin",
                             1L, FALSE, NULL), res.bin)

  expect_error(CONCUR_worker(NULL, Y, 2.3, 1.7, formula_c, covData, "cont",
                             "all", FALSE, NULL))

  expect_error(CONCUR_worker(test_cnv, NULL, 2.3, 1.7, formula_c, covData, "cont",
                             "all", FALSE, NULL), "`phenoData` must be a data.frame")
  expect_error(CONCUR_worker(test_cnv, Y, NULL, 1.7, formula_c, covData, "cont",
                             "all", FALSE, NULL), "`cutoff.dup` must be a scalar")
  expect_error(CONCUR_worker(test_cnv, Y, 2.3, NULL, formula_c, covData, "cont",
                             "all", FALSE, NULL), "`cutoff.del` must be a scalar")
  expect_error(CONCUR_worker(test_cnv, Y, 2.3, 1.7, formula_c, NULL, "cont",
                             "all", FALSE, NULL), "`covData` must be a data.frame")
  expect_error(CONCUR_worker(test_cnv, Y, 2.3, 1.7, NULL, covData, "cont",
                             "all", FALSE, NULL), "`formula` must be a formula of the form LHS ~ RHS")
  expect_error(CONCUR_worker(test_cnv, Y, 2.3, 1.7, formula_c, covData, "cont",
                             "all", NULL, NULL), "`use.big.memory` must be a logical")
  expect_error(CONCUR_worker(test_cnv, Y, 2.3, 1.7, formula_c, covData, "cont",
                             "all", TRUE, NULL), "`save_dir` must be specified")
})
