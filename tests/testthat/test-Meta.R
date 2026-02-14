test_that(".computeTSMeta works", {
  library(Matrix)

  study1 <- list(
    DUP = list(
      grid = c(150,200,300),
      SSA = c(1,0),
      SSB = c(1,0),
      HAHA = Matrix(c(1,0,0,0),nrow = 2,ncol = 2),
      HAHB = Matrix(c(1,0,0,0),nrow = 2,ncol = 2),
      HBHB = Matrix(c(1,0,0,0),nrow = 2,ncol = 2),
      HAX = Matrix(c(1,0,0,1),nrow = 2,ncol = 2),
      HBX = Matrix(c(1,0,0,0),nrow = 2,ncol = 2)),
    DEL = list(
      grid = c(150,200,300),
      SSA = c(0,2),
      SSB = c(0,0),
      HAHA = Matrix(c(0,0,0,2),nrow = 2,ncol = 2),
      HAHB = Matrix(rep(0,4),nrow = 2,ncol = 2),
      HBHB = Matrix(rep(0,4),nrow = 2,ncol = 2),
      HAX = Matrix(c(0,0,0,1),nrow = 2,ncol = 2),
      HBX = Matrix(c(0,0,0,0),nrow = 2,ncol = 2))
  )

  study2 <- list(
    DUP = list(grid = c(100,250,350),
                 SSA = c(2,1),
                 SSB = c(0,0),
                 HAHA = Matrix(c(0,1,0,1),nrow = 2,ncol = 2),
                 HAHB = Matrix(c(0,0,0,1),nrow = 2,ncol = 2),
                 HBHB = Matrix(c(0,0,0,1),nrow = 2,ncol = 2),
                 HAX = Matrix(c(0,1,0,1),nrow = 2,ncol = 2),
                 HBX = Matrix(c(0,0,0,1),nrow = 2,ncol = 2)),
    DEL = list(grid = c(100,250,350),
               SSA = c(0,1),
               SSB = c(0,1),
               HAHA = Matrix(c(2,0,0,0),nrow = 2,ncol = 2),
               HAHB = Matrix(rep(0,4),nrow = 2,ncol = 2),
               HBHB = Matrix(rep(0,4),nrow = 2,ncol = 2),
               HAX = Matrix(c(2,0,0,0),nrow = 2,ncol = 2),
               HBX = Matrix(c(0,0,0,0),nrow = 2,ncol = 2))
  )

  grid.meta <- c(100,150,200,250,300,350)

  study_SS <- list("DUP" = list(), "DEL" = list())

  study_SS$DUP <- rlist::list.append(study_SS$DUP, study1$DUP, study2$DUP)
  study_SS$DEL <- rlist::list.append(study_SS$DEL, study1$DEL, study2$DEL)
  
  meta.data1_DUP <- lapply(study_SS$DUP, .alignStudyToMeta, meta.grid = grid.meta)
  meta.data1_DEL <- lapply(study_SS$DEL, .alignStudyToMeta, meta.grid = grid.meta)
  meta.data1 <- list(meta = list(), "DUP" = meta.data1_DUP, "DEL" = meta.data1_DEL)

  meta.data2 <- meta.data1

  sum_SSA_DUP <- meta.data2$DUP[[1L]]$SSA + meta.data2$DUP[[2L]]$SSA
  sum_SSA_DEL <- meta.data2$DEL[[1L]]$SSA + meta.data2$DEL[[2L]]$SSA
  sum_SSB_DUP <- meta.data2$DUP[[1L]]$SSB + meta.data2$DUP[[2L]]$SSB
  sum_SSB_DEL <- meta.data2$DEL[[1L]]$SSB + meta.data2$DEL[[2L]]$SSB
  
  meta.data2$DUP[[1L]]$SSA <- NULL
  meta.data2$DUP[[1L]]$SSB <- NULL
  meta.data2$DUP[[2L]]$SSA <- NULL
  meta.data2$DUP[[2L]]$SSB <- NULL
  meta.data2$DEL[[1L]]$SSA <- NULL
  meta.data2$DEL[[1L]]$SSB <- NULL
  meta.data2$DEL[[2L]]$SSA <- NULL
  meta.data2$DEL[[2L]]$SSB <- NULL

  meta.data2$DUP$meta$TS <- .computeTS(sum_SSA_DUP, grid.meta) + .computeTS(sum_SSB_DUP, grid.meta)
  meta.data2$DEL$meta$TS <- .computeTS(sum_SSA_DEL, grid.meta)+ .computeTS(sum_SSB_DEL, grid.meta)
    
  expect_equal(.computeTSMeta(meta.data1$DUP, grid.meta), meta.data2$DUP)
  
  meta.data2$meta$TS <- .computeTS(sum_SSA_DEL, grid.meta) +
                            .computeTS(sum_SSB_DEL, grid.meta)
  
  
  expect_equal(.computeTSMeta(meta.data1$DEL, grid.meta), meta.data2$DEL)
  
})

test_that(".eigen.Matrix works",{
  library(Matrix)
  library(rlang)
  library(RSpectra)

  withr::with_seed(12345,{
    # create a random sparse matrix 25 by 25
    Phi <- rsparsematrix(nrow = 25, ncol = 25, density = 0.1)

    # Make it symmetric
    Phi_sym <- forceSymmetric(Phi)

    # Ensure it's in a sparse symmetric class
    Phi_sym <- as(Phi_sym, "generalMatrix")
  })

#  test1 <- RSpectra::eigs_sym(Phi_sym, 24, opts = list(retvec = FALSE))$values

  expect.values <- eigen(Phi_sym, only.values = TRUE)$values

  expect_equal(.eigen.Matrix(Phi_sym, .force_method = "base"), expect.values)


  expect_warning(vals <- .eigen.Matrix(Phi_sym, .force_method = "rspectra"),
                 "Phi is expected to be sparse, the input Phi is not sparse")

  expect_true(is.numeric(vals))

  expect_equal(vals, expect.values)
})

test_that(".eigen.FBM works", {
  library(Matrix)
  library(bigstatsr)

  withr::with_seed(12345,{
    # create a random sparse matrix 25 by 25
    Phi <- rsparsematrix(nrow = 25, ncol = 25, density = 0.1)

    # Make it symmetric
    Phi_sym <- forceSymmetric(Phi)
  })

  Phi1 <- as.matrix(Phi_sym)
  Phi2 <- as_FBM(Phi1)

  vals <- sort(.eigen.Matrix(Phi_sym), decreasing = TRUE)
  expect_equal(sort(.eigen.FBM(Phi2), decreasing = TRUE), vals, tolerance = 1e-6)

})

test_that(".getPValues works", {

  withr::with_seed(12345,{
    # 10 eigenvalues between -2 and 2
    eigvals <- runif(10, min = -2, max = 2)
    q <- runif(1, min = 10, max = 20)
  })
  exp.val <- c("davies" = CompQuadForm::davies(q, lambda = eigvals)$Qq,
               "liu" = CompQuadForm::liu(q, lambda = eigvals))

  expect_equal(.getPValues(q, eigvals), exp.val)

  expect_error(.getPValues(NULL, eigvals), "`q` must be a numeric vector")

  expect_error(.getPValues(q, NULL), "`lambda` must be a numeric vector")
})

test_that("Meta works",{
  library(Matrix)
  library(bigstatsr)

  withr::with_seed(12345,{
    # create a random sparse matrix 20 by 20
    HA1 <- rsparsematrix(nrow = 10, ncol = 10, density = 0.1)
    HB1 <- rsparsematrix(nrow = 10, ncol = 10, density = 0.1)
    HAB1 <- rsparsematrix(nrow = 10, ncol = 10, density = 0.1)
    grid <- cumsum(c(sample(101:200, 1),
                     sample(50:100, size = 10, replace = TRUE)))

    study1 <- list(
      DUP = list(grid = grid,
                 SSA = sample(1:10, 10, replace = TRUE),
                 SSB = sample(1:10, 10, replace = TRUE),
                 HAHA = (HA1 + t(HA1)) / 2,
                 HAHB = (HAB1 + t(HAB1)) / 2,
                 HBHB = (HB1 +t(HB1)) / 2,
                 HAX = sparseMatrix(i = sample(1:10, 5, replace = TRUE),
                                    j = sample(1:5, 5, replace = TRUE),
                                    x = rnorm(5), dims = c(10, 5)),
                 HBX = sparseMatrix(i = sample(1:10, 5, replace = TRUE),
                                    j = sample(1:5, 5, replace = TRUE),
                                    x = rnorm(5), dims = c(10, 5))),
      DEL = list(grid = grid,
                 SSA = sample(1:10, 10, replace = TRUE),
                 SSB = sample(1:10, 10, replace = TRUE),
                 HAHA = (HA1 + t(HA1)) / 2,
                 HAHB = (HAB1 + t(HAB1)) / 2,
                 HBHB = (HB1 +t(HB1)) / 2,
                 HAX = sparseMatrix(i = sample(1:10, 5, replace = TRUE),
                                    j = sample(1:5, 5, replace = TRUE),
                                    x = rnorm(5), dims = c(10, 5)),
                 HBX = sparseMatrix(i = sample(1:10, 5, replace = TRUE),
                                    j = sample(1:5, 5, replace = TRUE),
                                    x = rnorm(5), dims = c(10, 5)))
    )

    HA2 <- rsparsematrix(nrow = 10, ncol = 10, density = 0.1)
    HB2 <- rsparsematrix(nrow = 10, ncol = 10, density = 0.1)
    HAB2 <- rsparsematrix(nrow = 10, ncol = 10, density = 0.1)
    grid <- cumsum(c(sample(101:200, 1),
                     sample(50:100, size = 10, replace = TRUE)))
    study2 <- list(
      DUP = list(grid = grid,
                   SSA = sample(1:10, 10, replace = TRUE),
                   SSB = sample(1:10, 10, replace = TRUE),
                   HAHA = (HA2 + t(HA2)) / 2,
                   HAHB = (HAB2 + t(HAB2)) / 2,
                   HBHB = (HB2 +t(HB2)) / 2,
                   HAX = sparseMatrix(i = sample(1:10, 5, replace = TRUE),
                                      j = sample(1:5, 5, replace = TRUE),
                                      x = rnorm(5), dims = c(10, 5)),
                   HBX = sparseMatrix(i = sample(1:10, 5, replace = TRUE),
                                      j = sample(1:5, 5, replace = TRUE),
                                      x = rnorm(5), dims = c(10, 5))
    ),
    DEL = list(grid = grid,
                 SSA = sample(1:10, 10, replace = TRUE),
                 SSB = sample(1:10, 10, replace = TRUE),
                 HAHA = (HA2 + t(HA2)) / 2,
                 HAHB = (HAB2 + t(HAB2)) / 2,
                 HBHB = (HB2 +t(HB2)) / 2,
                 HAX = sparseMatrix(i = sample(1:10, 5, replace = TRUE),
                                    j = sample(1:5, 5, replace = TRUE),
                                    x = rnorm(5), dims = c(10, 5)),
                 HBX = sparseMatrix(i = sample(1:10, 5, replace = TRUE),
                                    j = sample(1:5, 5, replace = TRUE),
                                    x = rnorm(5), dims = c(10, 5))
    ))
  })

  meta_grid <- list(study1$DEL$grid,study2$DEL$grid) |> unlist() |> unique() |> sort()

  study_SS <- list("DUP" = list(), "DEL" = list())
  
  study_SS$DUP <- rlist::list.append(study_SS$DUP, study1$DUP, study2$DUP)
  study_SS$DEL <- rlist::list.append(study_SS$DEL, study1$DEL, study2$DEL)
  
  meta.data1_DUP <- lapply(study_SS$DUP, .alignStudyToMeta, meta.grid = meta_grid)
  meta.data1_DEL <- lapply(study_SS$DEL, .alignStudyToMeta, meta.grid = meta_grid)
  meta.data1 <- list(meta = list(), "DUP" = meta.data1_DUP, "DEL" = meta.data1_DEL)
  
  meta.data1$DUP <- .computeTSMeta(meta.data1$DUP, meta_grid)
  meta.data1$DUP <- .sumHAHA(meta.data1$DUP)
  meta.data1$DUP <- .sumHAHB(meta.data1$DUP)
  meta.data1$DUP <- .sumHBHB(meta.data1$DUP)
  meta.data1$DUP <- .rmStudies(meta.data1$DUP)

  meta.data1$DEL <- .computeTSMeta(meta.data1$DEL, meta_grid)
  meta.data1$DEL <- .sumHAHA(meta.data1$DEL)
  meta.data1$DEL <- .sumHAHB(meta.data1$DEL)
  meta.data1$DEL <- .sumHBHB(meta.data1$DEL)
  meta.data1$DEL <- .rmStudies(meta.data1$DEL)
  
  meta.data1$TSmeta <- meta.data1$DEL$meta$TS + meta.data1$DUP$meta$TS

  meta.data1$meta <- .combineDELDUP(meta.data1$DEL, meta.data1$DUP)
  
  Phi <- .combineBlocks(meta.data1$meta$meta, meta_grid)

  ee_meta <- .eigen.Matrix(Phi)

  exp.res <- .getPValues(q = meta.data1$meta$meta$TS, lambda = ee_meta, method = "meta")

  expect_equal(Meta(study1, study2), exp.res)
  withr::defer({
    rm(list = ls())
    gc()
  }, envir = environment())

})

test_that("Meta works",{
  library(Matrix)
  library(bigstatsr)
  
  withr::with_seed(12345,{
    # create a random sparse matrix 20 by 20
    HA1 <- rsparsematrix(nrow = 10, ncol = 10, density = 0.1)
    HB1 <- rsparsematrix(nrow = 10, ncol = 10, density = 0.1)
    HAB1 <- rsparsematrix(nrow = 10, ncol = 10, density = 0.1)
    grid <- cumsum(c(sample(101:200, 1),
                     sample(50:100, size = 10, replace = TRUE)))
    
    study1 <- list(
      DUP = list(grid = grid,
                 SSA = sample(1:10, 10, replace = TRUE),
                 SSB = sample(1:10, 10, replace = TRUE),
                 HAHA = (HA1 + t(HA1)) / 2,
                 HAHB = (HAB1 + t(HAB1)) / 2,
                 HBHB = (HB1 +t(HB1)) / 2,
                 HAX = sparseMatrix(i = sample(1:10, 5, replace = TRUE),
                                    j = sample(1:5, 5, replace = TRUE),
                                    x = rnorm(5), dims = c(10, 5)),
                 HBX = sparseMatrix(i = sample(1:10, 5, replace = TRUE),
                                    j = sample(1:5, 5, replace = TRUE),
                                    x = rnorm(5), dims = c(10, 5))),
      DEL = list(grid = grid,
                 SSA = sample(1:10, 10, replace = TRUE),
                 SSB = sample(1:10, 10, replace = TRUE),
                 HAHA = (HA1 + t(HA1)) / 2,
                 HAHB = (HAB1 + t(HAB1)) / 2,
                 HBHB = (HB1 +t(HB1)) / 2,
                 HAX = sparseMatrix(i = sample(1:10, 5, replace = TRUE),
                                    j = sample(1:5, 5, replace = TRUE),
                                    x = rnorm(5), dims = c(10, 5)),
                 HBX = sparseMatrix(i = sample(1:10, 5, replace = TRUE),
                                    j = sample(1:5, 5, replace = TRUE),
                                    x = rnorm(5), dims = c(10, 5)))
    )
    
    HA2 <- rsparsematrix(nrow = 10, ncol = 10, density = 0.1)
    HB2 <- rsparsematrix(nrow = 10, ncol = 10, density = 0.1)
    HAB2 <- rsparsematrix(nrow = 10, ncol = 10, density = 0.1)
    grid <- cumsum(c(sample(101:200, 1),
                     sample(50:100, size = 10, replace = TRUE)))
    study2 <- list(
      DUP = list(grid = grid,
                 SSA = sample(1:10, 10, replace = TRUE),
                 SSB = sample(1:10, 10, replace = TRUE),
                 HAHA = (HA2 + t(HA2)) / 2,
                 HAHB = (HAB2 + t(HAB2)) / 2,
                 HBHB = (HB2 +t(HB2)) / 2,
                 HAX = sparseMatrix(i = sample(1:10, 5, replace = TRUE),
                                    j = sample(1:5, 5, replace = TRUE),
                                    x = rnorm(5), dims = c(10, 5)),
                 HBX = sparseMatrix(i = sample(1:10, 5, replace = TRUE),
                                    j = sample(1:5, 5, replace = TRUE),
                                    x = rnorm(5), dims = c(10, 5))
      ),
      DEL = list(grid = grid,
                 SSA = sample(1:10, 10, replace = TRUE),
                 SSB = sample(1:10, 10, replace = TRUE),
                 HAHA = (HA2 + t(HA2)) / 2,
                 HAHB = (HAB2 + t(HAB2)) / 2,
                 HBHB = (HB2 +t(HB2)) / 2,
                 HAX = sparseMatrix(i = sample(1:10, 5, replace = TRUE),
                                    j = sample(1:5, 5, replace = TRUE),
                                    x = rnorm(5), dims = c(10, 5)),
                 HBX = sparseMatrix(i = sample(1:10, 5, replace = TRUE),
                                    j = sample(1:5, 5, replace = TRUE),
                                    x = rnorm(5), dims = c(10, 5))
      ))
  })
  
  meta_grid <- list(study1$DEL$grid,study2$DEL$grid) |> unlist() |> unique() |> sort()
  
  study_SS <- list("DUP" = list(), "DEL" = list())
  
  study_SS$DUP <- rlist::list.append(study_SS$DUP, study1$DUP, study2$DUP)
  study_SS$DEL <- rlist::list.append(study_SS$DEL, study1$DEL, study2$DEL)
  
  meta.data1_DUP <- lapply(study_SS$DUP, .alignStudyToMeta, meta.grid = meta_grid)
  meta.data1_DEL <- lapply(study_SS$DEL, .alignStudyToMeta, meta.grid = meta_grid)
  meta.data1 <- list(meta = list(), "DUP" = meta.data1_DUP, "DEL" = meta.data1_DEL)
  
  # test for FBM
  study_SS2 <- study_SS
  for (si in 1:2) {
    for (sj in c("HAHA", "HAHB", "HBHB", "HAX", "HBX")) {
      study_SS2$DEL[[si]][[sj]] <- as.matrix(study_SS2$DEL[[si]][[sj]])
      study_SS2$DUP[[si]][[sj]] <- as.matrix(study_SS2$DUP[[si]][[sj]])
      study_SS2$DEL[[si]][[sj]] <- as_FBM(study_SS2$DEL[[si]][[sj]])
      study_SS2$DUP[[si]][[sj]] <- as_FBM(study_SS2$DUP[[si]][[sj]])
      
      meta.data1$DEL[[si]][[sj]] <- as.matrix(meta.data1$DEL[[si]][[sj]])
      meta.data1$DUP[[si]][[sj]] <- as.matrix(meta.data1$DUP[[si]][[sj]])
      meta.data1$DEL[[si]][[sj]] <- as_FBM(meta.data1$DEL[[si]][[sj]])
      meta.data1$DUP[[si]][[sj]] <- as_FBM(meta.data1$DUP[[si]][[sj]])
    }
  }
  
  meta.data1$DUP <- .computeTSMeta(meta.data1$DUP, meta_grid)
  meta.data1$DUP <- .sumHAHA(meta.data1$DUP)
  meta.data1$DUP <- .sumHAHB(meta.data1$DUP)
  meta.data1$DUP <- .sumHBHB(meta.data1$DUP)
  meta.data1$DUP <- .rmStudies(meta.data1$DUP)
  
  meta.data1$DEL <- .computeTSMeta(meta.data1$DEL, meta_grid)
  meta.data1$DEL <- .sumHAHA(meta.data1$DEL)
  meta.data1$DEL <- .sumHAHB(meta.data1$DEL)
  meta.data1$DEL <- .sumHBHB(meta.data1$DEL)
  meta.data1$DEL <- .rmStudies(meta.data1$DEL)
  
  meta.data1$TSmeta <- meta.data1$DEL$meta$TS + meta.data1$DUP$meta$TS
  meta.data1$meta <- .combineDELDUP(meta.data1$DEL, meta.data1$DUP)
  
  Phi2 <- .combineBlocks(meta.data1$meta$meta, meta_grid)
  
  ee_meta2 <- .eigen.FBM(Phi2)
  
  exp.res2 <- .getPValues(q = meta.data1$meta$meta$TS, lambda = ee_meta2, method = "meta")
  
  expect_equal(Meta(study1,study2), exp.res2)
  
  expect_error(Meta(NULL))
})

test_that("load_chr_stats works with in-memory mode", {
  tmp <- withr::local_tempdir()
  chr_name <- "chr1"
  chr_dir <- file.path(tmp, chr_name)
  dir.create(chr_dir)

  # Fake stats
  stats <- list(
    grid = c(150,200,300),
    SSA = c(1,0,0,2),
    SSB = c(1,0,0,0),
    HAHA = matrix(1, 4, 4),
    HAHB = matrix(2, 4, 4),
    HBHB = matrix(3, 4, 4),
    HAX  = matrix(4, 4, 4),
    HBX  = matrix(5, 4, 4)
  )

  # Save using your save helper
  .save_chr_stats(stats, chr_name, use.big.memory = FALSE, save_dir = tmp)

  # Reload
  reloaded <- load_chr_stats(chr_name, tmp, use.big.memory = FALSE)

  # Check structure and values
  expect_equal(names(reloaded), names(stats))
  expect_equal(reloaded$grid, stats$grid)
  expect_equal(reloaded$SSA, stats$SSA)
  expect_equal(reloaded$SSB, stats$SSB)
  expect_equal(reloaded$HAHA, stats$HAHA)
})


test_that("load_chr_stats works with big-memory mode", {
  tmp <- withr::local_tempdir()
  chr_name <- "chr1"
  chr_dir <- file.path(tmp, chr_name)
  dir.create(chr_dir)

  # Create fake FBMs
  HAHA <- bigstatsr::FBM(4, 4, backingfile = file.path(chr_dir, "HAHA"))
  HAHA$save()
  HAHB <- bigstatsr::FBM(4, 4, backingfile = file.path(chr_dir, "HAHB"))
  HAHB$save()
  HBHB <- bigstatsr::FBM(4, 4, backingfile = file.path(chr_dir, "HBHB"))
  HBHB$save()
  HAX  <- bigstatsr::FBM(4, 4, backingfile = file.path(chr_dir, "HAX"))
  HAX$save()
  HBX  <- bigstatsr::FBM(4, 4, backingfile = file.path(chr_dir, "HBX"))
  HBX$save()

  stats <- list(
    grid = c(150,200,300),
    SSA = c(1,0,0,2),
    SSB = c(1,0,0,0),
    HAHA = HAHA, HAHB = HAHB, HBHB = HBHB, HAX = HAX, HBX = HBX
  )

  # Save using helper
  .save_chr_stats(stats, chr_name, use.big.memory = TRUE, save_dir = tmp)

  # Reload
  reloaded <- load_chr_stats(chr_name, tmp, use.big.memory = TRUE)

  # Check structure
  expect_equal(names(reloaded), names(stats))

  # Check small stats reload properly
  expect_equal(reloaded$SSA, stats$SSA)
  expect_equal(reloaded$SSB, stats$SSB)

  # Check FBMs reload as FBM class
  expect_s4_class(reloaded$HAHA, "FBM")
  expect_s4_class(reloaded$HBX, "FBM")
})
