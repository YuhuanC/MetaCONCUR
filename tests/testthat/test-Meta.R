test_that(".computeTSMeta works", {
  library(Matrix)

  study1 <- list(grid = c(150,200,300),
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

  study2 <- list(grid = c(100,250,350),
                 SSA = c(2,1,0,1),
                 SSB = c(0,0,0,1),
                 HAHA = rbind(cbind(Matrix(c(0,1,0,1),nrow = 2,ncol = 2),
                                    Matrix(rep(0,4),nrow = 2,ncol = 2)),
                              cbind(Matrix(rep(0,4),nrow = 2,ncol = 2),
                                    Matrix(c(2,0,0,0),nrow = 2,ncol = 2))),
                 HAHB = forceSymmetric(rbind(cbind(Matrix(c(0,0,0,1),nrow = 2,ncol = 2),
                                                   Matrix(rep(0,4),nrow = 2,ncol = 2)),
                                             cbind(Matrix(rep(0,4),nrow = 2,ncol = 2),
                                                   Matrix(rep(0,4),nrow = 2,ncol = 2)))),
                 HBHB = forceSymmetric(rbind(cbind(Matrix(c(0,0,0,1),nrow = 2,ncol = 2),
                                                   Matrix(rep(0,4),nrow = 2,ncol = 2)),
                                             cbind(Matrix(rep(0,4),nrow = 2,ncol = 2),
                                                   Matrix(rep(0,4),nrow = 2,ncol = 2)))),
                 HAX = Matrix(c(0,1,0,1,2,0,0,0),nrow = 4,ncol = 2),
                 HBX = Matrix(c(0,0,0,1,0,0,0,0),nrow = 4,ncol = 2)
  )

  grid.meta <- c(100,150,200,250,300,350)

  study_SS <- list()

  study_SS <- rlist::list.append(study_SS, study1)

  study_SS <- rlist::list.append(study_SS, study2)

  meta.data1 <- lapply(study_SS, .alignStudyToMeta, meta.grid = grid.meta)
  meta.data1$meta <- list()

  meta.data2 <- meta.data1

  sum_SSA <- meta.data2[[1L]]$SSA + meta.data2[[2L]]$SSA
  sum_SSB <- meta.data2[[1L]]$SSB + meta.data2[[2L]]$SSB

  meta.data2[[1L]]$SSA <- NULL
  meta.data2[[1L]]$SSB <- NULL
  meta.data2[[2L]]$SSA <- NULL
  meta.data2[[2L]]$SSB <- NULL

  meta.data2$meta$TS <- .computeTS(sum_SSA, grid.meta) +
    .computeTS(sum_SSB, grid.meta)

  expect_equal(.computeTSMeta(meta.data1, grid.meta), meta.data2)
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
    HA1 <- rsparsematrix(nrow = 20, ncol = 20, density = 0.1)
    HB1 <- rsparsematrix(nrow = 20, ncol = 20, density = 0.1)
    HAB1 <- rsparsematrix(nrow = 20, ncol = 20, density = 0.1)

    study1 <- list(grid = cumsum(c(sample(101:200, 1),
                                   sample(50:100, size = 10, replace = TRUE))),
                   SSA = sample(1:10, 20, replace = TRUE),
                   SSB = sample(1:10, 20, replace = TRUE),
                   HAHA = (HA1 + t(HA1)) / 2,
                   HAHB = (HAB1 + t(HAB1)) / 2,
                   HBHB = (HB1 +t(HB1)) / 2,
                   HAX = sparseMatrix(i = sample(1:20, 10, replace = TRUE),
                                      j = sample(1:5, 10, replace = TRUE),
                                      x = rnorm(10), dims = c(20, 5)),
                   HBX = sparseMatrix(i = sample(1:20, 10, replace = TRUE),
                                      j = sample(1:5, 10, replace = TRUE),
                                      x = rnorm(10), dims = c(20, 5))
    )

    HA2 <- rsparsematrix(nrow = 20, ncol = 20, density = 0.1)
    HB2 <- rsparsematrix(nrow = 20, ncol = 20, density = 0.1)
    HAB2 <- rsparsematrix(nrow = 20, ncol = 20, density = 0.1)

    study2 <- list(grid = cumsum(c(sample(101:200, 1),
                                   sample(50:100, size = 10, replace = TRUE))),
                   SSA = sample(1:10, 20, replace = TRUE),
                   SSB = sample(1:10, 20, replace = TRUE),
                   HAHA = (HA2 + t(HA2)) / 2,
                   HAHB = (HAB2 + t(HAB2)) / 2,
                   HBHB = (HB2 +t(HB2)) / 2,
                   HAX = sparseMatrix(i = sample(1:20, 10, replace = TRUE),
                                      j = sample(1:5, 10, replace = TRUE),
                                      x = rnorm(10), dims = c(20, 5)),
                   HBX = sparseMatrix(i = sample(1:20, 10, replace = TRUE),
                                      j = sample(1:5, 10, replace = TRUE),
                                      x = rnorm(10), dims = c(20, 5))
    )
  })

  meta_grid <- list(study1$grid,study2$grid) |> unlist() |> unique() |> sort()

  study_SS <- list()

  study_SS <- rlist::list.append(study_SS, study1)

  study_SS <- rlist::list.append(study_SS, study2)

  meta.data1 <- lapply(study_SS, .alignStudyToMeta, meta.grid = meta_grid)

  meta.data1$meta <- list()

  meta.data2 <- meta.data1

  meta.data1 <- .computeTSMeta(meta.data1, meta_grid)

  meta.data1 <- .sumHAHA(meta.data1)
  meta.data1 <- .sumHAHB(meta.data1)
  meta.data1 <- .sumHBHB(meta.data1)
  meta.data1 <- .rmStudies(meta.data1)

  Phi <- .combineBlocks(meta.data1$meta, meta_grid)

  ee_meta <- .eigen.Matrix(Phi)

  exp.res <- .getPValues(q = meta.data1$meta$TS, lambda = ee_meta, method = "meta")

  expect_equal(Meta(study_SS[[1L]],study_SS[[2L]]), exp.res)

  # test for FBM
  study_SS2 <- study_SS
  # Convert HAHA Matrix into a matrix
  study_SS2[[1L]][["HAHA"]] <- as.matrix(study_SS2[[1L]][["HAHA"]])
  # Convert HAHA matrix into an FBM
  study_SS2[[1L]][["HAHA"]] <- as_FBM(study_SS2[[1L]][["HAHA"]])

  # Convert HAHA Matrix into a matrix
  study_SS2[[2L]][["HAHA"]] <- as.matrix(study_SS2[[2L]][["HAHA"]])
  # Convert a matrix into an FBM
  study_SS2[[2L]][["HAHA"]] <- as_FBM(study_SS2[[2L]][["HAHA"]])

  # Convert HAHB Matrix into a matrix
  study_SS2[[1L]]$HAHB <- as.matrix(study_SS2[[1L]]$HAHB)
  # Convert HAHB matrix into an FBM
  study_SS2[[1L]]$HAHB <- as_FBM(study_SS2[[1L]]$HAHB)

  # Convert HAHB Matrix into a matrix
  study_SS2[[2L]]$HAHB <- as.matrix(study_SS2[[2L]]$HAHB)
  # Convert HAHB matrix into an FBM
  study_SS2[[2L]]$HAHB <- as_FBM(study_SS2[[2L]]$HAHB)

  # Convert HBHB Matrix into a matrix
  study_SS2[[1L]][["HBHB"]] <- as.matrix(study_SS2[[1L]][["HBHB"]])
  # Convert HBHB matrix into an FBM
  study_SS2[[1L]][["HBHB"]] <- as_FBM(study_SS2[[1L]][["HBHB"]])

  # Convert HBHB Matrix into a matrix
  study_SS2[[2L]][["HBHB"]] <- as.matrix(study_SS2[[2L]][["HBHB"]])
  # Convert HBHB matrix into an FBM
  study_SS2[[2L]][["HBHB"]] <- as_FBM(study_SS2[[2L]][["HBHB"]])

  # Convert HAX Matrix into a matrix
  study_SS2[[1L]]$HAX <- as.matrix(study_SS2[[1L]]$HAX)
  # Convert HAX matrix into an FBM
  study_SS2[[1L]]$HAX <- as_FBM(study_SS2[[1L]]$HAX)
  # Convert HBX Matrix into a matrix
  study_SS2[[1L]]$HBX <- as.matrix(study_SS2[[1L]]$HBX)
  # Convert HBX matrix into an FBM
  study_SS2[[1L]]$HBX <- as_FBM(study_SS2[[1L]]$HBX)

  # Convert HAX Matrix into a matrix
  study_SS2[[2L]]$HAX <- as.matrix(study_SS2[[2L]]$HAX)
  # Convert HAX matrix into an FBM
  study_SS2[[2L]]$HAX <- as_FBM(study_SS2[[2L]]$HAX)
  # Convert HBX Matrix into a matrix
  study_SS2[[2L]]$HBX <- as.matrix(study_SS2[[2L]]$HBX)
  # Convert HBX matrix into an FBM
  study_SS2[[2L]]$HBX <- as_FBM(study_SS2[[2L]]$HBX)

  # Convert HAHA Matrix into a matrix
  meta.data2[[1L]][["HAHA"]] <- as.matrix(meta.data2[[1L]][["HAHA"]])
  # Convert HAHA matrix into an FBM
  meta.data2[[1L]][["HAHA"]] <- as_FBM(meta.data2[[1L]][["HAHA"]])

  # Convert HAHA Matrix into a matrix
  meta.data2[[2L]][["HAHA"]] <- as.matrix(meta.data2[[2L]][["HAHA"]])
  # Convert a matrix into an FBM
  meta.data2[[2L]][["HAHA"]] <- as_FBM(meta.data2[[2L]][["HAHA"]])

  # Convert HAHB Matrix into a matrix
  meta.data2[[1L]]$HAHB <- as.matrix(meta.data2[[1L]]$HAHB)
  # Convert HAHB matrix into an FBM
  meta.data2[[1L]]$HAHB <- as_FBM(meta.data2[[1L]]$HAHB)

  # Convert HAHB Matrix into a matrix
  meta.data2[[2L]]$HAHB <- as.matrix(meta.data2[[2L]]$HAHB)
  # Convert HAHB matrix into an FBM
  meta.data2[[2L]]$HAHB <- as_FBM(meta.data2[[2L]]$HAHB)

  # Convert HBHB Matrix into a matrix
  meta.data2[[1L]][["HBHB"]] <- as.matrix(meta.data2[[1L]][["HBHB"]])
  # Convert HBHB matrix into an FBM
  meta.data2[[1L]][["HBHB"]] <- as_FBM(meta.data2[[1L]][["HBHB"]])

  # Convert HBHB Matrix into a matrix
  meta.data2[[2L]][["HBHB"]] <- as.matrix(meta.data2[[2L]][["HBHB"]])
  # Convert HBHB matrix into an FBM
  meta.data2[[2L]][["HBHB"]] <- as_FBM(meta.data2[[2L]][["HBHB"]])

  # Convert HAX Matrix into a matrix
  meta.data2[[1L]]$HAX <- as.matrix(meta.data2[[1L]]$HAX)
  # Convert HAX matrix into an FBM
  meta.data2[[1L]]$HAX <- as_FBM(meta.data2[[1L]]$HAX)
  # Convert HBX Matrix into a matrix
  meta.data2[[1L]]$HBX <- as.matrix(meta.data2[[1L]]$HBX)
  # Convert HBX matrix into an FBM
  meta.data2[[1L]]$HBX <- as_FBM(meta.data2[[1L]]$HBX)

  # Convert HAX Matrix into a matrix
  meta.data2[[2L]]$HAX <- as.matrix(meta.data2[[2L]]$HAX)
  # Convert HAX matrix into an FBM
  meta.data2[[2L]]$HAX <- as_FBM(meta.data2[[2L]]$HAX)
  # Convert HBX Matrix into a matrix
  meta.data2[[2L]]$HBX <- as.matrix(meta.data2[[2L]]$HBX)
  # Convert HBX matrix into an FBM
  meta.data2[[2L]]$HBX <- as_FBM(meta.data2[[2L]]$HBX)

  meta.data2 <- .computeTSMeta(meta.data2, meta_grid)

  meta.data2 <- .sumHAHA(meta.data2)
  meta.data2 <- .sumHAHB(meta.data2)
  meta.data2 <- .sumHBHB(meta.data2)
  meta.data2 <- .rmStudies(meta.data2)

  Phi2 <- .combineBlocks(meta.data2$meta, meta_grid)

  ee_meta2 <- .eigen.FBM(Phi2)

  exp.res2 <- .getPValues(q = meta.data2$meta$TS, lambda = ee_meta2, method = "meta")

  expect_equal(Meta(study_SS2[[1L]],study_SS2[[2L]]), exp.res2)

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
