test_that(".sumHxHx.Matrix works", {
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

  meta.data2$meta[["HAHA"]] <- meta.data2[[1L]][["HAHA"]] -
    as(tcrossprod(meta.data2[[1L]][["HAX"]]), "symmetricMatrix") +
    meta.data2[[2L]][["HAHA"]] -
    as(tcrossprod(meta.data2[[2L]][["HAX"]]), "symmetricMatrix")

  meta.data2[[1L]][["HAHA"]] <- NULL
  meta.data2[[2L]][["HAHA"]] <- NULL

  expect_equal(.sumHxHx.Matrix(meta.data1, "A"), meta.data2)

  # test for B
  meta.data3 <- meta.data1

  meta.data3$meta[["HBHB"]] <- meta.data3[[1L]][["HBHB"]] -
    as(tcrossprod(meta.data3[[1L]][["HBX"]]), "symmetricMatrix") +
    meta.data3[[2L]][["HBHB"]] -
    as(tcrossprod(meta.data3[[2L]][["HBX"]]), "symmetricMatrix")

  meta.data3[[1L]][["HBHB"]] <- NULL
  meta.data3[[2L]][["HBHB"]] <- NULL

  expect_equal(.sumHxHx.Matrix(meta.data1, "B"), meta.data3)

})

test_that(".sumHxHy.Matrix works", {
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

  meta.data2$meta$HAHB <- meta.data2[[1L]]$HAHB -
    tcrossprod(meta.data2[[1L]]$HAX,meta.data2[[1L]]$HBX) +
    meta.data2[[2L]]$HAHB -
    tcrossprod(meta.data2[[2L]]$HAX,meta.data2[[2L]]$HBX)

  meta.data2[[1L]]$HAHB <- NULL
  meta.data2[[2L]]$HAHB <- NULL

  expect_equal(.sumHxHy.Matrix(meta.data1), meta.data2)
})

test_that(".sumHxHx.FBM works",{
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

  # Convert HAHA Matrix into a matrix
  meta.data2[[1L]][["HAHA"]] <- as.matrix(meta.data2[[1L]][["HAHA"]])
  # Convert HAHA matrix into an FBM
  meta.data2[[1L]][["HAHA"]] <- as_FBM(meta.data2[[1L]][["HAHA"]])
  # Convert HAX Matrix into a matrix
  meta.data2[[1L]][["HAX"]] <- as.matrix(meta.data2[[1L]][["HAX"]])
  # Convert HAX matrix into an FBM
  meta.data2[[1L]][["HAX"]] <- as_FBM(meta.data2[[1L]][["HAX"]])

  # Convert HAHA Matrix into a matrix
  meta.data2[[2L]][["HAHA"]] <- as.matrix(meta.data2[[2L]][["HAHA"]])
  # Convert a matrix into an FBM
  meta.data2[[2L]][["HAHA"]] <- as_FBM(meta.data2[[2L]][["HAHA"]])

  # Convert HAX Matrix into a matrix
  meta.data2[[2L]][["HAX"]] <- as.matrix(meta.data2[[2L]][["HAX"]])
  # Convert HAX matrix into an FBM
  meta.data2[[2L]][["HAX"]] <- as_FBM(meta.data2[[2L]][["HAX"]])

  Matrix.res <- .sumHxHx.Matrix(meta.data1, "A")
  Matrix.res$meta[["HAHA"]] <- as.matrix(Matrix.res$meta[["HAHA"]])

  act.res <- .sumHxHx.FBM(meta.data2, "A")
  expect_s4_class(act.res$meta[["HAHA"]], "FBM")
  act.res$meta[["HAHA"]] <- act.res$meta[["HAHA"]][]

  # Test for "A"
  expect_equal(act.res$meta, Matrix.res$meta)

  # Convert HBHB Matrix into a matrix
  meta.data2[[1L]][["HBHB"]] <- as.matrix(meta.data2[[1L]][["HBHB"]])
  # Convert HBHB matrix into an FBM
  meta.data2[[1L]][["HBHB"]] <- as_FBM(meta.data2[[1L]][["HBHB"]])
  # Convert HBX Matrix into a matrix
  meta.data2[[1L]][["HBX"]] <- as.matrix(meta.data2[[1L]][["HBX"]])
  # Convert HBX matrix into an FBM
  meta.data2[[1L]][["HBX"]] <- as_FBM(meta.data2[[1L]][["HBX"]])

  # Convert HBHB Matrix into a matrix
  meta.data2[[2L]][["HBHB"]] <- as.matrix(meta.data2[[2L]][["HBHB"]])
  # Convert HBHB matrix into an FBM
  meta.data2[[2L]][["HBHB"]] <- as_FBM(meta.data2[[2L]][["HBHB"]])

  # Convert HBX Matrix into a matrix
  meta.data2[[2L]][["HBX"]] <- as.matrix(meta.data2[[2L]][["HBX"]])
  # Convert HBX matrix into an FBM
  meta.data2[[2L]][["HBX"]] <- as_FBM(meta.data2[[2L]][["HBX"]])

  Matrix.res <- .sumHxHx.Matrix(meta.data1, "B")
  Matrix.res$meta[["HBHB"]] <- as.matrix(Matrix.res$meta[["HBHB"]])

  act.res <- .sumHxHx.FBM(meta.data2, "B")
  expect_s4_class(act.res$meta[["HBHB"]], "FBM")
  act.res$meta[["HBHB"]] <- act.res$meta[["HBHB"]][]

  # Test for "B"
  expect_equal(act.res$meta, Matrix.res$meta)
})


test_that(".sumHxHy.FBM works",{
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

  # Convert HBHB Matrix into a matrix
  meta.data2[[1L]]$HAHB <- as.matrix(meta.data2[[1L]]$HAHB)
  # Convert HBHB matrix into an FBM
  meta.data2[[1L]]$HAHB <- as_FBM(meta.data2[[1L]]$HAHB)

  # Convert HBHB Matrix into a matrix
  meta.data2[[2L]]$HAHB <- as.matrix(meta.data2[[2L]]$HAHB)
  # Convert HBHB matrix into an FBM
  meta.data2[[2L]]$HAHB <- as_FBM(meta.data2[[2L]]$HAHB)

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

  Matrix.res <- .sumHxHy.Matrix(meta.data1)
  Matrix.res$meta$HAHB <- as.matrix(Matrix.res$meta$HAHB)

  act.res <- .sumHxHy.FBM(meta.data2)
  expect_s4_class(act.res$meta$HAHB, "FBM")
  act.res$meta$HAHB <- act.res$meta$HAHB[]

  expect_equal(act.res$meta, Matrix.res$meta)
})

test_that(".sumHAHA works",{
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

  expect_equal(.sumHAHA(meta.data1), .sumHxHx.Matrix(meta.data1, "A"))

  # Convert HAHA Matrix into a matrix
  meta.data2[[1L]][["HAHA"]] <- as.matrix(meta.data2[[1L]][["HAHA"]])
  # Convert HAHA matrix into an FBM
  meta.data2[[1L]][["HAHA"]] <- as_FBM(meta.data2[[1L]][["HAHA"]])
  # Convert HAX Matrix into a matrix
  meta.data2[[1L]][["HAX"]] <- as.matrix(meta.data2[[1L]][["HAX"]])
  # Convert HAX matrix into an FBM
  meta.data2[[1L]][["HAX"]] <- as_FBM(meta.data2[[1L]][["HAX"]])

  # Convert HAHA Matrix into a matrix
  meta.data2[[2L]][["HAHA"]] <- as.matrix(meta.data2[[2L]][["HAHA"]])
  # Convert a matrix into an FBM
  meta.data2[[2L]][["HAHA"]] <- as_FBM(meta.data2[[2L]][["HAHA"]])

  # Convert HAX Matrix into a matrix
  meta.data2[[2L]][["HAX"]] <- as.matrix(meta.data2[[2L]][["HAX"]])
  # Convert HAX matrix into an FBM
  meta.data2[[2L]][["HAX"]] <- as_FBM(meta.data2[[2L]][["HAX"]])

  act.res <- .sumHAHA(meta.data2)
  expect_s4_class(act.res$meta[["HAHA"]], "FBM")
  act.res$meta[["HAHA"]] <- act.res$meta[["HAHA"]][]

  exp.res <- .sumHxHx.FBM(meta.data2, "A")
  exp.res$meta[["HAHA"]] <- exp.res$meta[["HAHA"]][]

  expect_equal(act.res$meta, exp.res$meta)

})

test_that(".sumHAHB works",{
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
  expect_equal(.sumHAHB(meta.data1), .sumHxHy.Matrix(meta.data1))

  # Convert HBHB Matrix into a matrix
  meta.data2[[1L]]$HAHB <- as.matrix(meta.data2[[1L]]$HAHB)
  # Convert HBHB matrix into an FBM
  meta.data2[[1L]]$HAHB <- as_FBM(meta.data2[[1L]]$HAHB)

  # Convert HBHB Matrix into a matrix
  meta.data2[[2L]]$HAHB <- as.matrix(meta.data2[[2L]]$HAHB)
  # Convert HBHB matrix into an FBM
  meta.data2[[2L]]$HAHB <- as_FBM(meta.data2[[2L]]$HAHB)

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

  act.res <- .sumHAHB(meta.data2)
  expect_s4_class(act.res$meta$HAHB, "FBM")
  act.res$meta$HAHB <- act.res$meta$HAHB[]

  exp.res <- .sumHxHy.FBM(meta.data2)
  exp.res$meta$HAHB <- exp.res$meta$HAHB[]

  expect_equal(act.res$meta, exp.res$meta)

})

test_that(".sumHBHB works",{
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

  expect_equal(.sumHAHA(meta.data1), .sumHxHx.Matrix(meta.data1, "A"))

  # Convert HBHB Matrix into a matrix
  meta.data2[[1L]][["HBHB"]] <- as.matrix(meta.data2[[1L]][["HBHB"]])
  # Convert HBHB matrix into an FBM
  meta.data2[[1L]][["HBHB"]] <- as_FBM(meta.data2[[1L]][["HBHB"]])
  # Convert HBX Matrix into a matrix
  meta.data2[[1L]][["HBX"]] <- as.matrix(meta.data2[[1L]][["HBX"]])
  # Convert HBX matrix into an FBM
  meta.data2[[1L]][["HBX"]] <- as_FBM(meta.data2[[1L]][["HBX"]])

  # Convert HBHB Matrix into a matrix
  meta.data2[[2L]][["HBHB"]] <- as.matrix(meta.data2[[2L]][["HBHB"]])
  # Convert HBHB matrix into an FBM
  meta.data2[[2L]][["HBHB"]] <- as_FBM(meta.data2[[2L]][["HBHB"]])

  # Convert HBX Matrix into a matrix
  meta.data2[[2L]][["HBX"]] <- as.matrix(meta.data2[[2L]][["HBX"]])
  # Convert HBX matrix into an FBM
  meta.data2[[2L]][["HBX"]] <- as_FBM(meta.data2[[2L]][["HBX"]])

  act.res <- .sumHBHB(meta.data2)
  expect_s4_class(act.res$meta[["HBHB"]], "FBM")
  act.res$meta[["HBHB"]] <- act.res$meta[["HBHB"]][]

  exp.res <- .sumHxHx.FBM(meta.data2, "B")
  exp.res$meta[["HBHB"]] <- exp.res$meta[["HBHB"]][]

  expect_equal(act.res$meta, exp.res$meta)

})

test_that(".rmStudies works",{
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
  meta.data2[c(1L:2L)] <- NULL

  expect_equal(.rmStudies(meta.data1), meta.data2)
})
