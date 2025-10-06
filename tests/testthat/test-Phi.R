test_that(".combineBlocksSmall works", {
  library(Matrix)

  withr::with_seed(12345,{
    # create a random sparse matrix 20 by 20
    #HA <- rsparsematrix(nrow = 20, ncol = 20, density = 0.1)
    #HB <- rsparsematrix(nrow = 20, ncol = 20, density = 0.1)
    #HAB <- rsparsematrix(nrow = 20, ncol = 20, density = 0.1)
    # Test with only 1 nonzero value
    HA <- sparseMatrix(i = 5, j = 5, x = 10, dims = c(20, 20))
    HB <- sparseMatrix(i = 5, j = 5, x = 10, dims = c(20, 20))
    HAB <- sparseMatrix(i = 5, j = 5, x = 10, dims = c(20, 20))

    grid = cumsum(c(sample(101:200, 1),
                    sample(50:100, size = 10, replace = TRUE)))

    meta <- list(HAHA = (HA + t(HA)) / 2,
                 HAHB = (HAB + t(HAB)) / 2,
                 HBHB = (HB +t(HB)) / 2)
  })

  G_H <- as.matrix(rep(diff(grid), 2L))
  GHGHT <- tcrossprod(sqrt(G_H))
  GHGHT <- rbind(cbind(GHGHT, GHGHT),
                 cbind(GHGHT, GHGHT))

  Phi <- rbind(cbind(meta$HAHA, meta$HAHB),
               cbind(t(meta$HAHB), meta$HBHB))

  expect.phi <- as.matrix(GHGHT * Phi)
  expect.res <- expect.phi[c(5,25),c(5,25)]

  act.res <- .combineBlocksSmall(meta, grid)

  expect_equal(as.matrix(act.res), expect.res)

})

test_that(".combineBlocksBig works",{
  library(Matrix)
  library(bigstatsr)

  withr::with_seed(12345,{
    # create a random sparse matrix 20 by 20
    HA <- rsparsematrix(nrow = 20, ncol = 20, density = 0.1)
    HB <- rsparsematrix(nrow = 20, ncol = 20, density = 0.1)
    HAB <- rsparsematrix(nrow = 20, ncol = 20, density = 0.1)

    grid = cumsum(c(sample(101:200, 1),
                    sample(50:100, size = 10, replace = TRUE)))

    meta <- list(HAHA = (HA + t(HA)) / 2,
                 HAHB = (HAB + t(HAB)) / 2,
                 HBHB = (HB +t(HB)) / 2)
  })

  expect.res <- .combineBlocksSmall(meta, grid)
  expect.mat <- as.matrix(expect.res)

  meta2 <- meta

  meta2$HAHA <- as.matrix(meta2$HAHA)
  meta2$HAHA <- as_FBM(meta2$HAHA)

  meta2$HAHB <- as.matrix(meta2$HAHB)
  meta2$HAHB <- as_FBM(meta2$HAHB)

  meta2$HBHB <- as.matrix(meta2$HBHB)
  meta2$HBHB <- as_FBM(meta2$HBHB)

  act.res <- .combineBlocksBig(meta2, grid)
  act.mat <- act.res[]

  expect_equal(act.mat, expect.mat)
  expect_s4_class(act.res, "FBM")
})

test_that(".combineBlocks works",{
  library(Matrix)
  library(bigstatsr)

  withr::with_seed(12345,{
    # create a random sparse matrix 20 by 20
    HA <- rsparsematrix(nrow = 20, ncol = 20, density = 0.1)
    HB <- rsparsematrix(nrow = 20, ncol = 20, density = 0.1)
    HAB <- rsparsematrix(nrow = 20, ncol = 20, density = 0.1)

    grid = cumsum(c(sample(101:200, 1),
                    sample(50:100, size = 10, replace = TRUE)))

    meta <- list(HAHA = (HA + t(HA)) / 2,
                 HAHB = (HAB + t(HAB)) / 2,
                 HBHB = (HB +t(HB)) / 2)
  })

  expect.res <- .combineBlocksSmall(meta, grid)
  expect.mat <- as.matrix(expect.res)

  meta2 <- meta

  meta2$HAHA <- as.matrix(meta2$HAHA)
  meta2$HAHA <- as_FBM(meta2$HAHA)

  meta2$HAHB <- as.matrix(meta2$HAHB)
  meta2$HAHB <- as_FBM(meta2$HAHB)

  meta2$HBHB <- as.matrix(meta2$HBHB)
  meta2$HBHB <- as_FBM(meta2$HBHB)

  meta3 <- meta

  meta3$HAHA <- as.matrix(meta3$HAHA)
  meta3$HAHA <- as_FBM(meta3$HAHA)

  meta3$HAHB <- as.matrix(meta3$HAHB)
  meta3$HAHB <- as_FBM(meta3$HAHB)

  meta3$HBHB <- as.matrix(meta3$HBHB)
  meta3$HBHB <- as_FBM(meta3$HBHB)

  meta4 <- meta2

  meta4$HAHA <- big_copy(meta2$HAHA)
  meta4$HAHB <- big_copy(meta2$HAHB)
  meta4$HBHB <- big_copy(meta2$HBHB)

  act.res <- .combineBlocksBig(meta2, grid)
  act.mat <- act.res[]

  act2.res <- .combineBlocks(meta4, grid)
  act2.mat <- act2.res[]

  expect_equal(as.matrix(.combineBlocks(meta, grid)), expect.mat)
  expect_equal(act2.mat, act.mat)

  expect_error(.combineBlocks(NULL, grid),
               "`blocks` must be a list of lists of covariance matrix components")

  expect_error(.combineBlocks(meta, NULL),
               "`grid` must be a numeric vector")
})
