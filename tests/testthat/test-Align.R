test_that(".dotMatch works", {
  SS <- c(1,0,0,2)
  grid.study <- c(150,200,300)
  grid.meta <- c(100,150,200,250,300,350)
  expected <- c(0,1,0,0,0,0,0,2,2,0)

  result <- .dotMatch(SS,grid.study,grid.meta)

#  expect_length(result, 2*(length(grid.meta)-1))

  expect_equal(result, expected)

  expect_error(.dotMatch(SS = NULL, grid.study = grid.study, grid.meta = grid.meta),
               "`SS` must be a numeric vector")

  expect_error(.dotMatch(SS, grid.study = NULL, grid.meta = grid.meta),
               "`grid.study` must be a numeric vector")

  expect_error(.dotMatch(SS, grid.study = grid.study, grid.meta = NULL),
               "`grid.meta` must be a numeric vector")

  expect_error(.dotMatch(SS, grid.study = c(1000, 200, 300), grid.meta = grid.meta),
               "`grid.meta` must contain `grid.study`")
})

test_that(".dotHH.Matrix works",{
  library(Matrix)

  study.HH <- rbind(cbind(Matrix(c(1,0,0,0),nrow = 2,ncol = 2),
                          Matrix(rep(0,4),nrow = 2,ncol = 2)),
                    cbind(Matrix(rep(0,4),nrow = 2,ncol = 2),
                          Matrix(c(0,0,0,2),nrow = 2,ncol = 2)))

  grid.study <- c(150,200,300)
  grid.meta <- c(100,150,200,250,300,350)

  expected <- matrix(0, nrow = 10, ncol = 10)
  expected[2,2]<-1
  expected[8:9,8:9] <- matrix(rep(2,4),nrow = 2,ncol = 2)
  #expected <- as(expected, "sparseMatrix")

  #Convert a sparse Matrix to a matrix
  result <- as.matrix(.dotHH.Matrix(study.HH,grid.study,grid.meta))

  expect_equal(result, expected)

  expect_error(.dotHH.Matrix(study.HH, NULL, grid.meta),
               "`grid.study` must be a numeric vector")
  expect_error(.dotHH.Matrix(study.HH, grid.study, NULL),
               "`grid.meta` must be a numeric vector")
  expect_error(.dotHH.Matrix(study.HH, c(50, 450), grid.meta),
               "`grid.meta` must contain `grid.study`")

})

test_that(".dotHH.FBM works",{
  library(Matrix)
  library(bigstatsr)

  withr::with_seed(12345,{
    # create a random sparse matrix 20 by 20
    HH <- rsparsematrix(nrow = 20, ncol = 20, density = 0.1)
    # create a random study grids with length 11
    grid.study <- cumsum(c(sample(101:200, 1), sample(50:100, size = 10, replace = TRUE)))

    # create a random meta grids with length 21
    # add 1 before
    pre <- sample(50:100, 1)
    # add 2 between
    insert1 <- sample((grid.study[4] + 25):(grid.study[5] - 25), 1)
    insert2 <- sample((grid.study[8] + 25):(grid.study[9] - 25), 1)

    # add 7 after
    gaps_extra <- sample(50:100, size = 7, replace = TRUE)
    grid_extra <- grid.study[length(grid.study)] + cumsum(gaps_extra)

    grid.meta <- c(pre, grid.study[1:4], insert1, grid.study[5:8], insert2,
                   grid.study[9:11], grid_extra)
  })

  # Convert a Matrix into a matrix
  study.HH <- as.matrix(HH)
  # Convert a matrix into an FBM
  study.HH <- as_FBM(study.HH)

  # expected 40 by 40 sparse matrix
  expected <- .dotHH.Matrix(HH,grid.study,grid.meta)

  result <- .dotHH.FBM(study.HH,grid.study,grid.meta)

  # Check returns a FBM object
  expect_s4_class(result, "FBM")

  # Convert FBM object into a matrix
  result <- result[]
  # COnvert a matrix into a Matrix
  result <- Matrix(result)

  expect_equal(result, expected)

  expect_error(.dotHH.FBM(study.HH, NULL, grid.meta),
               "`grid.study` must be a numeric vector")
  expect_error(.dotHH.FBM(study.HH, grid.study, NULL),
               "`grid.meta` must be a numeric vector")
  expect_error(.dotHH.FBM(study.HH, c(20, 2000), grid.meta),
               "`grid.meta` must contain `grid.study`")

})

test_that(".dotHX.Matrix works",{
  library(Matrix)

  study.HX <- Matrix(c(1,0,0,1,0,0,0,1),nrow = 4,ncol = 2)

  grid.study <- c(150,200,300)
  grid.meta <- c(100,150,200,250,300,350)

  expected <- Matrix(0,nrow = 10,ncol = 2)
  expected[2,1] <- 1
  expected[8:9,] <- Matrix(rep(1,4),nrow = 2,ncol = 2)

  result <- as(.dotHX.Matrix(study.HX,grid.study,grid.meta),"Matrix")

  expect_equal(result, expected)

  expect_error(.dotHX.Matrix(study.HX, NULL, grid.meta),
               "`grid.study` must be a numeric vector")
  expect_error(.dotHX.Matrix(study.HX, grid.study, NULL),
               "`grid.meta` must be a numeric vector")
  expect_error(.dotHX.Matrix(study.HX, c(50, 450), grid.meta),
               "`grid.meta` must contain `grid.study`")
})

test_that(".dotHX.FBM works",{
  library(Matrix)
  library(bigstatsr)

  withr::with_seed(12345,{
  # create a random sparse matrix 20 by 2
  HX <- rsparsematrix(nrow = 20, ncol = 2, density = 0.2)

  # create a random study grids with length 11
  grid.study <- cumsum(c(sample(101:200, 1), sample(50:100, size = 10, replace = TRUE)))

  # create a random meta grids with length 21
  # add 1 before
  pre <- sample(50:100, 1)
  # add 2 between
  insert1 <- sample((grid.study[4] + 25):(grid.study[5] - 25), 1)
  insert2 <- sample((grid.study[8] + 25):(grid.study[9] - 25), 1)

  # add 7 after
  gaps_extra <- sample(50:100, size = 7, replace = TRUE)

  grid_extra <- grid.study[length(grid.study)] + cumsum(gaps_extra)

  grid.meta <- c(pre, grid.study[1:4], insert1, grid.study[5:8], insert2,
                 grid.study[9:11], grid_extra)
  })

  study.HX <- as.matrix(HX)
  study.HX <- as_FBM(study.HX)

  # expected 40 by 2 sparse matrix
  expected <- .dotHX.Matrix(HX,grid.study,grid.meta)

  result <- .dotHX.FBM(study.HX,grid.study,grid.meta)

  # Returns an FBM object
  expect_s4_class(result, "FBM")

  # Convert an FBM object into a matrix
  result <- result[]
  result <- as(result,"Matrix")

  expect_equal(result, expected)

  expect_error(.dotHX.FBM(study.HX, NULL, grid.meta),
               "`grid.study` must be a numeric vector")
  expect_error(.dotHX.FBM(study.HX, grid.study, NULL),
               "`grid.meta` must be a numeric vector")
  expect_error(.dotHX.FBM(study.HX, c(20, 2000), grid.meta),
               "`grid.meta` must contain `grid.study`")
})

test_that(".alignStudyToMeta works",{
  library(Matrix)
  library(bigstatsr)

  study.data <- list(grid = c(150,200,300),
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

  grid.meta <- c(100,150,200,250,300,350)

  expected <- list(grid = study.data$grid,
                   SSA = .dotMatch(study.data$SSA, study.data$grid, grid.meta),
                   SSB = .dotMatch(study.data$SSB, study.data$grid, grid.meta),
                   HAHA = .dotHH.Matrix(study.data$HAHA, study.data$grid, grid.meta),
                   HAHB = .dotHH.Matrix(study.data$HAHB, study.data$grid, grid.meta),
                   HBHB = .dotHH.Matrix(study.data$HBHB, study.data$grid, grid.meta),
                   HAX = .dotHX.Matrix(study.data$HAX, study.data$grid, grid.meta),
                   HBX = .dotHX.Matrix(study.data$HBX, study.data$grid, grid.meta)
                   )

  result <- .alignStudyToMeta(study.data, grid.meta)
  expect_equal(result, expected)

  expect_error(.alignStudyToMeta(NULL, grid.meta),
               "`study.data` must be a list of summary statistics in Worker-Output form")

  study.data$grid<- c(50,450)
  expect_error(.alignStudyToMeta(study.data, grid.meta),
               "`meta.grid` must contain `grid.study`")
  expect_error(.alignStudyToMeta(study.data, c(100,200)),
               "`meta.grid` must contain `grid.study`")

})
