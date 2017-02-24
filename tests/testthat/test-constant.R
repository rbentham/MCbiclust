
data(CCLE_small)
data(Mitochondrial_genes)

mito.loc <- which(row.names(CCLE_small) %in% Mitochondrial_genes)
CCLE.mito <- CCLE_small[mito.loc,]

set.seed(101)
CCLE.seed <- FindSeed(gem = CCLE.mito,seed.size = 10,
                      iterations = 100,messages = 50)

CCLE.sort <- SampleSort(gem = CCLE.mito,seed = CCLE.seed,sort.length = 11)

CCLE.samp.sort <- MCbiclust:::Vignette_sort[[1]]

CCLE.pc1 <- PC1VecFun(top.gem = CCLE.mito,
                      seed.sort = CCLE.samp.sort,
                      n = 10)

test_that("HclustGenesHiCor and CVEval work with matrices",{
  expect_error(as.numeric(HclustGenesHiCor(as.matrix(CCLE.mito),
                                                CCLE.seed,
                                                cuts = 58)), NA)
  
  expect_error(CVEval(gem.part = as.matrix(CCLE.mito),
         gem.all = as.matrix(CCLE_small),
         seed = CCLE.seed,
         splits = 58), NA)
  
})

test_that("PC1VecFun returns as expected",{
  expect_equal(length(CCLE.pc1), 500)
})

CCLE.cor.vec <-  CVEval(gem.part = CCLE.mito,
                        gem.all = CCLE_small,
                        seed = CCLE.seed,
                        splits = 10)

CCLE.bic <- ThresholdBic(cor.vec = CCLE.cor.vec,sort.order = CCLE.samp.sort,
                         pc1 = as.numeric(CCLE.pc1))

CCLE.pc1 <- PC1Align(gem = CCLE_small, pc1 = CCLE.pc1,
                     cor.vec = CCLE.cor.vec ,
                     sort.order = CCLE.samp.sort,
                     bic =CCLE.bic)

CCLE.fork <- ForkClassifier(CCLE.pc1, samp.num = length(CCLE.bic[[2]]))

test_that("Basic functioning remains constant",{
  expect_equal_to_reference(CCLE.seed,"ccle_seed.rds")
  expect_equal_to_reference(CCLE.sort,"ccle_sort.rds")
  expect_equal_to_reference(CCLE.cor.vec,"ccle_cv.rds")
  expect_equal_to_reference(CCLE.bic,"ccle_bic.rds")
  #expect_equal_to_reference(CCLE.pc1,"ccle_pc1.rds")
  expect_equal_to_reference(CCLE.fork,"ccle_fork.rds")
})


test_that("FindSeed returns as expected",{
  expect_equal(length(CCLE.seed), 10)
  expect_message(FindSeed(gem = CCLE.mito,seed.size = 10,
                          iterations = 100,messages = 50),"Iteration")
  expect_error(FindSeed(gem = CCLE.mito,seed.size = 501,
                        iterations = 100,messages = 50))
})

test_that("SampleSort returns as expected",{
  expect_equal(length(CCLE.sort), 11)
})

test_that("CVEval returns as expected",{
  expect_equal(length(CCLE.cor.vec), 1000)
})

test_that("ThresholdBic returns as expected",{
  expect_equal(length(CCLE.bic), 2)
  expect_equal(length(CCLE.bic[[1]]), 1000)
  expect_equal(length(CCLE.bic[[2]]), 58)
})

test_that("PC1Align returns as expected",{
  expect_equal(length(CCLE.pc1), 500)
})


CCLE.bic[[2]] <- CCLE.bic[[2]][1]

test_that("PC1Align works with bicluster with single sample",{
  expect_error(PC1Align(gem = CCLE_small, pc1 = CCLE.pc1,
                        cor.vec = CCLE.cor.vec ,
                        sort.order = CCLE.samp.sort,
                        bic =CCLE.bic), NA)
          }
          )

          