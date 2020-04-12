library(testthat)
library(MethylIT.utils)

context("countSignal tests")

test_that("countSignal test", {
    some.signal <- makeGRangesFromDataFrame(data.frame(chr = "chr1",
                                                      start = 1:15,
                                                      end = 1:15,
                                                      strand = '*'))
    
    some.regions <- makeGRangesFromDataFrame(data.frame(chr = "chr1",
                                                      start = c(2, 8),
                                                      end = c(7, 14),
                                                      strand = '*'))
    
    y <- countSignal(signal = some.signal, gr = some.regions)
    expect_equal(y$sites, c(6, 7))
})
