library(testthat)
library(MethylIT.utils)
context("MethylIT.utils unlistX tests")

test_that("unlistX function test", {
   gr1 <- GRanges(seqnames = "chr2", ranges = IRanges(3, 6),
               strand = "+", score = 5L, GC = 0.45)
   gr2 <- GRanges(seqnames = c("chr1", "chr1"),
               ranges = IRanges(c(7,13), width = 3),
               strand = c("+", "-"), score = 3:4, GC = c(0.3, 0.5))
   gr3 <- GRanges(seqnames = c("chr1", "chr2"),
               ranges = IRanges(c(1, 4), c(3, 9)),
               strand = c("-", "-"), score = c(6L, 2L), GC = c(0.4, 0.1))
   grl <- list("gr1" = gr1, "gr2" = gr2, "gr3" = gr3)
   class(grl) <-'InfDiv' # A trick
   grl <- unlist(grl)
   expect_true(length(grl) == 5)
})
