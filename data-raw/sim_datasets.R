# ============================================================================ #
#
# ======== Script used to generate the datasets used in the examples ========= #
#
# ============================================================================ #

### -------------------------------------------------------------------------- #
## ------- GRanges object with hypothetical signal matrix in the metacolums
### -------------------------------------------------------------------------- #

set.seed(123)
## An auxiliary function to generate simulated hypothetical values from a
## variable with normal distribution

hypDT <- function(mean, sd, n, num.pos, noise = 20) {
    h <- hist(rnorm(n, mean = mean, sd = sd), breaks = num.pos, plot = FALSE)
    hyp <- h$density * 60 + runif(length(h$density)) * noise
    return(hyp)
}

mean <- 12
sd <- 2

## To add some noise
noise <- c(4, 10)
noise2 <- list(c(5, 5), c(6, 6))

## To generate a matrix of values with variations introduced by noise
hyp <- lapply(1:2, function(k) {
     h <- hypDT(mean = mean, sd = sd, n = 10^5,
                num.pos = 8000, noise = noise[k])
    h1 <- h + runif(length(h)) * noise2[[k]][1]
    h2 <- h + runif(length(h)) * noise2[[k]][2]
    h <- h + runif(length(h)) * noise2[[k]][1]
    return(cbind(h, h1, h2))
})

## A GRanges object is built, which will carries the previous matrix on its
## meta-columns
min.length <- min(unlist(lapply(hyp, nrow)))
hyp <- lapply(hyp, function(h) h[1:min.length,])
hyp <- do.call(cbind, hyp)
starts <- seq(0, 30000, 3)[1:min.length]
ends <- starts + 2
GRMatrix <- GRanges(seqnames = "chr1", ranges = IRanges(start = starts,
                end = ends))
mcols(GRMatrix) <- data.frame(hyp = hyp)
colnames(mcols(GRMatrix)) <- c("CT1", "CT2", "CT3", "TT1", "TT2", "TT3")

usethis::use_data(GRMatrix, overwrite = TRUE )

