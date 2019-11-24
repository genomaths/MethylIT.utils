#' @rdname simulateCounts
#' @title Simulate read counts of methylated and unmethylated cytosines
#' @description Auxiliary function to simulate read counts of methylated and
#'     unmethylated cytosines
#'@details Methylation coverages (minimum 10) are generated from a Negative
#'     Binomial distribution with function \code{\link[MASS]{rnegbin}} from R
#'     package MASS. This function uses the representation of the Negative
#'     Binomial distribution as a continuous mixture of Poisson distributions
#'     with Gamma distributed means. Prior methylation levels are randomly
#'     generated with beta distribution using \code{\link[stats]{Beta}}
#'     function from R package “stats” and posterior methylation levels are
#'     generated according Bayes' theorem. The read of methylation counts are
#'     obtained as the product of coverage by the posterior methylation level.
#'     
#'     Counts on regions are generated from a Negative Binomial distribution
#'     with function \code{\link[MASS]{rnegbin}} with mean mu and variance:
#'      mu + mu^2/theta. 
#' @param num.samples Number of samples to generate.
#' @param sites Number of cytosine sites for each sample.
#' @param alpha Alpha parameter of beta distribution. Parameter shape1 from
#'     \code{\link[stats]{Beta}} function.
#' @param beta Beta parameter of beta distribution. Parameter shape2 from
#'     \code{\link[stats]{Beta}} function.
#' @param size number of trials (11 or more). Expected cytosine coverage.
#' @param theta Parameter theta from \code{\link[MASS]{rnegbin}}
#'     (overdispersion parameter).
#' @param sample.ids Names for the samples.
#' @param chromosome A character string naming the chromosome. Default "1".
#' @param start An nteger vector with the start positions for each cytosine 
#'     site. Default start = NULL.
#' @param end An integer vector with the end position for each cytosine site.
#'     Default end = start. Notice that if end > start, each counts will cover
#'     a region.
#' @param strand One of the characters '*', '+', or '-' denoting the DNA strand.
#'     Default is '*'.
#' @param type One of the string 'on_sites' or 'on_regions'. Default is 
#'     'on_sites'. If type == 'on_sites', then the counts are intended on 
#'     single bases, otherwise the counts are covering regions.
#' @param regions,min_width,max_width,minCountPerIndv,minCountPerIndv,mu 
#'     Arguments to provide when type == 'on_regions':
#'     
#'     \describe{
#'         \item{regions}{Number of regions carrying counts}
#'         \item{min_width}{Minimum size for a region}
#'         \item{max_width}{Maximum size for a region}
#'         \item{minCountPerIndv}{Each region must have more than 
#'             'minCountPerIndv' counts (on average) per individual
#'             (if 'mu' is not provided)}
#'         \item{maxCountPerIndv}{Each region must have less than 
#'             'maxCountPerIndv' counts (on average) per individual (if 'mu' is
#'             not provided)}
#'         \item{mu}{The expected value of the data generated with Negative
#'             Binomial distribution. This a vector of means. Short vectors are
#'             recycled. Default is NULL and, in this case, simulation is
#'             performed with mu = seq(minCountPerIndv, maxCountPerIndv).}
#'     }
#' @param noise A single number from the interval [0, 1] or a numeric vector of
#'     lengh(sites) with all its elements from the interval [0, 1]. Adds noise
#'     to the read counts. The noise is added to the methylation levels 'p',
#'     which are used to compute the coverage.  If for some site 'p' + noise >
#'     1', then the noise is not added to the site. Default is zero.
#' @param seed seed a single value, interpreted as an integer, or NULL, to
#'     set a seed for Random Number Generation. Default is seed = 123.    
#' @importFrom stats rbeta rbinom
#' @importFrom MASS rnegbin
#' @importFrom MethylIT sortBySeqnameAndStart
#' @return A list of GRanges objects with the methylated and unmethylated counts
#'     in its metacolumn.
#' @export
#' @author Robersy Sanchez (\url{https://github.com/genomaths}).
#' @examples
#' ## *** Simulate samples with expected average of difference of methylation
#' ## levels equal to 0.0427.
#' ## === Expected mean of methylation levels ===
#' bmean <- function(alpha, beta) alpha/(alpha + beta)
#' bmean(0.03, 0.5) - bmean(0.007, 0.5) #' Expected difference = 0.04279707
#'
#' ## === The number of cytosine sitesto generate ===
#' sites = 5000
#'
#' ## === Simulate samples ===
#' ref = simulateCounts(num.samples = 1, sites = sites, alpha = 0.007,
#'                     beta = 0.5, size = 50, theta = 4.5, sample.ids = "C1")
#' treat = simulateCounts(num.samples = 2, sites = sites, alpha = 0.03,
#'                     beta = 0.5, size = 50, theta = 4.5,
#'                     sample.ids = c("T1", "T2"))
#' ### === Simulate counts on regions ====
#' simulateCounts(num.samples = 7, 
#'                sample.ids = c(paste0("C",1:4), paste0("T", 1:3)),
#'                type = "on_regions", theta = 2.5, regions = 10)
#'                
simulateCounts <- function(num.samples, sites = NULL, alpha = NULL, 
                           beta = NULL, size = NULL, theta,
                           sample.ids = NULL, chromosome = "1",
                           start = NULL, end = NULL, strand = "*", 
                           type = c("on_sites", "on_regions"),
                           regions =  10, min_width = 1000,
                           max_width = 5000, minCountPerIndv = 8,
                           maxCountPerIndv = 300, mu = NULL, 
                           noise = 0, seed = 123){
   
   type <- match.arg(type)
   missed <- sapply(list(sites, alpha, beta, size), is.null)
   if (type == "on_sites" && any(missed)) 
       stop(paste0("\n *** Agument ", 
                   c("sites ", "alpha ", "beta ", "size ")[missed],
                   "must be provided"))
   
   if (length(sample.ids) != num.samples) 
       stop("*** The number of samples must be equal to length(sample.ids)")
   
   set.seed(seed)
   
   if (type == "on_sites") {
      coverage <- rnegbin(n = sites, mu = 100, theta = theta)
      chromosome <- chromosome[1] # only one chromosome per simulation
      if (is.null(start)) {start <- seq_len(sites); end <- start}
      if (is.null(end)) end <- start
      # fix sites w/o reads
      coverage <- ifelse(coverage < 10, 10, coverage)
      LR = list()
      for(k in seq_len(num.samples)) {
           p = rbeta(n = sites, shape1 = alpha, shape2 = beta)
           shape1 <- theta * p
           shape2 <- theta * (1 - p)
           p0 <- rbinom(n = sites, size=size,
                       prob = rbeta(n=sites, shape1=shape1, shape2=shape2))/size
           p <- (p0 + noise)
           idx <- which(p >  1)
           p[ idx ] <- p0[ idx ]; rm(p0, idx)
           mC = ceiling(coverage * p)
           uC = coverage - mC
           LR[[k]] <- makeGRangesFromDataFrame(data.frame(chr = chromosome, 
                                                       start = start,
                                                       end = end,
                                                       strand = strand, 
                                                       mC = mC, uC = uC), 
                                               keep.extra.columns = TRUE)
      }
      if (!is.null(sample.ids)) names(LR) <- sample.ids
   } else {
      if (is.null(start) || is.null(end)) {
         widths <- round(seq(from = min_width, to = max_width,
                          length.out = regions))
         LR <- GRanges(seqnames = chromosome, 
                       IRanges(start = seq_len(length(widths)) * widths,
                               width = widths),
                       strand = strand)
      } else LR <- makeGRangesFromDataFrame(data.frame(chr = chromosome,
                                                       start = start,
                                                       end = end, 
                                                       strand = strand))
      if (is.null(mu)) mu <- seq(minCountPerIndv, maxCountPerIndv)
      counts <- t(sapply(seq_len(length(widths)), 
                           function(...) rnegbin(n = num.samples, 
                                               mu = sample(mu, 1),
                                               theta = theta)))
      colnames(counts) <- sample.ids
      mcols(LR) <- counts
      LR <- sortBySeqnameAndStart(LR)
   }
   return(LR)
}



