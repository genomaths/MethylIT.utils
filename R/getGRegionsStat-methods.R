#' @name getGRegionsStat-methods
#' @rdname getGRegionsStat-methods
#' @title Statistic of Genomic Regions
#' @description A function to estimate summarized measures of a specified
#'     variable given in a GRanges object (a column from the metacolums of the
#'     GRanges object) after split the GRanges object into intervals. A faster
#'     alternative would be \code{\link{getGRegionsStat2}}.
#' @details This function split a Grange object into intervals genomic regions
#'     (GR) of fixed size (as given in function 'tileMethylCounts2' R package
#'     methylKit, with small changes). A summarized statistic (mean, median,
#'     geometric mean or sum) is calculated for the specified variable values
#'     from each region. Notice that if win.size == step.size, then
#'     non-overlapping windows are obtained.
#' @param GR A GRange object or a list of GRanges object with the variable of 
#'     interest in the GRanges metacolumn.
#' @param win.size An integer for the size of the windows/regions size of the
#'     intervals of genomics regions.
#' @param step.size Interval at which the regions/windows must be defined
#' @param grfeatures A GRanges object corresponding to an annotated genomic
#'     feature. For example, gene region, transposable elements, exons,
#'     intergenic region, etc. If provided, then parameters 'win.size' and
#'     step.size are ignored and the statistics are estimated for 'grfeatures'.
#' @param stat Statistic used to estimate the summarized value of the variable
#'     of interest in each interval/window. Posible options are: 
#' \describe{
#'   \item{\strong{'mean':}}{The mean of values inside each region.}
#'   \item{\strong{'gmean':}}{The geometric mean of values inside each region.}
#'   \item{\strong{'median':}}{The median of values inside each region.}
#'   \item{\strong{'density':}}{The density of values inside each region. That
#'          is, the sum of values found in each region divided by the width of 
#'          the region.}
#'   \item{\strong{'count':}}{Compute the number/count of positions with values
#'          greater than zero inside each regions.}
#'   \item{\strong{'denCount':}}{The number of sites with value > 0 inside each 
#'          region divided by the width of the region.}
#'   \item{\strong{'sum':}}{The sum of values inside each region.}
#' }
#' @param absolute Optional. Logic (default: FALSE). Whether to use the absolute
#'     values of the variable provided. For example, the difference of
#'     methylation levels could take negative values (TV) and we would be
#'     interested on the sum of abs(TV), which is sum of the total variation
#'     distance.
#' @param select.strand Optional. If provided,'+' or '-', then the summarized
#'     statistic is computed only for the specified DNA chain.
#' @param column Integer number denoting the column where the variable of
#'     interest is located in the metacolumn of the GRanges object.
#' @param prob Logic. If TRUE and the variable of interest has values between
#'     zero and 1, then the summarized statistic is comuputed using Fisher's
#'     transformation.
#' @param entropy Logic. Whether to compute the entropy at each site from the 
#'     specified regions. All the values from the selected column must belong to
#'     the interval [0, 1]. Next, the requested statistics for the entropy 
#'     values at each site inside the regions is computed.
#' @param maxgap,minoverlap,type See ?findOverlaps in the IRanges package for a
#'     description of these arguments.
#' @param ignore.strand When set to TRUE, the strand information is ignored in
#'     the overlap calculations.
#' @param scaling integer (default 1). Scaling factor to be used when
#'     stat = 'density'. For example, if scaling = 1000, then density * scaling
#'     denotes the sum of values in 1000 bp.
#' @param logbase A positive number: the base with respect to which logarithms
#'     are computed when parameter 'entropy = TRUE' (default: logbase = 2).
#' @param missings Whether to write '0' or 'NA' on regions where there is not
#'     data to compute the statistic.
#' @param maxgap,minoverlap,type,select,ignore.strand Used to find overlapped 
#'     regions. See \code{\link[IRanges]{findOverlaps-methods}} in the 
#'     \strong{IRanges} package for a description of these arguments.
#' @param na.rm Logical value. If TRUE, the NA values will be removed
#' @param naming Logical value. If TRUE, the rows GRanges object will be 
#'     given the names(GR). Default is FALSE.
#' @param output A string. Setting output = 'all' will return all the regions
#'     given in 'grfeatures'. Default is output = 'hits'.
#' @param num.cores,tasks Paramaters for parallele computation using package
#'     \code{\link[BiocParallel]{BiocParallel-package}}: the number of cores to
#'     use, i.e. at most how many child processes will be run simultaneously
#'     (see \code{\link[BiocParallel]{bplapply}} and the number of tasks per job
#'     (only for Linux OS).
#' @param verbose Logical. Default is TRUE. If TRUE, then the progress of the
#'     computational tasks is given.
#' @return An object of the same class of \emph{GR} with the new genomic regions
#'     and their corresponding summarized statistic.
#' @examples
#' library(GenomicRanges)
#' gr <- GRanges(seqnames = Rle( c('chr1', 'chr2', 'chr3', 'chr4'),
#'             c(5, 5, 5, 5)),
#'             ranges = IRanges(start = 1:20, end = 1:20),
#'             strand = rep(c('+', '-'), 10),
#'             GC = seq(1, 0, length = 20))
#' grs <- getGRegionsStat(gr, win.size = 4, step.size = 4)
#' grs
#'
#' ## Selecting the positive strand
#' grs <- getGRegionsStat(gr, win.size = 4, step.size = 4, select.strand = '+')
#' grs
#'
#' ## Selecting the negative strand
#' grs <- getGRegionsStat(gr, win.size = 4, step.size = 4, select.strand = '-')
#' grs
#'
#' ## Operating over a list of GRanges objects
#' gr2 <- GRanges(seqnames = Rle( c('chr1', 'chr2', 'chr3', 'chr4'),
#'                             c(5, 5, 5, 5)),
#'                 ranges = IRanges(start = 1:20, end = 1:20),
#'                 strand = rep(c('+', '-'), 10),
#'                 GC = runif(20))
#'
#' grs <- getGRegionsStat(list(gr1 = gr, gr2 = gr2), win.size = 4, step.size=4)
#' 
#' ## Compute the density of entropy inside each region
#' gr$GC <- runif(20) 
#' getGRegionsStat(gr, win.size = 4, step.size = 4, entropy = TRUE, 
#'                 stat = "density")
#'                 
#' @importFrom GenomeInfoDb seqnames seqlengths
#' @importFrom GenomicRanges GRanges GRangesList findOverlaps
#' @importFrom IRanges IRanges
#' @importFrom data.table data.table
#' @importFrom stats median
#' @importFrom S4Vectors subjectHits queryHits DataFrame mcols 
#' @importFrom S4Vectors mcols<-
#' @importFrom MethylIT sortBySeqnameAndStart
#' @importFrom BiocGenerics strand start end
#' @importFrom BiocParallel MulticoreParam SnowParam bplapply bpstart
#' @seealso \code{\link{getGRegionsStat2}}.
#' @export
#' @author Robersy Sanchez (\url{https://github.com/genomaths}).
#'
#' @aliases getGRegionsStat
setGeneric("getGRegionsStat", 
            function(
                  GR, 
                  win.size = 1,
                  step.size = 1,
                  grfeatures = NULL, 
                  stat = c("sum", "mean", "gmean", "median",
                            "density", "count", "denCount"),
                  absolute = FALSE, 
                  select.strand = NULL, 
                  column = 1L, 
                  prob = FALSE, 
                  entropy = FALSE, 
                  maxgap = -1L, 
                  minoverlap = 0L,
                  scaling = 1000L, 
                  logbase = 2, 
                  missings = 0,
                  type = c("any", "start", "end", "within", "equal"),
                  ignore.strand = FALSE, 
                  na.rm = TRUE,  
                  naming = FALSE, 
                  output = c("hits", "all"), 
                  num.cores = 1L, 
                  tasks = 0,
                  verbose = TRUE, ...) standardGeneric("getGRegionsStat"))

#' @aliases getGRegionsStat
setMethod("getGRegionsStat", signature(GR = "GRanges"), function(GR, 
    win.size = 350, step.size = 350, grfeatures = NULL, stat = c("sum", 
        "mean", "gmean", "median", "density", "count", "denCount"), 
    absolute = FALSE, select.strand = NULL, column = 1L, prob = FALSE, 
    entropy = FALSE, maxgap = -1L, minoverlap = 0L, scaling = 1000L, 
    logbase = 2, missings = 0, type = c("any", "start", "end", "within", 
        "equal"), ignore.strand = FALSE, na.rm = TRUE, naming = FALSE, 
    output = c("hits", "all")) {
   
    ## These NULL quiet: no visible binding for global variable 'x2'
    ent <- statistic <- NULL
    if (!inherits(GR, "GRanges")) 
       stop("GR object must inherits from a GRanges class")
    if (!is.null(grfeatures) && !inherits(grfeatures, "GRanges")) {
        stop("* 'grfeatures', if provided, must be a GRanges object")
    }
    stat <- match.arg(stat)
    output <- match.arg(output)
    if (length(column) > 1) 
        stop("*** Argument of 'colum' must be an integer number")
    
    if (!is.element(missings, c(0, NA))) missings <- NA
    
    type <- match.arg(type, c("any", "start", "end", "within", "equal"))
    
    ## === Some functions to use ===
    statist <- function(x, stat = c()) {
        x <- switch(stat, 
                    count = sum(x != 0, na.rm = na.rm), 
                    sum = sum(x, na.rm = na.rm), 
                    mean = mean(x, na.rm = na.rm), 
                    gmean = exp(sum(log(x[x > 0]), na.rm = na.rm)/length(x)),
                    median = median(x, na.rm = na.rm), 
                    density = sum(x, na.rm = na.rm), 
                    denCount = sum(x != 0, na.rm = na.rm), entropy = x)
    }
    
    ## =============================== ##
    if (!is.null(select.strand)) {
        ## possible values '-', '+', NULL
        if (!is.element(select.strand, unique(strand(GR)))) {
            stop("* The GRanges object does not have strand named ", 
                "'", select.strand, "'")
        }
        idx <- which(as.character(strand(GR)) == select.strand)
        GR <- GR[idx]
    }
    
    GR <- GR[, column]
    if (absolute) 
        mcols(GR) <- data.frame(abs(as.matrix(mcols(GR))))
    
    chrs <- as.character(unique(seqnames(GR)))
    
    ## === If genomic features are not specified ===
    if (is.null(grfeatures)) {
        if (length(GR) < win.size || length(GR) < step.size) {
            stop("* 'GR'length is lesser of 'win.size' or 'step.size'")
        }
        grfeatures <- GRanges()
        for (k in seq_along(chrs)) {
            ## get max length of chromosome
            max.length <- max(IRanges::end(GR[seqnames(GR) == chrs[k], ]))
            
            ## get start chromosome coordinate
            start0 <- min(IRanges::start(GR[seqnames(GR) == chrs[k], ]))
            
            ## get sliding windows
            numTiles <- floor((max.length - (win.size - step.size))/step.size) + 
                1
            ranges <- IRanges(start = (start0 + 0:(numTiles - 1) * step.size), 
                width = rep(win.size, numTiles))
            
            ends <- IRanges::end(ranges)
            if (max(ends) > max.length) {
                idx <- (which(ends <= max.length))
                idx <- c(idx, which.max(idx) + 1)
                ranges <- ranges[idx]
            }
            
            temp.wins <- GRanges(seqnames = rep(chrs[k], length(ranges)), 
                                ranges = ranges)
            grfeatures <- suppressWarnings(c(grfeatures, temp.wins))
        }
        
        ## sites of interest inside of the windows
        if (output == "all") 
            GR0 <- grfeatures
        Hits <- findOverlaps(GR, grfeatures, maxgap = maxgap, 
                            minoverlap = minoverlap, 
                            ignore.strand = ignore.strand, type = type)
        if (length(Hits) > 0) {
            m <- ncol(mcols(GR))
            if (m > 1) {
                mcols(grfeatures) <- matrix(missings, nrow = length(grfeatures), 
                  ncol = m)
            } else mcols(grfeatures) <- missings
            grfeatures <- grfeatures[subjectHits(Hits)]
            GR <- GR[queryHits(Hits)]
            mcols(grfeatures) <- mcols(GR)
            chr <- seqnames(grfeatures)
            
            ## Variable to mark the limits of each GR
            text <- paste(chr, start(grfeatures), end(grfeatures), 
                sep = "_")
            cluster.id <- data.frame(cluster.id = text)
            GR <- grfeatures
            rm(text, grfeatures)
            gc()
            colnames(mcols(GR)) <- "statistic"
        } else {
            GR <- grfeatures
            mcols(GR) <- missings
            colnames(mcols(GR)) <- "statistic"
        }
    } else {
        ## sites of interest inside of the windows
        if (output == "all") 
            GR0 <- grfeatures
        Hits <- findOverlaps(GR, grfeatures, maxgap = maxgap,
                            minoverlap = minoverlap, 
                            ignore.strand = ignore.strand, type = type)
        if (length(Hits) > 0) {
            m <- ncol(mcols(GR))
            if (m > 1) {
                mcols(grfeatures) <- matrix(missings, nrow = length(grfeatures), 
                  ncol = m)
            } else mcols(grfeatures) <- missings
            grfeatures <- grfeatures[subjectHits(Hits)]
            GR <- GR[queryHits(Hits)]
            mcols(grfeatures) <- mcols(GR)
            
            chr <- seqnames(grfeatures)
            if (class(names(grfeatures)) == "character") {
                cluster.id <- data.frame(cluster.id = names(grfeatures))
            } else {
                text <- paste(chr, start(grfeatures), end(grfeatures), 
                  strand(grfeatures), sep = "_")
                cluster.id <- data.frame(cluster.id = text)
                rm(text)
            }
            GR <- grfeatures
            rm(grfeatures)
            gc()
            colnames(mcols(GR)) <- "statistic"
        } else {
            GR <- grfeatures
            mcols(GR) <- missings
            colnames(mcols(GR)) <- "statistic"
        }
    }
    
    if (abs(sum(GR$statistic, na.rm = TRUE)) > 0) {
        mcols(GR) <- DataFrame(cluster.id, mcols(GR))
        GR <- data.table(as.data.frame(GR))
        if (length(column) < 2) {
            colnames(GR) <- c("seqnames", "start", "end", "width", 
                "strand", "cluster.id", "statistic")
        }
        
        if (prob || entropy) {
            ## Apply Fisher transformation
            cond1 <- all(GR$statistic < 1)
            cond2 <- all(GR$statistic > 0)
            if (!cond1 && !cond2) 
                stop("*** \nAll the values from a probability vector must", 
                  "belong to the interval [0,1]")
            GR$statistic <- atanh(GR$statistic)
        }
        grn <- c("seqnames", "start", "end")
        ## =========== Compute statistic for regions =====================
        if (!entropy) {
            GR <- GR[, list(seqnames = unique(seqnames), start = min(start), 
                end = max(end), statistic = statist(statistic, stat)), 
                by = cluster.id]
            cluster.id <- GR$cluster.id
            GR <- data.frame(GR)[, c(grn, "statistic")]
        }
        
        if (entropy) {
            GR$statistic <- shannonEntr(GR$statistic, logbase = logbase)
            GR <- GR[, list(seqnames = unique(seqnames), start = min(start),
                           end = max(end), 
                           statistic = statist(statistic, stat)), 
                        by = cluster.id]
            cluster.id <- GR$cluster.id
            GR <- data.frame(GR)[, c(grn, "statistic")]
        }
        
        
        if (!is.null(select.strand)) 
            GR$strand <- select.strand
        GR <- makeGRangesFromDataFrame(GR, keep.extra.columns = TRUE)
        if (stat == "density" && !prob && !entropy) {
            widths <- width(GR)
            GR$statistic <- (scaling * GR$statistic/widths)
        }
        if (stat == "denCount" && !prob && !entropy) {
            widths <- width(GR)
            GR$statistic <- (scaling * GR$statistic/widths)
        }
        if (!is.na(missings)) {
            idx <- is.na(GR$statistic)
            if (any(idx)) 
                GR$statistic[idx] <- 0
        }
    } else cluster.id <- NULL
    if (naming) 
        names(GR) <- cluster.id
    if (output == "all") {
        mcols(GR0) <- integer(length(GR0))
        colnames(mcols(GR0)) <- "statistic"
        if (naming && is.null(names(GR0))) 
            names(GR0) <- paste(chr, start(GR0), end(GR0), sep = "_")
        GR0 <- GR0[-subjectHits(Hits)]
        GR0 <- unique(GR0)
        if (ignore.strand) 
            strand(GR0) <- "*"
        GR <- c(GR, GR0)
    }
    return(GR)
})

# ==================== Function to operate on lists ====================== #
getGRegionsStats <- function(GR, win.size = 350, step.size = 350, 
    grfeatures = NULL, stat = c("sum", "mean", "gmean", "median", 
        "density", "count", "denCount"), absolute = FALSE, select.strand = NULL, 
    column = 1L, prob = FALSE, entropy = FALSE, maxgap = -1L, minoverlap = 0L, 
    scaling = 1000L, logbase = 2, missings = 0, type = c("any", "start", 
        "end", "within", "equal"), ignore.strand = FALSE, na.rm = TRUE, 
    naming = FALSE, output = c("hits", "all"), num.cores = 1L, tasks = 0, 
    verbose = TRUE, ...) {
    
    if (verbose) 
        progressbar <- TRUE else progressbar <- FALSE
    if (inherits(GR, "list") && !(inherits(GR, "InfDiv") || inherits(GR, 
        "pDMP"))) 
        GR <- try(as(GR, "GRangesList"))
    if (Sys.info()["sysname"] == "Linux") {
        bpparam <- MulticoreParam(workers = num.cores, tasks = tasks, 
            progressbar = progressbar)
    } else {
        bpparam <- SnowParam(workers = num.cores, type = "SOCK",
                            progressbar = progressbar)
        BiocParallel::register(bpstart(bpparam))
    }
    
    if (is.character(names(GR))) 
        nams <- names(GR) else nams <- NULL
    
    GR <- bplapply(GR, getGRegionsStat, win.size, step.size, grfeatures, 
        stat, absolute, select.strand, column, prob, entropy, maxgap, 
        minoverlap, scaling, logbase, missings, type, ignore.strand, 
        na.rm, naming, output, BPPARAM = bpparam)
    if (!is.null(nams)) 
        names(GR) <- nams
    return(GR)
}

setClass("InfDiv")
setClass("pDMP")

# ----------------------------------------------------------------------------
# #

#' @aliases getGRegionsStat, list-method
#' @importFrom GenomeInfoDb seqnames seqlengths
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom data.table data.table
#' @importFrom BiocParallel MulticoreParam SnowParam bplapply
setMethod("getGRegionsStat", signature(GR = "list"), function(GR, 
    win.size = 350, step.size = 350, grfeatures = NULL, stat = c("sum", 
        "mean", "gmean", "median", "density", "count", "denCount"), 
    absolute = FALSE, select.strand = NULL, column = 1L, prob = FALSE, 
    entropy = FALSE, maxgap = -1L, minoverlap = 0L, scaling = 1000L, 
    logbase = 2, missings = 0, type = c("any", "start", "end", "within", 
        "equal"), ignore.strand = FALSE, na.rm = TRUE, naming = FALSE, 
    output = c("hits", "all"), num.cores = 1L, tasks = 0, verbose = TRUE, 
    ...) getGRegionsStats(GR, win.size, step.size, grfeatures, stat, 
    absolute, select.strand, column, prob, entropy, maxgap, minoverlap, 
    scaling, logbase, missings, type, ignore.strand, na.rm, naming, 
    output, num.cores, tasks, verbose, ...))


#' @aliases getGRegionsStat, InfDiv-method
#' @importFrom GenomeInfoDb seqnames seqlengths
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom data.table data.table
#' @importFrom BiocParallel MulticoreParam SnowParam bplapply
setMethod("getGRegionsStat", signature(GR = "InfDiv"), function(GR, 
    win.size = 350, step.size = 350, grfeatures = NULL, stat = c("sum", 
        "mean", "gmean", "median", "density", "count", "denCount"), 
    absolute = FALSE, select.strand = NULL, column = 1L, prob = FALSE, 
    entropy = FALSE, maxgap = -1L, minoverlap = 0L, scaling = 1000L, 
    logbase = 2, missings = 0, type = c("any", "start", "end", "within", 
        "equal"), ignore.strand = FALSE, na.rm = TRUE, naming = FALSE, 
    output = c("hits", "all"), num.cores = 1L, tasks = 0, verbose = TRUE, 
    ...) getGRegionsStats(GR, win.size, step.size, grfeatures, stat, 
    absolute, select.strand, column, prob, entropy, maxgap, minoverlap, 
    scaling, logbase, missings, type, ignore.strand, na.rm, naming, 
    output, num.cores, tasks, verbose, ...))


#' @aliases getGRegionsStat, pDMP-method
#' @importFrom GenomeInfoDb seqnames seqlengths
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom data.table data.table
#' @importFrom BiocParallel MulticoreParam SnowParam bplapply
setMethod("getGRegionsStat", signature(GR = "pDMP"), function(GR, 
    win.size = 350, step.size = 350, grfeatures = NULL, stat = c("sum", 
        "mean", "gmean", "median", "density", "count", "denCount"), 
    absolute = FALSE, select.strand = NULL, column = 1L, prob = FALSE, 
    entropy = FALSE, maxgap = -1L, minoverlap = 0L, scaling = 1000L, 
    logbase = 2, missings = 0, type = c("any", "start", "end", "within", 
        "equal"), ignore.strand = FALSE, na.rm = TRUE, naming = FALSE, 
    output = c("hits", "all"), num.cores = 1L, tasks = 0, verbose = TRUE, 
    ...) getGRegionsStats(GR, win.size, step.size, grfeatures, stat, 
    absolute, select.strand, column, prob, entropy, maxgap, minoverlap, 
    scaling, logbase, missings, type, ignore.strand, na.rm, naming, 
    output, num.cores, tasks, verbose, ...))


#' @aliases getGRegionsStat, GRangesList-method
#' @importFrom GenomeInfoDb seqnames seqlengths
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom data.table data.table
#' @importFrom BiocParallel MulticoreParam SnowParam bplapply
setMethod("getGRegionsStat", signature(GR = "GRangesList"), function(GR, 
    win.size = 350, step.size = 350, grfeatures = NULL, stat = c("sum", 
        "mean", "gmean", "median", "density", "count", "denCount"), 
    absolute = FALSE, select.strand = NULL, column = 1L, prob = FALSE, 
    entropy = FALSE, maxgap = -1L, minoverlap = 0L, scaling = 1000L, 
    logbase = 2, missings = 0, type = c("any", "start", "end", "within", 
        "equal"), ignore.strand = FALSE, na.rm = TRUE, naming = FALSE, 
    output = c("hits", "all"), num.cores = 1L, tasks = 0, verbose = TRUE, 
    ...) getGRegionsStats(GR, win.size, step.size, grfeatures, stat, 
    absolute, select.strand, column, prob, entropy, maxgap, minoverlap, 
    scaling, logbase, missings, type, ignore.strand, na.rm, naming, 
    output, num.cores, tasks, verbose, ...))

