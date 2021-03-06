#' @rdname getGRegionsStat2
#' @title Statistic of Genomic Regions
#' @description A function to estimate the summarized measures of a specified
#'     variable given in a GRanges object (a column from the metacolums of the
#'     GRanges object) after split the GRanges object into intervals.
#' @details This function split a Grange object into intervals genomic regions
#'     (GRs) of fixed size A summarized statistic (mean, median, geometric mean
#'     or sum) is calculated for the specified variable values from each region.
#'     Notice that if win.size == step.size, then non-overlapping windows are
#'     obtained.
#' @param GR A GRange object carying the variables of interest in the 
#'     GRanges metacolumn.
#' @param win.size An integer for the size of the windows/regions size of the
#'     intervals of genomics regions.
#' @param step.size Interval at which the regions/windows must be defined
#' @param grfeatures A GRanges object corresponding to an annotated genomic
#'     feature. For example, gene region, transposable elements, exons,
#'     intergenic region, etc. If provided, then parameters 'win.size' and
#'     step.size are ignored and the statistics are estimated for 'grfeatures'.
#' @param stat Statistic used to estimate the summarized value of the variable
#'     of interest in each interval/window. Posible options are:
#'     
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
#' 
#' If \strong{GR} have zero metacolum, then it is set \emph{stat = "count"} and
#' all the sites are included in the computation.
#' @param column Integer number denoting the column where the variable of
#'     interest is located in the metacolumn of the GRanges object.
#' @param absolute Optional. Logic (default: FALSE). Whether to use the absolute
#'     values of the variable provided. For example, the difference of
#'     methylation levels could take negative values (TV) and we would be
#'     interested on the sum of abs(TV), which is sum of the total variation
#'     distance.
#' @param select.strand Optional. If provided,'+' or '-', then the summarized
#'     statistic is computed only for the specified DNA chain.
#' @param maxgap,minoverlap,type See 
#'     \code{\link[IRanges]{findOverlaps-methods}} in the
#'     \strong{IRanges} package for a description of these arguments.
#' @param ignore.strand When set to TRUE, the strand information is ignored in
#'     the overlap calculations.
#' @param scaling integer (default 1). Scaling factor to be used when
#'     stat = 'density'. For example, if scaling = 1000, then density * scaling
#'     denotes the sum of values in 1000 bp.
#' @param logbase A positive number: the base with respect to which logarithms
#'     are computed when parameter 'entropy = TRUE' (default: logbase = 2).
#' @param missings Whether to write '0' or 'NA' on regions where there is not
#'     data to compute the statistic.
#' @param naming Logical value. If TRUE, the rows GRanges object will be 
#'     given the names(grfeatures). Default is FALSE.
#' @param na.rm Logical value. If TRUE, the NA values will be removed.
#' @param verbose Logical. Default is TRUE. If TRUE, then the progress of the
#'     computational tasks is given.
#' @return A GRanges object with the new genomic regions and their corresponding
#'     summarized statistic.
#' @examples
#' library(GenomicRanges)
#' set.seed(1)
#' gr <- GRanges(seqnames = Rle( c('chr1', 'chr2', 'chr3', 'chr4'),
#'             c(5, 5, 5, 5)),
#'             ranges = IRanges(start = 1:20, end = 1:20),
#'             strand = rep(c('+', '-'), 10),
#'             A = seq(1, 0, length = 20))
#' gr$B <- runif(20)
#' grs <- getGRegionsStat2(gr, win.size = 4, step.size = 4)
#' grs
#' 
#' ## Selecting the positive strand
#' grs <- getGRegionsStat2(gr, win.size = 4, step.size = 4, select.strand = '+')
#' grs
#' 
#' ## Selecting the negative strand
#' grs <- getGRegionsStat2(gr, win.size = 4, step.size = 4, select.strand = '-')
#' grs
#' 
#' @importFrom GenomeInfoDb seqnames seqlengths seqlevels
#' @importFrom GenomicRanges GRanges findOverlaps
#' @importFrom IRanges IRanges width
#' @importFrom stats median
#' @importFrom dplyr group_by summarise_all '%>%'
#' @importFrom S4Vectors subjectHits queryHits DataFrame mcols 
#' @importFrom S4Vectors mcols<-
#' @importFrom BiocGenerics strand start end
#' @importFrom MethylIT uniqueGRanges sortBySeqnameAndStart unlist
#' @importFrom utils txtProgressBar
#' @seealso \code{\link{getGRegionsStat}}
#' @author Robersy Sanchez (\url{https://github.com/genomaths}).
#' @export
getGRegionsStat2 <- function(GR, win.size = 1, step.size = 1, 
                            grfeatures = NULL, 
                            stat = c("sum", "mean", "gmean", "median",
                                    "density", "count", "denCount"), 
                            column = NULL, absolute = FALSE, 
                            select.strand = NULL, maxgap = -1L, 
                            minoverlap = 0L, select = "all", 
                            ignore.strand = TRUE, 
                            type = c("within", "start", "end", "equal", "any"),
                            scaling = 1000L, logbase = 2, missings = 0,
                            naming = FALSE, na.rm = TRUE, verbose = TRUE, ...) {
    
    ## These NULL quiet: no visible binding for global variable 'x2'
    if (!inherits(GR, "GRanges")) 
        stop("GR object must inherits from a GRanges class")
    if (!is.null(grfeatures) && !inherits(grfeatures, "GRanges")) {
        stop("* 'grfeatures', if provided, must be a GRanges object")
    }
    stat <- match.arg(stat, c("sum", "mean", "gmean", "median", "density", 
        "count", "denCount"))
    
    if (!is.element(missings, c(0, NA)))  missings <- NA
    
    type <- match.arg(type)
    
    if (!is.null(column)) GR <- GR[, column]
    
    ## === Some functions to use ===
    stats <- function(x, stat = c(), absolute, na.rm) {
        if (absolute) 
            x <- abs(x)
        x <- switch(stat, 
                    count = sum(x != 0, na.rm = na.rm), 
                    sum = sum(x, na.rm = na.rm), 
                    mean = mean(x, na.rm = na.rm), 
                    gmean = exp(sum(log(x[x > 0]), na.rm = na.rm)/length(x)),
                    median = median(x, na.rm = na.rm),
                    density = sum(x, na.rm = na.rm), 
                    denCount = sum(x != 0, na.rm = na.rm))
    }
    fn <- function(x) stats(x, stat = stat, absolute = absolute, na.rm = na.rm)
    
    ## =============================== ##
    if (!is.null(select.strand)) {
        ## possible values '-', '+', NULL
        if (!is.element(select.strand, unique(strand(GR)))) {
            stop("* The GRanges object does not have strand named ", 
                "'", select.strand, "'")
        }
        idx <- which(as.character(strand(GR)) == select.strand)
        GR <- GR[idx]
        ignore.strand <- FALSE
    }

    ## Progress bar
    if (verbose) {
        # setup progress bar
        pb <- txtProgressBar(max = 100, style = 3)
        on.exit(close(pb))  # on exit, close progress bar
    }
    
    ## === If genomic features are not specified ===
    if (is.null(grfeatures)) {
        if (length(GR) < win.size || length(GR) < step.size) 
            stop("* 'GR'length is lesser of 'win.size' or 'step.size'")
        if (verbose) 
            setTxtProgressBar(pb, 1)  # update progress bar
        
        all.wins <- slidingGRanges(GR = GR, win.size = win.size,
                                    step.size = step.size)
        
        ## sites of interest inside of the windows
        Hits <- findOverlaps(GR, all.wins, maxgap = maxgap, 
                            minoverlap = minoverlap, 
                            select = select, 
                            ignore.strand = ignore.strand,
                            type = type)
        if (length(Hits) > 0) {
            m <- ncol(mcols(GR))
            all.wins <- all.wins[subjectHits(Hits)]
            GR <- GR[queryHits(Hits)]
            strand(all.wins) <- strand(GR)
            if (m > 1) {
                mcols(all.wins) <- matrix(missings, nrow = length(all.wins),
                                        ncol = m)
                mcols(all.wins) <- mcols(GR)
                cn <- colnames(mcols(GR))
            } else {
                mcols(all.wins) <- 1
                stat <- "count"
                cn <- stat
            }
            chr <- seqnames(all.wins)
            if (!ignore.strand) strands <- strand(all.wins)
            else strands <- rep("*", length(all.wins))
                
            ## Variable to mark the limits of each GR
            all.wins$cluster.id <- paste(chr, start(all.wins), end(all.wins), 
                                        strands, sep = "_")
            GR <- all.wins
            rm(all.wins); gc()
            GR <- as.data.frame(GR)
            GR <- GR[, -c(1:5)]
            if (verbose) 
                setTxtProgressBar(pb, 25)  # update progress bar
            
            GR <- GR %>% group_by(cluster.id) %>% summarise_all(list(fn))
            
            if (verbose) setTxtProgressBar(pb, 75)  # update progress bar
            
            strands <- matrix(unlist(strsplit(GR$cluster.id, "_")),
                            ncol = 4, byrow = TRUE)
            GR <- data.frame(GR[, -1], chr = strands[, 1], 
                            start = as.numeric(strands[, 2]),
                            end = as.numeric(strands[, 3]),
                            strand = as.character(strands[, 4]))
            GR <- makeGRangesFromDataFrame(GR, keep.extra.columns = TRUE)
            if (stat == "density" || stat == "denCount") {
                widths <- width(GR)
                mcols(GR) <- (scaling * as.matrix(mcols(GR))/widths)
            }
            colnames(mcols(GR)) <- cn
        } else {
            m <- ncol(mcols(GR))
            if (m > 1) {
                mcols(all.wins) <- matrix(missings, nrow = length(all.wins), 
                                            ncol = m)
                colnames(mcols(all.wins)) <- colnames(mcols(GR))
            } else mcols(all.wins) <- missings
            GR <- all.wins
            warnings("There is not overlap between the 'GR' & the regions")
        }
    } else {
        ## sites of interest inside of the windows
        if (verbose) 
            setTxtProgressBar(pb, 1)  # update progress bar
        
        Hits <- findOverlaps(GR, grfeatures, maxgap = maxgap, select = select,
                            minoverlap = minoverlap, 
                            ignore.strand = ignore.strand, 
                            type = type)
        if (length(Hits) > 0) {
            m <- ncol(mcols(GR))
            grfeatures <- grfeatures[subjectHits(Hits)]
            GR <- GR[queryHits(Hits)]
            if (m > 1) {
                mcols(grfeatures) <- matrix(missings, 
                                            nrow = length(grfeatures), 
                                            ncol = m)
                mcols(grfeatures) <- mcols(GR)
            } else {
                mcols(grfeatures) <- 1
                stat <- "count"
            }
            chr <- seqnames(grfeatures)
            if (!ignore.strand) strands <- strand(grfeatures)
            else strands <- rep("*", length(grfeatures))
            
            if (class(names(grfeatures)) == "character" && naming) 
                grfeatures$cluster.id <- paste(chr, start(grfeatures), 
                                                end(grfeatures),
                                                strands,
                                                names(grfeatures), sep = "_") 
            else grfeatures$cluster.id <- paste(chr, start(grfeatures), 
                                                end(grfeatures), strands,
                                                sep = "_")
            
            GR <- grfeatures
            rm(grfeatures)
            gc()
            GR <- as.data.frame(GR)
            GR <- GR[, -c(1:5)]
            if (verbose) setTxtProgressBar(pb, 25)  # update progress bar
            
            GR <- GR %>% group_by(cluster.id) %>% summarise_all(list(fn))
            
            if (verbose) 
                setTxtProgressBar(pb, 75)  # update progress bar
            if (naming) 
                strands <- matrix(unlist(strsplit(GR$cluster.id, "_")),
                                    ncol = 5, byrow = TRUE) 
            else strands <- matrix(unlist(strsplit(GR$cluster.id, "_")),
                                    ncol = 4, byrow = TRUE)
            GR <- data.frame(GR[, -1], chr = strands[, 1], 
                            start = as.numeric(strands[, 2]), 
                            end = as.numeric(strands[, 3]),
                            strand = as.character(strands[, 4]))
            GR <- makeGRangesFromDataFrame(GR, keep.extra.columns = TRUE)
            if (naming) names(GR) <- strands[, 5]
            if (stat == "density") {
                widths <- width(GR)
                mcols(GR) <- (scaling * as.matrix(mcols(GR))/widths)
            }
            colnames(mcols(GR)) <- cn
        } else {
            m <- ncol(mcols(GR))
            if (m > 1) {
                mcols(grfeatures) <- matrix(missings, nrow = length(grfeatures), 
                                            ncol = m)
            } else mcols(grfeatures) <- missings
            colnames(mcols(grfeatures)) <- colnames(mcols(GR))
            GR <- grfeatures
            warnings("There is not overlap between the 'GR' & the regions")
        }
    }
    GR <- sortBySeqnameAndStart(GR)
    if (verbose) 
        setTxtProgressBar(pb, 100)  # update progress bar
    return(GR)
}

