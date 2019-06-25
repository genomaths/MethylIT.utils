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
#' @param GR Preferibly a GRange object or a list of GRanges objects with the 
#'     variable of interest in the GRanges metacolumn.
#' @param win.size An integer for the size of the windows/regions size of the
#'     intervals of genomics regions.
#' @param step.size Interval at which the regions/windows must be defined
#' @param grfeatures A GRanges object corresponding to an annotated genomic
#'     feature. For example, gene region, transposable elements, exons,
#'     intergenic region, etc. If provided, then parameters 'win.size' and
#'     step.size are ignored and the statistics are estimated for 'grfeatures'.
#' @param stat Statistic used to estimate the summarized value of the variable
#'     of interest in each interval/window. Posible options are: "mean",
#'     geometric mean ("gmean"), "median", "density", "count" and "sum"
#'     (default). Here, we define "density" as the sum of values from the
#'     variable of interest in the given region devided by the length of the
#'     region. The option 'count' compute the number/count of positions in the 
#'     specified regions with values greater than zero in the selected 'column'. 
#' @param absolute Optional. Logic (default: FALSE). Whether to use the absolute
#'     values of the variable provided. For example, the difference of
#'     methylation levels could take negative values (TV) and we would be
#'     interested on the sum of abs(TV), which is sum of the total variation
#'     distance.
#' @param select.strand Optional. If provided,"+" or "-", then the summarized
#'     statistic is computed only for the specified DNA chain.
#' @param maxgap,minoverlap,type See ?findOverlaps in the IRanges package for a
#'     description of these arguments.
#' @param ignore.strand When set to TRUE, the strand information is ignored in
#'     the overlap calculations.
#' @param scaling integer (default 1). Scaling factor to be used when
#'     stat = "density". For example, if scaling = 1000, then density * scaling
#'     denotes the sum of values in 1000 bp.
#' @param logbase A positive number: the base with respect to which logarithms
#'     are computed when parameter 'entropy = TRUE' (default: logbase = 2).
#' @param missings Whether to write '0' or 'NA' on regions where there is not
#'     data to compute the statistic.
#' @param naming Logical value. If TRUE, the rows GRanges object will be 
#'     given the names(grfeatures). Default is FALSE.
#' @param na.rm Logical value. If TRUE, the NA values will be removed
#' @param num.cores,tasks Paramaters for parallele computation using package
#'     \code{\link[BiocParallel]{BiocParallel-package}}: the number of cores to
#'     use, i.e. at most how many child processes will be run simultaneously
#'     (see \code{\link[BiocParallel]{bplapply}} and the number of tasks per job
#'     (only for Linux OS).
#' @param verbose Logical. Default is TRUE. If TRUE, then the progress of the
#'     computational tasks is given.
#' @param ... Argumetns to pass to \code{\link[MethylIT]{uniqueGRanges}} 
#'     function if \emph{GR} is a list of GRanges objects.
#' @return A GRanges object with the new genomic regions and their corresponding
#'     summarized statistic.
#' @examples
#' library(GenomicRanges)
#' gr <- GRanges(seqnames = Rle( c("chr1", "chr2", "chr3", "chr4"),
#'             c(5, 5, 5, 5)),
#'             ranges = IRanges(start = 1:20, end = 1:20),
#'             strand = rep(c("+", "-"), 10),
#'             A = seq(1, 0, length = 20))
#' gr$B <- runif(20)
#' grs <- getGRegionsStat2(gr, win.size = 4, step.size = 4)
#' grs
#' 
#' ## Selecting the positive strand
#' grs <- getGRegionsStat2(gr, win.size = 4, step.size = 4, select.strand = "+")
#' grs
#' 
#' ## Selecting the negative strand
#' grs <- getGRegionsStat2(gr, win.size = 4, step.size = 4, select.strand = "-")
#' grs
#' 
#' @importFrom GenomeInfoDb seqnames seqlengths seqlevels
#' @importFrom GenomicRanges GRanges findOverlaps
#' @importFrom IRanges IRanges
#' @importFrom stats median
#' @importFrom dplyr group_by summarise summarize '%>%'
#' @importFrom S4Vectors subjectHits queryHits DataFrame mcols 
#' @importFrom S4Vectors mcols<-
#' @importFrom BiocGenerics strand start end
#' @importFrom MethylIT uniqueGRanges
#' @importFrom utils txtProgressBar
#' @seealso \code{\link{getGRegionsStat}}
#' @export
#' @author Robersy Sanchez

getGRegionsStat2 <- function(GR, win.size=350, step.size=350, grfeatures=NULL,
            stat = c("sum", "mean", "gmean", "median", "density", "count"),
            absolute = FALSE, select.strand = NULL, maxgap =-1L, 
            minoverlap = 0L, scaling = 1000L, logbase = 2, missings = 0,
            type = c("any", "start", "end", "within", "equal"), 
            ignore.strand = FALSE, na.rm = TRUE, naming = FALSE, 
            verbose = TRUE, ...) {
   
   if (inherits(GR, "list")) {
      GR <- uniqueGRanges(GR, ...)
   }
   
   ## These NULL quiet: no visible binding for global variable 'x2'
   if (class( GR ) != "GRanges") stop( "object must be a GRanges object!")
   if (!is.null(grfeatures) && !inherits(grfeatures,"GRanges")) {
       stop("* 'grfeatures', if provided, must be a GRanges object")
   }
   stat <- match.arg(stat, c("sum", "mean", "gmean", "median", "density", 
                               "count"))
       
   if (!is.element(missings, c(0, NA))) missings <- NA
       
   type <- match.arg(type, c("any", "start", "end", "within", "equal"))

   ## === Some functions to use ===
   stats <- function(x, stat = c(), absolute, na.rm) {
           if (absolute) x = abs(x)
           x <- switch(stat,
                       count = sum(x > 0, na.rm = na.rm),
                       sum = sum(x, na.rm = na.rm),
                       mean = mean(x, na.rm = na.rm),
                       gmean = exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x)),
                       median = median(x, na.rm = na.rm),
                       density = sum(x, na.rm = na.rm))
   }
   fn <- function(x) stats(x, stat = stat, absolute = absolute, na.rm = na.rm)
   sortBySeqnameAndStart <- function(gr) {
       seqlevels(gr) <- sort(seqlevels(gr))
       return(gr[order(as.factor(seqnames(gr)), start(gr)), ])
   }
   ## =============================== ##
   if (!is.null(select.strand)) {
       ## possible values "-", "+", NULL
       if (!is.element(select.strand, unique(strand(GR)))) {
           stop("* The GRanges object does not have strand named ", "'",
               select.strand, "'")
       }
       idx <- which(as.character(strand(GR)) == select.strand)
       GR <- GR[idx]
   }
   chrs <- as.character(unique(seqnames(GR)))

   
   ## Progress bar
   if(verbose) {
      # setup progress bar
      pb <- txtProgressBar(max = 100, style = 3) 
      on.exit(close(pb)) # on exit, close progress bar
   }
   
   ## === If genomic features are not specified ===
   if (is.null(grfeatures)) {
       if (length(GR) < win.size || length(GR) < step.size) 
           stop("* 'GR'length is lesser of 'win.size' or 'step.size'")
       if(verbose) setTxtProgressBar(pb, 1) # update progress bar
      
       all.wins <- sapply(1:length(chrs), function(k) {
                   ## get max length of chromosome
                   max.length <- max(IRanges::end(GR[seqnames(GR) == chrs[k],]))
                   ## get sliding windows
                   numTiles <- floor((max.length -
                                   (win.size - step.size)) / step.size) + 1
                   ranges <- IRanges(start=(1 + 0:(numTiles - 1) * step.size),
                                   width=rep(win.size, numTiles))
                   temp.wins <- GRanges(seqnames = rep(chrs[k], numTiles),
                                   ranges = ranges)
                   return(temp.wins)} 
       )
           
       all.wins <- MethylIT::unlist(all.wins)
               
       ## sites of interest inside of the windows
       Hits <- findOverlaps(GR, all.wins, maxgap = maxgap,
                           minoverlap = minoverlap,
                           ignore.strand = ignore.strand, type = type)
       if (length(Hits) > 0) {
           m <- ncol(mcols(GR))
           if (m  > 1) {
               mcols(all.wins) <- matrix(missings, nrow = length(all.wins),
                                           ncol = m)
           } else mcols(all.wins) <- missings
           mcols(all.wins[subjectHits(Hits)]) <- mcols(GR[queryHits(Hits)])
           colnames(mcols(all.wins)) <- colnames(mcols(GR))
           chr <- seqnames(all.wins)

           ## Variable to mark the limits of each GR
           all.wins$cluster.id <- paste(chr, start(all.wins),
                                        end(all.wins), sep = "_")
           GR <- all.wins; rm(all.wins); gc()
           GR <- as.data.frame(GR)
           GR <- GR[, -c(1:5)]
           if(verbose) setTxtProgressBar(pb, 25) # update progress bar
           
           GR <- GR %>% group_by(cluster.id) %>% summarise_all(list(fn))
           
           if(verbose) setTxtProgressBar(pb, 75) # update progress bar
           
           strands <- matrix(unlist(strsplit(GR$cluster.id, "_")), 
                               ncol = 3, byrow = TRUE)
           GR <- data.frame(GR[, -1], chr = strands[,1], 
                           start = as.numeric(strands[, 2]),
                           end = as.numeric(strands[, 3]), strand =  "*")
           GR <- makeGRangesFromDataFrame(GR, keep.extra.columns = TRUE)
       } else {
               m <- ncol(mcols(GR))
               if (m  > 1) {
                   mcols(all.wins) <- matrix(missings, nrow = length(all.wins), 
                                           ncol = m)
               } else mcols(all.wins) <- missings
               colnames(mcols(all.wins)) <- colnames(mcols(GR))
               GR <- all.wins
               warnings("There is not overlap between the 'GR' & the regions")
       }
   } else {
       ## sites of interest inside of the windows
       if(verbose) setTxtProgressBar(pb, 1) # update progress bar
      
       Hits <- findOverlaps(GR, grfeatures, maxgap=maxgap,
                           minoverlap=minoverlap, ignore.strand=ignore.strand,
                           type=type)
       if (length(Hits) > 0) {
           m <- ncol(mcols(GR))
           if (m > 1) {
               mcols(grfeatures) <- matrix(missings, nrow = length(grfeatures), 
                                           ncol = m)
           }
           else mcols(grfeatures) <- missings
           mcols(grfeatures[subjectHits(Hits)]) <- mcols(GR[queryHits(Hits)])
           colnames(mcols(grfeatures)) <- colnames(mcols(GR))
           chr <- seqnames(grfeatures)
               
           if (class(names(grfeatures)) == "character" && naming) 
               grfeatures$cluster.id <- paste(chr, start(grfeatures),
                                               end(grfeatures), 
                                               names(grfeatures), sep = "_")
           else grfeatures$cluster.id  <- paste(chr, start(grfeatures),
                                               end(grfeatures), sep = "_")
               
           GR <- grfeatures; rm(grfeatures); gc()
           GR <- as.data.frame(GR)
           GR <- GR[, -c(1:5)]
           if(verbose) setTxtProgressBar(pb, 25) # update progress bar
           
           GR <- GR %>% group_by(cluster.id) %>% summarise_all(list(fn))
           
           if(verbose) setTxtProgressBar(pb, 75) # update progress bar
           if (naming) strands <- matrix(unlist(strsplit(GR$cluster.id, "_")), 
                                           ncol = 4, byrow = TRUE)
           else strands <- matrix(unlist(strsplit(GR$cluster.id, "_")), 
                                  ncol = 3, byrow = TRUE)
           GR <- data.frame(GR[, -1], chr = strands[,1], 
                           start = as.numeric(strands[, 2]),
                           end = as.numeric(strands[, 3]), strand =  "*")
           GR <- makeGRangesFromDataFrame(GR, keep.extra.columns = TRUE)
           if (naming) names(GR) <- strands[, 4]
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
   if(verbose) setTxtProgressBar(pb, 100) # update progress bar
   return(GR)
}

