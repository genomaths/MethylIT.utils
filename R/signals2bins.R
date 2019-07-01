# ================================= signals2bins ============================= #
#' @rdname signals2bins 
#' @title Genomic Signals to Summarized Bins 
#' @description This function summarizes a genomic signal (variable) split into 
#'     bins (intervals). The signal must be provided in the metacolumn of a
#'     \code{\link[GenomicRanges]{GRanges-class}} object. 
#' @details This function is useful, for example, to get the  profile of the 
#'     metylation signal around genes regions: gene-body plus 2kb upstream of
#'     the TSS and 2kb downtream of the TES. The intensity of the signal profile
#'     would vary depending on the sample conditions. If a given treatment has
#'     an effect on methylation then the intesity of the signal profile for the
#'     treatment would go over or below the control samples. 
#'     
#'     This function does the same as function \code{\link{signal2bins}},
#'     except for that it is significantly faster than
#'     \code{\link{signal2bins}} function and small variation on the signal
#'     profiles. These variations came from the way to split the regions into 
#'     bins, for which there is not an exact algorithm to perform it. Function
#'     \code{\link{signal2bins}} uses \code{\link[base]{cut}}, while current
#'     function uses \emph{tile} function
#'     (\code{\link[IRanges]{IPosRanges-class}}).
#' @param signal Preferibly a single GRanges object with genomic signals in
#'     the meta-columns (each colum carrying a signal) or a list of GRanges
#'     objects, each GRanges carrying a signal in the meta-column. For example,
#'     methylation levels, any variable regularly measuring some genomic
#'     magnitude. This GRanges object can be created by using function
#'     \code{\link[MethylIT]{uniqueGRanges}} from \emph{MethylIT} R package.
#' @param regions A GRanges carrying the genomic region where a summarized 
#'     statistic can be computed. For example, annotated gene coordinates.
#' @param stat Statistic used to estimate the summarized value of the variable
#'     of interest in each interval/window. Posible options are: "mean",
#'     geometric mean ("gmean"), "median", "density", "count" and "sum"
#'     (default). Here, we define "density" as the sum of values from the
#'     variable of interest in the given region devided by the length/width of
#'     the region. The option 'count' compute the number/count of positions in
#'     the specified regions with values greater than zero in the selected
#'     'column'.
#' @param nbins,nbinsUP,nbinsDown An integer denoting the number of bins used to 
#'     split the \emph{regions}, upstream the main regions, and downstream the 
#'     main \emph{regions}, respectively.
#' @param streamUp,streamDown  An interger denonting how many base-pairs 
#'     up- and down-stream the provided \emph{regions} must be include in the
#'     computation. Default is NULLL. 
#' @param absolute Optional. Logic (default: FALSE). Whether to use the absolute
#'     values of the variable provided. For example, the difference of
#'     methylation levels could take negative values (TV) and we would be
#'     interested on the sum of abs(TV), which is sum of the total variation
#'     distance.
#' @param na.rm Logical value. If TRUE, the NA values will be removed
#' @param missings Whether to write '0' or 'NA' on regions where there is not
#'     data to compute the statistic.
#' @param region.size An integer. The minimun size of a region to be included in
#'     the computation. Default 300 (bp).  
#' @param verbose Logical. Default is TRUE. If TRUE, then the progress of the
#'     computational tasks is given.
#' @param ... Arguments to pass to 
#'     \code{\link[GenomicRanges]{findOverlaps-methods}} function.
#' @importFrom GenomeInfoDb seqlevels
#' @importFrom GenomicRanges makeGRangesFromDataFrame GRanges findOverlaps
#' @importFrom S4Vectors subjectHits queryHits mcols mcols<- 
#' @importFrom stats median
#' @importFrom dplyr group_by summarise_all '%>%'
#' @importFrom utils txtProgressBar
#' @importFrom MethylIT uniqueGRanges
#' @return A data.frame object carrying the bin coordinates: \emph{binCoord} 
#'     and, for each sample, the signal summarized in the requested statistic:
#'     \emph{statSumary}. Notice that the bin coordinates are relative to
#'     original coordinates given in the \emph{GR} objeect. For example, if the
#'     \emph{GR} object carries genome-wide metylation signals (from several
#'     samples) and we are interested in to get the methylation signal profile
#'     around the genes regions, then we must provide the gene annotated
#'     coordinates in the argument \emph{regions}, and set up the amount of bp
#'     upstream of TSS and dowstream of TES, say, \emph{streamUp} = 2000 and
#'     \emph{streamDown} = 2000, repectively. Next, if we set nbins = 20L,
#'     nbinsUP = 20L, nbinsDown = 20L, then the first and the last 20 bins of
#'     the returned signal profile represent 2000 bp each of them. Since
#'     gene-body sizes vary genome-wide, there is not a specific number of bp
#'     represented by the 20 bins covering the gene-body regions.
#' @seealso \code{\link{signal2bins}}.
#' @export
#' @author Robersy Sanchez. \url{https://genomaths.com}

signals2bins <- function(signal, regions, stat = "mean", nbins = 20L,
                       nbinsUP = 20L, nbinsDown = 20L, streamUp = NULL, 
                       streamDown = NULL, absolute = FALSE, na.rm = TRUE,
                       missings = 0, region.size = 300, scaling = 1000L, 
                       verbose = TRUE, ...) {
   t1 <- Sys.time()
   if (!inherits(regions, "GRanges")) 
       stop("*** 'regions' argument must be a GRanges object")

   if (inherits(signal, "list")) {
      signal <- uniqueGRanges(signal, ...)
   }
   
   signal.chr <- NULL
   if (inherits(signal, "GRanges")) {
       if (region.size < nbins) 
           stop("* Minimum 'region.size' must be greater than 'nbins'")
       widths <- width(regions)
       regions <- regions[widths > region.size]
       if (length(regions) == 0) 
           stop("* There is not region with width > region.size")
       
       chrs <- as.character(seqnames(regions))
       signal.chr <- unique(as.character(seqnames(signal)))
       idx <- match(signal.chr, unique(chrs))
     
       if (all(is.na(idx))) 
           stop("*** chromosomes in the signal did not", 
                   " match regions chromosomes")
       else {
           idx <- unique(chrs)[na.omit(idx)]
           seqlevels(signal, pruning.mode = "coarse") <- idx
       }
   }
   
   if (absolute) mcols(signal) <- abs(as.matrix(mcols(signal)))
   
   if (verbose) cat("* Computing bins for the main regions, ... \n")
   bdr <- binbuilder(regions = regions, num.bins = nbins, verbose = verbose)
   
   if (verbose) cat("* Computing summarized statistic for main regions, ... \n")
   bdr.stat <- statRegion(signal = signal, regions = bdr, stat = stat,
                           na.rm = na.rm, missings = missings, 
                           scaling = scaling, ...)
   if (verbose) 
       cat("--- System elapsed time", format.difftime(Sys.time() - t1), "\n\n")

   if (!is.null(streamUp)) {
       upr <- GeneUpDownStream(GR = regions, upstream = streamUp, onlyUP = TRUE)
       cat("* Computing bins for upstream regions, ... \n")
       upr <- binbuilder(regions = upr, num.bins = nbinsUP, verbose = verbose)
       if (verbose) 
           cat("* Computing summarized statistic for upstream regions, ... \n")
       upr.stat <- statRegion(signal = signal, regions = upr, stat = stat,
                               na.rm = na.rm, missings = missings,
                               scaling = scaling, ...)
       if (verbose) 
          cat("--- System elapsed time", format.difftime(Sys.time() - t1),
               "\n\n")
   }
   
   if (!is.null(streamDown)) {
       dwr <- GeneUpDownStream(GR = regions, downstream = streamDown, 
                               onlyDown = TRUE)
       cat("* Computing bins for downstream regions, ... \n")
       dwr <- binbuilder(regions = dwr, num.bins = nbinsDown, verbose = verbose)
       if (verbose) 
           cat("* Computing summarized statistic for downstream regions",
               "... \n")
       dwr.stat <- statRegion(signal = signal, regions = dwr, stat = stat,
                               na.rm = na.rm, missings = missings,
                               scaling = scaling, ...)
   }

   if (!is.null(streamDown) && !is.null(streamUp)) {
      statSumary <- rbind(upr.stat, bdr.stat, dwr.stat)
      m <- nbins + nbinsUP + nbinsDown
   }
   
   if (is.null(streamDown) && !is.null(streamUp)) {
      statSumary <- rbind(upr.stat, bdr.stat)
      m <- nbins + nbinsUP
   }
   
   if (!is.null(streamDown) && is.null(streamUp)) {
      statSumary <- rbind(bdr.stat, dwr.stat)
      m <- nbins + nbinsDown
   }
   
   if (is.null(streamDown) && is.null(streamUp)) {
      statSumary <- bdr.stat
      m <- nbins
   }
   
   if (verbose) 
       cat("--- System elapsed time", format.difftime(Sys.time() - t1), "\n\n")
   return(data.frame(binCoord = seq(1, m, 1), statSumary[, -1]))
}

# ====================== Auxiliary function to build bins ======================

binbuilder <- function(regions, num.bins, verbose) {
   bdr <- tile(regions, n = num.bins)
   bdr <- unlist(bdr)
   bdr$bins <- rep(1:num.bins,length(regions))
   return(bdr)
}

# ====== Auxiliary function to compute summarized statistic for regions ========

statRegion <- function(signal, regions, stat, na.rm, missings, 
                        maxgap, minoverlap, ignore.strand, type,
                        scaling, ...) {
   
   ## ------------------------- Stistics to use --------------------------------
   stats <- function(x, stat = c(), na.rm) {
      x <- switch(stat,
                  count = sum(x > 0, na.rm = na.rm),
                  sum = sum(x, na.rm = na.rm),
                  mean = mean(x, na.rm = na.rm),
                  gmean = exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x)),
                  median = median(x, na.rm = na.rm),
                  density = sum(x, na.rm = na.rm))
   }
   fn <- function(x) stats(x, stat = stat, na.rm = na.rm)
   # ------------------------------------------------------------------------- #
   
   hits <- findOverlaps(signal, regions, ...)
   if (length(hits) > 0) {
       bins <- regions$bins
       m <- ncol(mcols(signal))
       if (m > 1) {
           mcols(regions) <- matrix(missings, nrow = length(regions), ncol = m)
       }
       else mcols(regions) <- missings
       mcols(regions[subjectHits(hits)]) <- mcols(signal[queryHits(hits)])
       colnames(mcols(regions)) <- colnames(mcols(signal))
       if (stat == "density") {
           widths <- width(regions)
           mcols(regions) <- (scaling * as.matrix(mcols(regions))/widths)
       }

       signal <- regions
       signal$bins <- bins; rm(regions, bins); gc()
       names(signal) <- NULL
       signal <- as.data.frame(signal)
       signal <- signal[, -c(1:5)]

       signal <- signal %>% group_by(bins) %>% summarise_all(list(fn))
      
       return(as.data.frame(signal))
   }
}

