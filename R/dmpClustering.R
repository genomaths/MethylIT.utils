#' @rdname dmpClustering
#' @name dmpClustering
#' @title DMP Clustering
#' @description Given a 'pDMP' object carrying DMPs obtained in Methyl-IT
#'  downstream analysis, function \strong{\emph{'dmpClustering'}} build
#'  clusters of DMPs, which can be further tested to identify differentially
#'  methylated regions (DMRs) with \code{\link{countTest2}} function.
#' @details The number of DMPs reported in each cluster corresponds to the 
#'  numbers of sites inside the cluster where DMPs were found in at least one
#'  of the samples (from control or from treatment). That is,
#'  \strong{dmpClustering} is just a tool to locate regions with high density of
#'  DMPs from all the samples. It does not detect DMRs. It is assumed that only
#'  DMP coordinates are given in the \emph{'dmps'} object. That is, all the
#'  sites provided are considered in the computation.
#' @param dmps An object from \strong{\emph{'pDMP'}} class, which is returned by
#'  \code{\link{selectDIMP}} function or simply a GRanges object carrying DMP
#'  coordinates.
#' @param win.size An integer. The size of the windows/intervals genomics.
#'  Default: \eqn{win.size = 1}.
#' @param step.size Interval at which the regions/windows must be defined. 
#'  Default: \eqn{step.size = 1}.
#' @param minNumDMPs Minimum number of DMPs inside of each cluster.
#'  Default: \eqn{minNumDMPs = 1}.
#' @param maxClustDist Clusters separated by a distance lesser than
#'   'maxClustDist' positions are merged. Default: \eqn{maxClustDist = NULL}.
#' @param method Two different approaches are implemented to clustering DMPs:
#' \describe{
#'   \item{\strong{"relaxed":}}{DMP ranges which are separated by a distance  
#'          less than \emph{'maxClustDist'} are merged and ranges with less 
#'          than \emph{'minNumDMPs'} are removed.}
#'   \item{\strong{"fixed.int":}}{A partition of the ranges covered by the DMPs
#'          is built at fixed intervals \emph{'win.size'} and at fixed step 
#'          \emph{'step.size'}. next, ranges which are separated by a distance  
#'          less than \emph{'maxClustDist'} are merged and ranges with less than
#'          \emph{'minNumDMPs'} are removed.}
#' } 
#' 
#' @param ignore.strand Same as in 
#'  \code{\link[GenomicRanges]{findOverlaps-methods}.}
#' @param verbose if TRUE, prints the function log to stdout 	
#' @return A GRanges object carrying the coordinates of DMP clusters from all 
#'  the samples and the number of DMPs on each of them. 
#' @examples
#' ## Creates a GRanges object carrying DMPs. Notice that only the DMP 
#' ## coordinates are needed.
#' gr <- GRanges(seqnames = Rle( c('chr1', 'chr2', 'chr3', 'chr4'),
#'             c(5, 5, 5, 5)),
#'             ranges = IRanges(start = 1:20, end = 1:20),
#'             strand = rep(c('+', '-'), 10))
#' 
#' ## Simple DMP clustering ignoring the DNA strand
#' dmpClustering(gr, win.size = 4,  step.size = 4, minNumDMPs = 2,
#'               method = "fixed.int")
#' 
#' ## Now, the information on the DNA strand is included in the clustering
#' dmpClustering(dmps = gr, win.size = 4,  step.size = 4, minNumDMPs = 2,
#'               method = "fixed.int", ignore.strand = FALSE)
#' 
#' ## Next, as before adding that clusters separated by a distance lesser than
#' ## 'maxClustDist = 2' will be merged
#' dmpClustering(dmps = gr, win.size = 4,  step.size = 4, minNumDMPs = 2,
#'               method = "fixed.int", maxClustDist = 2, ignore.strand = FALSE)
#' 
#' ## Finally, the relaxed approach. Notice that only two parameter values are 
#' ## needed
#' dmpClustering(gr, minNumDMPs = 2, maxClustDist = 2,
#'               method = "relaxed", ignore.strand = FALSE)
#' @author Robersy Sanchez (\url{https://github.com/genomaths}).
#' @export
dmpClustering <- function(dmps, win.size = 1, 
                        step.size = 1, minNumDMPs = 1,
                        maxClustDist = NULL, 
                        method = c("relaxed","fixed.int"), 
                        ignore.strand = TRUE, verbose = FALSE) {
  
    if (!(inherits(dmps, "pDMP") || inherits(dmps, "GRanges")))  
        stop('\n*** "dmps" objec must inherits from "pDMP" or "GRanges-class"')
  
    method <- match.arg(method)
    
    if (inherits(dmps, "pDMP"))   
        dmps <- uniqueGRanges(dmps, columns = 9L, type = "equal", 
                            ignore.strand = ignore.strand, verbose = FALSE)
    
    
    if (ncol(mcols(dmps)) > 0) mcols(dmps) <- NULL
    
    if (method == "fixed.int") {
        ## Progress bar
        if (verbose) {
            # setup progress bar
            pb <- txtProgressBar(max = 100, style = 3)
            on.exit(close(pb))  # on exit, close progress bar
            setTxtProgressBar(pb, 1)  # update progress bar
        }
        dmrs <- getGRegionsStat2(GR = dmps, win.size = win.size,
                                step.size = step.size, stat = "count",
                                ignore.strand = ignore.strand, 
                                verbose = FALSE)
        colnames(mcols(dmrs)) <- "dmps"
        
        if (verbose) setTxtProgressBar(pb, 25)  # update progress bar
        
        if (!is.null(maxClustDist)) {
            dmrs <- reduce(dmrs, min.gapwidth = maxClustDist,
                            ignore.strand = ignore.strand)
            
            dmrs <- countSignal(signal = dmps, gr = dmrs, 
                                ignore.strand = ignore.strand, 
                                verbose = FALSE)
            
            if (verbose) setTxtProgressBar(pb, 50)  # update progress bar
            
            dmrs <- dmrs[, "sites"]
            colnames(mcols(dmrs)) <- "dmps"
        }
        if (verbose) setTxtProgressBar(pb, 75)  # update progress bar
        
        dmrs <- dmrs[ dmrs$dmps >= minNumDMPs ]

        if (verbose) setTxtProgressBar(pb, 100)  # update progress bar
    } 
    
    if (method == "relaxed") {
        if (is.null(maxClustDist)) 
            stop('\nIf method = "relaxed", then a value for "maxClustDist"',  
                ' must be provided')
        dmrs <- reduce(dmps, min.gapwidth = maxClustDist, 
                        ignore.strand = ignore.strand)
        dmrs <- countSignal(signal = dmps, gr = dmrs, 
                            ignore.strand = ignore.strand, verbose = verbose)
        dmrs <- dmrs[ dmrs $sites >= minNumDMPs, "sites" ]
        colnames(mcols(dmrs)) <- "dmps"
    }
    return(dmrs[, "dmps"])
}
