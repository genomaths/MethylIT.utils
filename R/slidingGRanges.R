#' @rdname slidingGRanges
#' @name slidingGRanges
#' @title Generates intervals for a GRanges objects
#' @description \strong{slidingGRanges} generates sliding intervals of a 
#' specified width.
#' @details This function split a GRange object into intervals regions  of
#'     fixed size. If win.size == step.size, then non-overlapping windows are
#'     obtained. \strong{slidingGRanges} function generates sliding windows
#'     within each range of IRL, according to width and step, returning a
#'     \code{\link[IRanges]{IRangesList-class}}. If the sliding windows do not
#'     exactly cover a range in IRL, the last window is partial.
#' @param GR A \code{\link[GenomicRanges]{GRanges-class}} object or a list 
#'  of \code{\link[GenomicRanges]{GRanges-class}} objects or an object that can
#'  be coerced to a list of \code{\link[GenomicRanges]{GRanges-class}} objects.
#' @param win.size An integer. The size of the windows/intervals genomics.
#' @param step.size Interval at which the regions/windows must be defined
#' @param num.cores,tasks Parameters for parallel computation using package
#'     \code{\link[BiocParallel]{BiocParallel-package}}: the number of cores to
#'     use, i.e. at most how many child processes will be run simultaneously
#'     (see \code{\link[BiocParallel]{bplapply}} and the number of tasks per job
#'     (only for Linux OS).
#' @param verbose Logical. Default is TRUE. If TRUE, then the progress of the
#'     computational tasks is given.
#' @param ... Not in use.
#' @return An object from \code{\link[IRanges]{IRangesList-class}}.
#' @export
#' @importFrom IRanges IRanges 
#' @importFrom GenomicRanges GRanges
#' @importFrom BiocParallel MulticoreParam bpmapply SnowParam
#' @seealso \code{\link{slidingWindowSignals}}
#' @author Robersy Sanchez (\url{https://genomaths.com}).

#' @aliases slidingGRanges
setGeneric("slidingGRanges", function(
                                    GR, 
                                    win.size = 1, 
                                    step.size = 1,
                                    num.cores = 1L,                                       
                                    tasks = 0, 
                                    verbose = FALSE, ...)
    standardGeneric("slidingGRanges")
)

#' @aliases slidingGRanges
setMethod("slidingGRanges",  signature(GR = "GRanges"), 
        function(GR, win.size = 1, step.size = 1) {
    if (length(GR) < win.size || length(GR) < step.size)
            stop("* 'GR'length is lesser of 'win.size' or 'step.size'")
          
    gr <- GRanges()
    chrs <- as.character(unique(seqnames(GR)))
    for (k in seq_along(chrs)) {
        ## get max length of chromosome
        max.length <- max(IRanges::end(GR[seqnames(GR) == chrs[k], ]))
        
        ## get start chromosome coordinate
        start0 <- min(IRanges::start(GR[seqnames(GR) == chrs[k], ]))

        ## get sliding windows
        numTiles <- floor((max.length - (win.size - step.size))/step.size) + 1
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
        gr <- suppressWarnings(c(gr, temp.wins))
    }
    return(gr)
})

# ==================== Function to operate on lists ====================== #
slidingGR <- function(
                    GR,
                    win.size = 1,
                    step.size = 1,
                    num.cores = 1L,
                    tasks = 0, 
                    verbose = FALSE) {
  
    if (verbose)  progressbar <- TRUE else progressbar <- FALSE
    
    if (Sys.info()["sysname"] == "Linux") {
        bpparam <- MulticoreParam(workers = num.cores, tasks = tasks, 
                                  progressbar = progressbar)
    } else {
        bpparam <- SnowParam(workers = num.cores, type = "SOCK",
                            progressbar = progressbar)
        BiocParallel::register(bpstart(bpparam))
    }
      
    if (is.character(names(GR))) nams <- names(GR) else nams <- NULL
        
    GR <- bplapply(GR, slidingGRanges, 
                  win.size = win.size,
                  step.size = step.size, 
                  BPPARAM = bpparam)
    if (!is.null(nams)) names(GR) <- nams
    return(GR)
}

#' @aliases slidingGRanges, list-method
#' @importFrom BiocParallel MulticoreParam SnowParam bplapply
setMethod("slidingGRanges",  signature(GR = "pDMP"), 
          function(
            GR,
            win.size = 1,
            step.size = 1,
            num.cores = 1L,
            tasks = 0, 
            verbose = FALSE) 
            slidingGR(GR, win.size, step.size, num.cores, tasks, verbose))


#' @aliases slidingGRanges, list-method
#' @importFrom BiocParallel MulticoreParam SnowParam bplapply
setMethod("slidingGRanges",  signature(GR = "InfDiv"), 
          function(
            GR,
            win.size = 1,
            step.size = 1,
            num.cores = 1L,
            tasks = 0, 
            verbose = FALSE) 
            slidingGR(GR, win.size, step.size, num.cores, tasks, verbose))

#' @aliases slidingGRanges, list-method
#' @importFrom BiocParallel MulticoreParam SnowParam bplapply
setMethod("slidingGRanges",  signature(GR = "list"), 
          function(
                  GR,
                  win.size = 1,
                  step.size = 1,
                  num.cores = 1L,
                  tasks = 0, 
                  verbose = FALSE) 
            slidingGR(GR, win.size, step.size, num.cores, tasks, verbose))


