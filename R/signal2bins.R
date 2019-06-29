# ================================= regions2bins ============================= #
#' @rdname signal2bins 
#' @title Genomic Signal to Summarized Bins 
#' @description This function summarizes a genomic signal (variable) split into 
#'     bins (intervals). The signal must be provided in the metacolumn of a
#'     \code{\link[GenomicRanges]{GRanges-class}} object.
#' @details This function is useful, for example, to get the  profile of the 
#'     metylation signal around genes regions: gene-body plus 2kb upstream of
#'     the TSS and 2kb downtream of the TES. The intensity of the signal profile
#'     would vary depending on the sample conditions. If a given treatment has
#'     an effect on methylation then the intesity of the signal profile for the
#'     treatment would go over or below the control samples.
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
#' @param num.cores,tasks Paramaters for parallele computation using package
#'     \code{\link[BiocParallel]{BiocParallel-package}}: the number of cores to
#'     use, i.e. at most how many child processes will be run simultaneously
#'     (see \code{\link[BiocParallel]{bplapply}} and the number of tasks per job
#'     (only for Linux OS).
#' @param missings Whether to write '0' or 'NA' on regions where there is not
#'     data to compute the statistic.
#' @param verbose Logical. Default is TRUE. If TRUE, then the progress of the
#'     computational tasks is given.
#' @param ... Argumetns to pass to \code{\link[MethylIT]{uniqueGRanges}} 
#'     function if \emph{GR} is a list of GRanges objects.
#' @importFrom BiocParallel MulticoreParam SnowParam
#' @importFrom GenomeInfoDb seqlevels
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom S4Vectors mcols
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
#' @seealso A faster version: \code{\link{signals2bins}}.
#' @export
#' @author Robersy Sanchez. \url{https://genomaths.com}

signal2bins <- function(signal, regions, stat = "mean", nbins = 20L,
                        nbinsUP = 20L, nbinsDown = 20L, streamUp = NULL, 
                        streamDown = NULL, absolute = FALSE, na.rm = TRUE,
                        missings = 0, num.cores = 4L, tasks = 0L,
                        verbose = TRUE, ...) {
  t1 <- Sys.time()
  if (!inherits(regions, "GRanges")) 
    stop("*** 'regions' argument must be a GRanges object")
  
  if (inherits(signal, "list")) {
    signal <- uniqueGRanges(signal, ...)
  }
  
  signal.chr <- NULL
  if (inherits(signal, "GRanges")) {
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
  
  if (Sys.info()['sysname'] == "Linux") {
    bpparam <- MulticoreParam(workers = num.cores, tasks = 0L)
  } else {
    bpparam <- SnowParam(workers = num.cores, type = "SOCK")
  }
  
  if (absolute) mcols(signal) <- abs(as.matrix(mcols(signal)))
  
  nams <- colnames(mcols(signal))
  
  if (verbose) cat("* Computing bins for the main regions, ... \n")
  bdr <- binBuilder(regions = regions, num.bins = nbins, num.cores = num.cores,
                    verbose = verbose)
  
  if (verbose) cat("* Computing summarized statistic for main regions, ... \n")
  bdr.stat <- statRegions(signal = signal, regions = bdr, stat = stat, 
                          na.rm = na.rm, missings = missings,
                          num.cores = num.cores, verbose = verbose)
  m <- length(nams)
  if (verbose) {
    cat("* Summarizing the main regions, ... \n")
    cat("--- System elapsed time", format.difftime(Sys.time() - t1), "\n\n")
  }
  
  ## Progress bar
  if(verbose) {
    # setup progress bar
    pb <- txtProgressBar(max = m + 1, style = 3) 
  }
  
  if(verbose) setTxtProgressBar(pb, 1) # update progress bar
  statR <- matrix(NA, nrow = m, ncol = nbins)
  for (k in 1:m) {
    st <- plyr::ldply(bdr.stat, function(x) x[, k])
    st <- colMeans(st, na.rm = na.rm)
    st[is.na(st)] <- 0
    statR[k, ] <- st
    if(verbose) setTxtProgressBar(pb, k + 1) # update progress bar
  }
  close(pb)
  
  if (!is.null(streamUp)) {
    upr <- GeneUpDownStream(GR = regions, upstream = streamUp, onlyUP = TRUE)
    cat("* Computing bins for upstream regions, ... \n")
    upr <- binBuilder(regions = upr, num.bins = nbinsUP, 
                      num.cores = num.cores, verbose = verbose)
    
    if (verbose) 
      cat("* Computing summarized statistic for upstream regions, ... \n")
    upr.stat <- statRegions(signal = signal, regions = upr,
                            stat = stat, na.rm = na.rm, missings = missings,
                            num.cores = num.cores, verbose = verbose)
    
    if (verbose) {
      cat("* Summarizing upstream regions, ... \n")
      cat("--- System elapsed time", format.difftime(Sys.time() - t1),
          "\n\n")
    }
    m <- length(nams)
    
    ## Progress bar
    if(verbose) {
      # setup progress bar
      pb <- txtProgressBar(max = m + 1, style = 3) 
    }
    
    if(verbose) setTxtProgressBar(pb, 1) # update progress bar
    statUp <- matrix(NA, nrow = m, ncol = nbinsUP)
    for (k in 1:m) {
      st <- plyr::ldply(upr.stat, function(x) x[, k])
      st <- colMeans(st, na.rm = na.rm)
      st[is.na(st)] <- 0
      statUp[k, ] <- st
      if(verbose) setTxtProgressBar(pb, k + 1) # update progress bar
    }
    close(pb)
  }
  
  if (!is.null(streamDown)) {
    dwr <- GeneUpDownStream(GR = regions, downstream = streamDown, 
                            onlyDown = TRUE)
    cat("* Computing bins for downstream regions, ... \n")
    dwr <- binBuilder(regions = dwr, num.bins = nbinsDown, 
                      num.cores = num.cores, verbose = verbose)
    if (verbose) 
      cat("* Computing summarized statistic for downstream regions",
          "... \n")
    dwr.stat <- statRegions(signal = signal, regions = dwr,
                            stat = stat, na.rm = na.rm, missings = missings,
                            num.cores = num.cores, verbose = verbose)
    
    if (verbose) {
      cat("* Summarizing downstream regions, ... \n")
      cat("--- System elapsed time", format.difftime(Sys.time() - t1),
          "\n\n")
    }
    m <- length(nams)
    
    ## Progress bar
    if(verbose) {
      # setup progress bar
      pb <- txtProgressBar(max = m + 1, style = 3) 
    }
    
    if(verbose) setTxtProgressBar(pb, 1) # update progress bar
    statDw <- matrix(NA, nrow = m, ncol = nbinsDown)
    for (k in 1:m) {
      st <- plyr::ldply(dwr.stat, function(x) x[, k])
      st <- colMeans(st, na.rm = na.rm)
      st[is.na(st)] <- 0
      statDw[k, ] <- st
      if(verbose) setTxtProgressBar(pb, k + 1) # update progress bar
    }
    close(pb)
  }
  
  statSumary <- t(cbind(statUp, statR, statDw))
  colnames(statSumary) <- nams
  m <- nbins + nbinsUP + nbinsDown
  if (verbose) 
    cat("--- System elapsed time", format.difftime(Sys.time() - t1), "\n\n")
  return(data.frame(binCoord = seq(1, m, 1), statSumary))
}

# ====================== Auxiliary function to build bins ======================

binBuilder <- function(regions, num.bins, num.cores, verbose) {
  
  starts <- start(regions)
  ends <- end(regions)
  chrs <- as.character(seqnames(regions))
  strands <- as.character(strand(regions))
  max.end <- max(ends, na.rm = TRUE)
  min.start <- min(starts, na.rm = TRUE)
  
  if (verbose) progressbar = TRUE else progressbar = FALSE
  if (Sys.info()['sysname'] == "Linux") {
    bpparam <- MulticoreParam(workers = num.cores, tasks=0L,
                              progressbar = progressbar)
  } else {
    bpparam <- SnowParam(workers = num.cores, type = "SOCK",
                         progressbar = progressbar)
  }
  
  lgr <- bplapply(1:length(regions), function(k) {
    s <- seq(starts[k], ends[k], 1)
    labs <- levels(cut(s, num.bins))
    gr <- cbind(start = as.numeric( sub("\\((.+),.*", "\\1", labs) ),
                end = as.numeric( sub("[^,]*,([^]]*)\\]", "\\1", labs) ))
    
    gr <- DataFrame(gr, chr = chrs[k], strand = strands[k])
    gr$start <- gr$start + 1
    gr <- makeGRangesFromDataFrame(gr)
    return(gr)
  }, BPPARAM = bpparam) 
  return(lgr)
} 


# ====== Auxiliary function to compute summarized statistic for regions ========

statRegions <- function(signal, regions, stat, na.rm, missings, 
                        num.cores, verbose) {
  
  if (verbose) progressbar = TRUE else progressbar = FALSE
  if (Sys.info()['sysname'] == "Linux") {
    bpparam <- MulticoreParam(workers = num.cores, tasks = 0L,
                              progressbar = progressbar)
  } else {
    bpparam <- SnowParam(workers = num.cores, type = "SOCK",
                         progressbar = progressbar)
  }
  
  l <- length(regions)
  st <- bplapply(1:l, function(k) {
    st <- getGRegionsStat2(GR = signal, grfeatures = regions[[k]],
                           missings = missings, stat = stat, verbose = FALSE)
    st <- as.matrix(mcols(st))
    
    return(st)
  }, BPPARAM = bpparam)
  
  return(st)
}
