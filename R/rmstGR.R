#' @rdname rmstGR
#'
#' @title Root Mean Square Test for Methylation Analysis
#' @description Count data in MethylIT pipeline is carried in GRanges objects.
#'     This function provides a shortcut to apply the parametric Bootstrap of
#'     2x2 Contingency independence test, which is implemented in function
#'     \code{\link[MethylIT.utils]{bootstrap2x2}}. The  goodness of fit
#'     statistic is the root-mean-square statistic (RMST) or Hellinger
#'     divergence, as proposed by Perkins et al. [1, 2]. Hellinger divergence
#'     (HD) is computed as proposed in [3].
#' @details Samples from each group are pooled according to the statistic
#'     selected (see parameter pooling.stat) and a unique GRanges object is
#'     created with the methylated and unmethylated read counts for each group
#'     (control and treatment) in the metacolumn. So, a contingency table can be
#'     built for range from GRanges object.
#' @param LR A list of GRanges, a GRangesList, a CompressedGRangesList object,
#'     or an object from Methyl-IT downstream analyses: 'InfDiv' or 'pDMP'
#'     object. Each GRanges object from the list must have two columns:
#'     methylated (mC) and unmethylated (uC) counts. The name of each element
#'     from the list must coincide with a control or a treatment name.
#' @param count.col 2d-vector of integers with the indexes of the read count
#'     columns. If not given, then it is asssumed that the methylated and
#'     unmethylated read counts are located in columns 1 and 2 of each GRanges
#'     metacolumns. If object LR is the output of Methyl-IT function
#'     \code{\link[MethylIT]{estimateDivergence}}, then columns 1:4 are the read
#'     count columns: columns 1 and 2 are methylated and unmethylated read
#'     counts from the reference group, while columns 3 and 4 are methylated and
#'     unmethylated read counts from the treatment group, respectively. In this
#'     case, if the requested comparison is reference versus treatment, then no
#'     specification is needed for count.col. The comparison control versus
#'     treatment can be obtained by setting count.col = 3:4 and providing
#'     control.names and treatment.names.
#' @param control.names,treatment.names Names/IDs of the control samples, which
#'     must be included in the variable GR at the metacolumn. Default is NULL.
#'     If NULL, then it is assumed that each GRanges object in LR has four
#'     columns of counts. The first two columns correspond to the methylated and
#'     unmethylated counts from control/reference and the other two columns are
#'     the methylated and unmethylated counts from treatment, respectively.
#' @param stat Statistic to be used in the testing: 'rmst' (root mean square
#'     test) or 'hdiv' (Hellinger divergence test).
#' @param pooling.stat statistic used to estimate the methylation pool: row sum,
#'     row mean or row median of methylated and unmethylated read counts across
#'     individuals. If the number of control samples is greater than 2 and
#'     pooling.stat is not NULL, then they will pooled. The same for treatment.
#'     Otherwise, all the pairwise comparisons will be done.
#' @param tv.cut A cutoff for the total variation distance (TVD; absolute value
#'     of methylation levels differences) estimated at each site/range as the
#'     difference of the group means of methylation levels. If tv.cut is
#'     provided, then sites/ranges k with abs(TV_k) < tv.cut are removed before
#'     to perform the regression analysis. Its value must be NULL or a number
#'     0 < tv.cut < 1.
#' @param hdiv.cut An optional cutoff for the Hellinger divergence (*hdiv*). If
#'     the LR object derives from the previous application of function
#'     \code{\link{estimateDivergence}}, then a column with the *hdiv* values is
#'     provided. If combined with tv.cut, this permits a more effective
#'     filtering of the signal from the noise. Default is NULL.
#' @param hdiv.col Optional. Columns where *hdiv* values are located in each
#'     GRange object from LR. It must be provided if together with *hdiv.cut*.
#'     Default is NULL.
#' @param num.permut Number of permutations.
#' @param pAdjustMethod method used to adjust the results; default: BH
#' @param pvalCutOff cutoff used when a p-value adjustment is performed
#' @param saveAll if TRUE all the temporal results are returned
#' @param num.cores,tasks Paramaters for parallele computation using package
#'     \code{\link[BiocParallel]{BiocParallel-package}}: the number of cores to
#'     use, i.e. at most how many child processes will be run simultaneously
#'     (see \code{\link[BiocParallel]{bplapply}} and the number of tasks per job
#'     (only for Linux OS).
#' @param verbose if TRUE, prints the function log to stdout
#' @param progressbar logical(1). Enable progress bar
#' @param ... Additional parameters for function
#'     \code{\link[uniqueGRanges]{MethylIT}}.
#'
#' @importFrom BiocParallel MulticoreParam bplapply SnowParam
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom MethylIT validateClass
#' @return A GRanges object with the original sample counts, bootstrap p-value
#'     probability, total variation (difference of methylation levels), and
#'     p-value adjusment.
#'
#' @examples
#' #' A list of GRanges
#' set.seed(123)
#' sites = 15
#' data <- list(
#'   C1 = data.frame(chr = 'chr1', start = 1:sites,
#'                   end = 1:sites,strand = '*',
#'                   mC = rnbinom(size = 8, mu = 3, n = sites),
#'                   uC = rnbinom(size = 50, mu = 10, n = sites)),
#'   C2 = data.frame(chr = 'chr1', start = 1:sites,
#'                   end = 1:sites, strand = '*',
#'                   mC = rnbinom(size = 8, mu = 3, n = sites),
#'                   uC = rnbinom(size = 50, mu = 10, n = sites)),
#'   T1 = data.frame(chr = 'chr1', start = 1:sites,
#'                   end = 1:sites,strand = '*',
#'                   mC = rnbinom(size = 50, mu = 10, n = sites),
#'                   uC = rnbinom(size = 10, mu = 10, n = sites)),
#'   T2 = data.frame(chr = 'chr1', start = 1:sites,
#'                   end = 1:sites, strand = '*',
#'                   mC = rnbinom(size = 50, mu = 10, n = sites),
#'                   uC = rnbinom(size = 5, mu = 10, n = sites)))
#' #' Transforming the list of data frames into a list of GRanges objects
#' data = lapply(data,
#'               function(x)
#'                 makeGRangesFromDataFrame(x, keep.extra.columns = TRUE))
#'
#' rmstGR(LR = data, control.names = c('C1', 'C2'),
#'        treatment.names = c('T1', 'T2'), pooling.stat = 'sum',
#'        tv.cut = 0.25, num.permut = 100, pAdjustMethod='BH',
#'        pvalCutOff = 0.05, num.cores = 4L, verbose=TRUE)
#' @references
#'     \enumerate{
#'         \item Perkins W, Tygert M, Ward R. Chi-square and Classical Exact
#'         Tests Often Wildly Misreport Significance; the Remedy Lies in
#'         Computers.
#'         [Internet]. Uploaded to ArXiv. 2011. Report No.: arXiv:1108.4126v2.
#'         \item Perkins, W., Tygert, M. & Ward, R. Computing the confidence
#'         levels for a root-mean square test of goodness-of-fit. 217, 9072-9084
#'         (2011).
#'         \item Basu, A., Mandal, A. & Pardo, L. Hypothesis testing for two
#'         discrete populations based on the Hellinger distance. Stat. Probab.
#'         Lett. 80, 206-214 (2010).
#'     }
#'
#' @seealso \code{\link[MethylIT]{FisherTest}}
#' @export
rmstGR <- function(LR, count.col = 1:2, control.names = NULL, treatment.names = NULL, 
    stat = "rmst", pooling.stat = NULL, tv.cut = NULL, hdiv.cut = NULL, 
    hdiv.col = NULL, num.permut = 100, pAdjustMethod = "BH", pvalCutOff = 0.05, 
    saveAll = FALSE, num.cores = 1L, tasks = 0L, verbose = TRUE, progressbar = TRUE, 
    ...) {
    if (any(!unlist(lapply(LR, function(GR) inherits(GR, "GRanges"))))) 
        stop("At least one element from 'LR' is not a 'GRanges' object")
    
    if (inherits(LR, c("InfDiv", "pDMP"))) {
        validateClass(LR)
    }
    if (!is.null(control.names) && !is.null(treatment.names)) 
        LR <- try(LR[c(control.names, treatment.names)], silent = TRUE)
    if (inherits(LR, "try-error")) 
        stop("List's names does not match control & treatment names")
    
    # === Auxiliar function to perform RMST ===
    rmst <- function(GR, num.cores, tasks) {
        count.matrix <- as.matrix(mcols(GR))
        p1 <- count.matrix[, 1:2]
        p1 <- p1[, 1]/rowSums(p1)
        p2 <- count.matrix[, 3:4]
        p2 <- p2[, 1]/rowSums(p2)
        TV <- p2 - p1
        rm(p1, p2)
        gc()
        
        ind <- FALSE
        idx <- c()
        if (!is.null(tv.cut)) 
            idx.tv <- which(abs(TV) >= tv.cut)
        if (!is.null(hdiv.cut) && !is.null(hdiv.col)) 
            idx.hdiv <- which(mcols(GR[, hdiv.col])[, 1] >= hdiv.cut)
        
        if (!is.null(tv.cut) && (!is.null(hdiv.cut) && !is.null(hdiv.col))) {
            idx <- unique(c(idx.tv, idx.hdiv))
            ind <- !is.na(idx)
        } else {
            if (!is.null(tv.cut)) {
                idx <- idx.tv
                ind <- !is.na(idx)
            }
            if (!is.null(hdiv.cut) && !is.null(hdiv.col)) {
                idx <- idx.hdiv
                ind <- !is.na(idx)
            }
        }
        
        if (sum(ind) > 0) {
            idx <- idx[ind]
            count.matrix <- as.matrix(mcols(GR[idx]))
        } else count.matrix <- as.matrix(mcols(GR))
        
        if (nrow(count.matrix) > 1) {
            count.matrix <- count.matrix[, 1:4]
            sites <- nrow(count.matrix)
            GR$TV <- TV
            GR$pvalue <- rep(1, length(GR))
            GR$adj.pval <- rep(1, length(GR))
            count.matrix <- count.matrix[, 1:4]
            count.matrix <- split(count.matrix, row(count.matrix))
            if (verbose) 
                cat("*** Performing", stat, "test... \n
                   # of sites after filtering: ", 
                  sites, "\n")
            if (Sys.info()["sysname"] == "Linux") {
                bpparam <- MulticoreParam(workers = num.cores, tasks = tasks, 
                  progressbar = progressbar)
            } else bpparam <- SnowParam(workers = num.cores, type = "SOCK", 
                progressbar = progressbar)
            
            pvals <- unname(unlist(bplapply(count.matrix, function(v) {
                m <- matrix(as.integer(v), 2, byrow = TRUE)
                bootstrap2x2(x = m, stat = stat, num.permut = num.permut)
            }, BPPARAM = bpparam)))
        } else {
            count.matrix <- count.matrix[1:4]
            pvals <- bootstrap2x2(x = count.matrix, stat = stat, num.permut = num.permut)
        }
        
        if ((!is.null(tv.cut) | !is.null(hdiv.cut)) && !saveAll && 
            length(idx) > 0) {
            GR <- GR[idx]
            GR$pvalue <- pvals
            GR$adj.pval <- p.adjust(pvals, method = pAdjustMethod)
        } else {
            if ((!is.null(tv.cut) | !is.null(hdiv.cut)) && saveAll && 
                length(idx) > 0) {
                GR$pvalue[idx] <- pvals
                GR$adj.pval[idx] <- p.adjust(pvals, method = pAdjustMethod)
            }
            if (is.null(tv.cut) && is.null(hdiv.cut) && saveAll) {
                GR$pvalue <- pvals
                GR$adj.pval <- p.adjust(pvals, method = pAdjustMethod)
            }
        }
        
        if (!is.null(pvalCutOff) && !saveAll) {
            GR <- GR[GR$adj.pval < pvalCutOff]
        }
        return(GR)
    }
    
    if (is.null(control.names) || is.null(treatment.names)) {
        res <- lapply(LR, function(GR) rmst(GR, num.cores = num.cores, 
            tasks = tasks))
    }
    
    if (!is.null(control.names) && !is.null(treatment.names)) {
        ctrl <- LR[control.names]
        treat <- LR[treatment.names]
        
        if (!is.null(pooling.stat)) {
            ctrl <- poolFromGRlist(LR = ctrl, stat = pooling.stat, 
                num.cores = num.cores, verbose = verbose)
            colnames(mcols(ctrl)) <- c("c1", "t1")  # control counts
            
            treat <- poolFromGRlist(treat, stat = pooling.stat, num.cores = num.cores, 
                verbose = verbose)
            colnames(mcols(treat)) <- c("c2", "t2")  # treatment counts
            
            GR <- uniqueGRanges(list(ctrl, treat), verbose = verbose)
            res <- list(rmst(GR, num.cores = num.cores, tasks = tasks))
        } else {
            treat <- lapply(treat, function(GR) {
                GR <- GR[, count.col]
                colnames(mcols(GR)) <- c("c2", "t2")  # treatment counts
                return(GR)
            })
            ctrl <- lapply(ctrl, function(GR) {
                GR <- GR[, count.col]
                colnames(mcols(GR)) <- c("c1", "t1")  # Control counts
                return(GR)
            })
            res <- list()
            i <- 1
            test.name <- c()
            for (j in 1:length(ctrl)) {
                for (k in 1:length(treat)) {
                  test.name <- c(test.name, paste0(control.names[j], 
                    ".", treatment.names[k]))
                  GR <- uniqueGRanges(list(ctrl[[j]], treat[[k]]), 
                    verbose = verbose, ...)
                  if (verbose) 
                    cat("*** Testing", paste0(control.names[k], " versus ", 
                      treatment.names[j]), "\n")
                  res[[i]] <- rmst(GR, num.cores = num.cores, tasks = tasks)
                  i <- i + 1
                }
                names(res) <- test.name
            }
        }
    }
    if (!is.list(res)) 
        res <- list(groupComparison = res)
    if (inherits(LR, "InfDiv") || inherits(LR, "pDMP")) 
        cl <- class(LR) else cl <- class(res)
    res <- structure(res, class = c(cl, "testDMP"))
    return(res)
}
