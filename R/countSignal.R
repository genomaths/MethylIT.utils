#' @rdname countSignal
#' @name countSignal
#' @title Count the Sites Carrying an Arbitrary Signal Inside Genomic Regions
#' @description A simple function to count the number of sites carrying an 
#' arbitrary signal, which are inside given genomic regions.
#' @details Given a GRanges object 'signal' carrying an arbitrary signal on each 
#' range, this function counts the number of sites inside the the given genomic 
#' regions 'gr'.
#' @param signal A GRange object carrying the sites of interest.
#' @param gr A GRange object carrying the regions of interest.
#' @param ignore.strand When set to TRUE, the strand information is ignored in
#'     the overlap calculations.
#' @param verbose Logical. Default is TRUE. If TRUE, then the progress of the
#'     computational tasks is given.
#' @return A GRanges object with the genomic regions that carry signals and the
#' number of sites on them.
#' @importFrom GenomeInfoDb seqnames  
#' @importFrom GenomicRanges GRanges findOverlaps 
#' @importFrom GenomicRanges reduce makeGRangesFromDataFrame
#' @importFrom IRanges IRanges width
#' @importFrom S4Vectors subjectHits queryHits  
#' @importFrom BiocGenerics strand start end
#' @importFrom data.table data.table
#' @export
#' @seealso \code{\link{getGRegionsStat}}
#' @author Robersy Sanchez. \url{https://genomaths.com}
#' @examples 
#' some.signal <- makeGRangesFromDataFrame(data.frame(chr = "chr1", 
#'                                                    start = 1:15,
#'                                                    end = 1:15,
#'                                                    strand = '*'))
#' 
#' some.regions <- makeGRangesFromDataFrame(data.frame(chr = "chr1", 
#'                                                     start = c(2, 8),
#'                                                     end = c(7, 14),
#'                                                     strand = '*'))
#' 
#' countSignal(signal = some.signal, gr = some.regions)
#' 
countSignal <- function(signal, gr, maxDist = NULL, ignore.strand = TRUE,
                        verbose = FALSE) {
    if (!is.null(maxDist)) gr <- reduce(gr, min.gapwidth = maxDist)
    gr$region <- paste(seqnames(gr), start(gr), end(gr), sep = "_")
    signal$region <- NA
    
    # To assign the signal ids to the signal
    Hits <- findOverlaps(signal, gr, ignore.strand = ignore.strand,
                        type = "within")
    signal$region[queryHits(Hits)] <- gr$region[subjectHits(Hits)]
    signal <- unique(signal)
    
    # To count the number of sites inside each region
    if (verbose)  message("*** Counting sites in clusters ...")
    signal <- data.table(as.data.frame(signal))
    signal <- signal[!is.na(signal$region), list(seqnames = unique(seqnames),
                                    start = min(start), end = max(end), 
                                    sites = length(start)),
                by = region]
    signal <- makeGRangesFromDataFrame(data.frame(signal), 
                                        keep.extra.columns = TRUE)
    return(signal)
}
