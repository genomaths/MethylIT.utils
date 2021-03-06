% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getGRegionsStat2.R
\name{getGRegionsStat2}
\alias{getGRegionsStat2}
\title{Statistic of Genomic Regions}
\usage{
getGRegionsStat2(
  GR,
  win.size = 1,
  step.size = 1,
  grfeatures = NULL,
  stat = c("sum", "mean", "gmean", "median", "density", "count", "denCount"),
  column = NULL,
  absolute = FALSE,
  select.strand = NULL,
  maxgap = -1L,
  minoverlap = 0L,
  select = "all",
  ignore.strand = TRUE,
  type = c("within", "start", "end", "equal", "any"),
  scaling = 1000L,
  logbase = 2,
  missings = 0,
  naming = FALSE,
  na.rm = TRUE,
  verbose = TRUE,
  ...
)
}
\arguments{
\item{GR}{A GRange object carying the variables of interest in the 
GRanges metacolumn.}

\item{win.size}{An integer for the size of the windows/regions size of the
intervals of genomics regions.}

\item{step.size}{Interval at which the regions/windows must be defined}

\item{grfeatures}{A GRanges object corresponding to an annotated genomic
feature. For example, gene region, transposable elements, exons,
intergenic region, etc. If provided, then parameters 'win.size' and
step.size are ignored and the statistics are estimated for 'grfeatures'.}

\item{stat}{Statistic used to estimate the summarized value of the variable
    of interest in each interval/window. Posible options are:
    
\describe{
  \item{\strong{'mean':}}{The mean of values inside each region.}
  \item{\strong{'gmean':}}{The geometric mean of values inside each region.}
  \item{\strong{'median':}}{The median of values inside each region.}
  \item{\strong{'density':}}{The density of values inside each region. That
         is, the sum of values found in each region divided by the width of 
         the region.}
  \item{\strong{'count':}}{Compute the number/count of positions with values
         greater than zero inside each regions.}
  \item{\strong{'denCount':}}{The number of sites with value > 0 inside each 
         region divided by the width of the region.}
  \item{\strong{'sum':}}{The sum of values inside each region.}
}

If \strong{GR} have zero metacolum, then it is set \emph{stat = "count"} and
all the sites are included in the computation.}

\item{column}{Integer number denoting the column where the variable of
interest is located in the metacolumn of the GRanges object.}

\item{absolute}{Optional. Logic (default: FALSE). Whether to use the absolute
values of the variable provided. For example, the difference of
methylation levels could take negative values (TV) and we would be
interested on the sum of abs(TV), which is sum of the total variation
distance.}

\item{select.strand}{Optional. If provided,'+' or '-', then the summarized
statistic is computed only for the specified DNA chain.}

\item{maxgap, minoverlap, type}{See 
\code{\link[IRanges]{findOverlaps-methods}} in the
\strong{IRanges} package for a description of these arguments.}

\item{ignore.strand}{When set to TRUE, the strand information is ignored in
the overlap calculations.}

\item{scaling}{integer (default 1). Scaling factor to be used when
stat = 'density'. For example, if scaling = 1000, then density * scaling
denotes the sum of values in 1000 bp.}

\item{logbase}{A positive number: the base with respect to which logarithms
are computed when parameter 'entropy = TRUE' (default: logbase = 2).}

\item{missings}{Whether to write '0' or 'NA' on regions where there is not
data to compute the statistic.}

\item{naming}{Logical value. If TRUE, the rows GRanges object will be 
given the names(grfeatures). Default is FALSE.}

\item{na.rm}{Logical value. If TRUE, the NA values will be removed.}

\item{verbose}{Logical. Default is TRUE. If TRUE, then the progress of the
computational tasks is given.}
}
\value{
A GRanges object with the new genomic regions and their corresponding
    summarized statistic.
}
\description{
A function to estimate the summarized measures of a specified
    variable given in a GRanges object (a column from the metacolums of the
    GRanges object) after split the GRanges object into intervals.
}
\details{
This function split a Grange object into intervals genomic regions
    (GRs) of fixed size A summarized statistic (mean, median, geometric mean
    or sum) is calculated for the specified variable values from each region.
    Notice that if win.size == step.size, then non-overlapping windows are
    obtained.
}
\examples{
library(GenomicRanges)
set.seed(1)
gr <- GRanges(seqnames = Rle( c('chr1', 'chr2', 'chr3', 'chr4'),
            c(5, 5, 5, 5)),
            ranges = IRanges(start = 1:20, end = 1:20),
            strand = rep(c('+', '-'), 10),
            A = seq(1, 0, length = 20))
gr$B <- runif(20)
grs <- getGRegionsStat2(gr, win.size = 4, step.size = 4)
grs

## Selecting the positive strand
grs <- getGRegionsStat2(gr, win.size = 4, step.size = 4, select.strand = '+')
grs

## Selecting the negative strand
grs <- getGRegionsStat2(gr, win.size = 4, step.size = 4, select.strand = '-')
grs

}
\seealso{
\code{\link{getGRegionsStat}}
}
\author{
Robersy Sanchez (\url{https://github.com/genomaths}).
}
