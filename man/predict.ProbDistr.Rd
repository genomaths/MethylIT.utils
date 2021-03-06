% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predictProbDistr.R
\name{predict.ProbDistr}
\alias{predict.ProbDistr}
\alias{predict.ProbDistrList}
\title{Predict function for probability distributions in Methyl-IT}
\usage{
\method{predict}{ProbDistr}(nlm, pred = "quant", q = 0.95, dist.name)

\method{predict}{ProbDistrList}(nlm, pred = "quant", q = 0.95, dist.name, num.cores = 1L, tasks = 0L)
}
\arguments{
\item{nlm}{An object carrying the best nonlinear fit for a distribution model
obtained with function \code{\link{nonlinearFitDist}}.}

\item{pred}{Type of prediction resquested: *density* ('dens'),*quantiles*
('quant'), *random number* ('rnum') or *probabilities* ('prob').}

\item{q}{numeric vector of quantiles, probabilities or an interger if
pred = 'rnum'.}

\item{dist.name}{name of the distribution to fit: Weibull2P (default:
'Weibull2P'), Weibull three-parameters (Weibull3P), gamma with
three-parameter (Gamma3P), gamma with two-parameter (Gamma2P),
generalized gamma with three-parameter ('GGamma3P') or four-parameter
('GGamma4P').}

\item{num.cores, tasks}{Paramaters for parallele computation using package
\code{\link[BiocParallel]{BiocParallel-package}}: the number of cores to
use, i.e. at most how many child processes will be run simultaneously
(see \code{\link[BiocParallel]{bplapply}} and the number of tasks per job
(only for Linux OS).}
}
\description{
This is an utility function to get predictions from the
    probability distributions models used in Methyl-IT: Weibull, Gamma, and
    generalized Gamma. Some times, after the nonlinear fit of any of the
    mentioned modelsm we would like to evaluate the model output.
}
\details{
Predictions are based on the best model fit returned by function
    \code{\link{nonlinearFitDist}}. The possible prediction are: *density*,
    *quantiles*, *random number* or *probabilities*.
}
\examples{
set.seed(1)
num.points <- 1000
HD <- makeGRangesFromDataFrame(
  data.frame(chr = 'chr1', start = 1:num.points, end = 1:num.points,
            strand = '*',
            hdiv = rweibull(1:num.points, shape = 0.75, scale = 1)),
  keep.extra.columns = TRUE)
nlms <- nonlinearFitDist(list(HD), column = 1, verbose = FALSE)

x=seq(0.1, 10, 0.05)
y <- predict(nlms[[1]], pred='dens', q = x,
                dist.name='Weibull2P')
y1 <- dweibull(x, shape = 0.75, scale = 1)
# The maximum difference between the 'theoretical' and estimated densities
max(abs(round(y, 2) - round(y1, 2)))

}
