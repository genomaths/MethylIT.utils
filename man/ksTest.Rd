% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ksTest.R
\name{ksTest}
\alias{ksTest}
\title{Kolmogorov-Smirnov statistics}
\usage{
ksTest(
  x,
  CDF = "Weibull",
  pars,
  num.sampl = 999,
  sample.size,
  numcores = 1,
  verbose = TRUE,
  ...
)
}
\arguments{
\item{x}{numerical vector to perform the goodness of fit}

\item{CDF}{the name of the cummulative distribution function (CDF)}

\item{pars}{vector of parameters to evaluate the CDF:
4P GG distribution: c(shape=value, scale=value, mu=value, psi=value)
3P GG distribution: c(shape=value, scale=value, psi=value)
3P Weibull distribution: c(shape=value, scale=value, mu=value)
2P Weibull distribution: c(shape=value, scale=value)}

\item{num.sampl}{number of elements to be sampled}

\item{sample.size}{number of permutations. If sample.size < length(x), then
the test becomes a Monte Carlo test}

\item{numcores}{number of cores}

\item{verbose}{If TRUE, prints the function log to stdout}

\item{...}{other parameters}
}
\value{
gamma distribution CDF
}
\description{
Permutation test for Kolmogorov-Smirnov statistics
}
\examples{
num.samples <- 1000
x <- rweibull(num.samples, shape = 1.01, scale = 1.01)
ksTest(x, pars = c(shape = 1, scale = 1))

}
\references{
Alastair Sanderson. Using R to analyse data statistical and
    numerical data analysis with R
    http://www.sr.bham.ac.uk/~ajrs/R/r-analyse_data.html
}
\author{
Robersy Sanchez - 02/29/2016
}
