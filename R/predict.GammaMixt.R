#' @rdname predict.GammaMixt
#' @aliases predict.GammaMixt
#' @title Predict function for the DMP's Mixtures of Gamma Distributions model
#' @description This is an utility function to get posterior probability 
#'     predictions based on a given DMP's Mixtures of Gamma Distributions (GMD)
#'     model, obtained with functio \code{\link{gammaMixtCut}}.  
#' @details Predictions are based on the best model fit returned by function
#'     \code{\link{nonlinearFitDist}}. The possible prediction are: *density*,
#'     *quantiles*, *random number* or *probabilities*.
#' @param gmd An object carrying the best nonlinear fit for a distribution model
#'     obtained with function \code{\link{nonlinearFitDist}} ('GammaMixt' 
#'     class).
#' @param q numeric vector of quantiles, probabilities or an interger if 
#'     pred = "rnum", or A "pDMP"or "InfDiv" object obtained with functions
#'     \code{\link[MethylIT]{getPotentialDIMP}} or
#'     \code{\link[MethylIT]{estimateDivergence}}. These are list of GRanges
#'     objects, where each GRanges object from the list must have at least two
#'     columns: a column containing the total variation of methylation level
#'     (TV, difference of methylation levels) and a column containing a
#'     divergence of methylation levels (it could be TV or  Hellinger
#'     divergence).
#' @param pred Type of prediction resquested: *density* ("dens"),*quantiles*
#'     ("quant"), *random number* ("rnum"), *probabilities* ("prob"), or
#'     classification *posterior probability* ("postPrb").
#' @param dist.name name of the distribution to fit: Weibull2P (default:
#'     "Weibull2P"), Weibull three-parameters (Weibull3P), gamma with
#'     three-parameter (Gamma3P), gamma with two-parameter (Gamma2P),
#'     generalized gamma with three-parameter ("GGamma3P") or four-parameter
#'     ("GGamma4P").
#' @param num.cores,tasks Paramaters for parallele computation using package
#'     \code{\link[BiocParallel]{BiocParallel-package}}: the number of cores to
#'     use, i.e. at most how many child processes will be run simultaneously
#'     (see \code{\link[BiocParallel]{bplapply}} and the number of tasks per job
#'     (only for Linux OS).
#' @importFrom S4Vectors DataFrame
#' @examples
#' 
predict.GammaMixt <- function(gmd, ...) UseMethod("predict", gmd)

predict.GammaMixt <- function(gmd, pred="quant", q=0.95, div.col=NULL,
                              interval = NULL) {
  
   if (!is.null(q) && (inherits(LR, "pDMP") || inherits(LR, "InfDiv"))) {
       if (is.null(div.col)) stop("Please provide a divergence column")
       q = unlist(q)
       q = q[, div.col]; q <- q$hdiv
   }
  
   # === Auxiliary functions ===
   dens <- function(x, pars, lambda) {
       alpha = pars[1,]
       beta = pars[2,]
       d <- cbind(dgamma(x, shape = alpha[1], scale = beta[1]),
               dgamma(x, shape = alpha[2], scale = beta[2]))
       d = t(lambda * t(d))
       colnames(d) <- c("CT", "TT")
       return(DataFrame(div = x, d))
   }
   
   PostPrb <- function(x, pars, lambda) {
       res <- dens(x = x, pars = pars, lambda = lambda)[,2:3]
       res <- as.matrix(res)
       return(DataFrame(div = x, res/rowSums(res)))
   }
  
   P <- function(x, pars, lambda) {
       alpha = pars[1,]
       beta = pars[2,]
       p <- cbind(pgamma(x, shape = alpha[1], scale = beta[1]),
               pgamma(x, shape = alpha[2], scale = beta[2]))
       p = rowSums(t(lambda * t(p)))
       return((DataFrame(div = x, Prob = p)))
   }
  
   quantileGMD <- function(p, pars, lambda, interval) {
       G <- function(x) P(pars = pars, lambda = lambda)$Prob - p
       root <- try(uniroot(G, interval)$root, silent = TRUE)
       if (inherits(root, "try-error")) {
          txt <- paste0("The quatile cannot be estimated in the interval",
                        "provided. Please provide a suitable interval")
          stop(txt)
       }
       return() 
   }
  
   rGMD <- function(q, pars, lambda) {
       q <- as.integer(q)
       if (q == 0) stop ("q must be an integer greater than zero")
       alpha = pars[1,]
       beta = pars[2,]
       U <- runif(q); l <- min(lambda); idx <- which.min(lambda)
       n1 <- sum(U < l)
       n2 <- sum(U > l)
       rand.samples <- c(rgamma(n=n1, shape = alpha[idx], scale = beta[idx]),
                         rgamma(n=n2, shape = alpha[-idx], scale = beta[-idx]))
       return(rand.samples)
   }

   pars <- gmd$gammaPars
   lambda <- gmd$lambda
  
   if (pred == "quant" && is.null(interval)) {
       interval <- c(min(div, na.rm = TRUE), max(div, na.rm = TRUE))
       rmin <- P(interval[1], pars=pars, lambda=lambda)$Prob - q
       rmax <- P(interval[2], pars=pars, lambda=lambda)$Prob - q
       if (!(rmin < 0 && rmax > 0)) 
           txt <- paste0("The 'interval' parameter cannot be estimated from ",
                         "the data. Please provide  'interval' to compute",
                         "the quantile")
           stop(txt)
   }
   
   # ------------------------------------------------------------------------- #  
   res <- switch(pred,
               dens = dens(x=q, pars=pars, lambda=lambda),
               postPrb = PostPrb(x=q, pars=pars, lambda=lambda),
               prob = P(x=q, pars=pars, lambda=lambda),
               quant = quantileGMD(p=q, pars=pars, lambda=lambda,
                                   interval=interval),
               rnum = rGMD(q, pars=pars, lambda=lambda)
   )
   return(res)
}


