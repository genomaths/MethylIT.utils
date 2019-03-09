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
#' @param LR Data used to generate the model. A "pDMP"or "InfDiv" object 
#'     obtained with functions \code{\link[MethylIT]{getPotentialDIMP}} or
#'     \code{\link[MethylIT]{estimateDivergence}}. These are list of GRanges
#'     objects, where each GRanges object from the list must have at least two
#'     columns: a column containing the total variation of methylation level
#'     (TV, difference of methylation levels) and a column containing a
#'     divergence of methylation levels (it could be TV or  Hellinger
#'     divergence).
#' @param pred Type of prediction resquested: *density* ("dens"),*quantiles*
#'     ("quant"), *random number* ("rnum"), *probabilities* ("prob"), or
#'     classification *posterior probability* ("postPrb").
#' @param q numeric vector of quantiles, probabilities or an interger if
#'     pred = "rnum".
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

predict.GammaMixt <- function(gmd, LR, pred="quant", q=0.95, div.col=NULL) {
   if (!inherits(LR, "pDMP") && !inherits(LR, "InfDiv"))
       stop("* LR must an object from class 'pPDM' or 'InfDiv'")
  
   if (is.null(div.col)) stop("Please provide a divergence column")
  
   # === Auxiliary functions ===
   dens <- function(x, pars, lambda) {
       alpha = pars[1,]
       beta = pars[2,]
       d <- cbind(dgamma(x, shape = alpha[1], scale = beta[1]),
               dgamma(x, shape = alpha[2], scale = beta[2]))
       d = t(lambda * t(d))
       colnames(d) <- c("CT", "TT")
       return(DataFrame(div=div, d))
   }
   
   PostPrb <- function(x, pars, lambda) {
       res <- dens(x=div, pars=pars, lambda=lambda)[,2:3]
       res <- as.matrix(res)
       return(DataFrame(div=div, res/rowSums(res)))
   }
  
   P <- function(x, pars, lambda) {
       alpha = pars[1,]
       beta = pars[2,]
       p <- cbind(pgamma(x, shape = alpha[1], scale = beta[1]),
               pgamma(x, shape = alpha[2], scale = beta[2]))
       p = sum(t(lambda * t(p)))
       return((DataFrame(div=div, p)))
   }
  
   quantileGMD <- function(x, p, pars, lambda) {
       interval <- c(min(x), max(x))
       G <- function(x) P(x, pars, lambda) - p
       return(uniroot(G, interval)$root) 
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

   LR = unlist(LR)
   div = LR[, div.col]; div <- div$hdiv
   pars <- gmd$gammaPars
   lambda <- gmd$lambda
  
   # ------------------------------------------------------------------------- #  
   res <- switch(pred,
               dens = dens(x=div, pars=pars, lambda=lambda),
               postPrb = PostPrb(x=div, pars=pars, lambda=lambda),
               prob = P(x=div, pars=pars, lambda=lambda),
               quant = quantileGMD(x=div, p=q, pars=pars, lambda=lambda),
               rnum = rGMD(q, pars=pars, lambda=lambda)
   )
   return(res)
}


