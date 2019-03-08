#' @rdname gammaMixtCut
#' @title Cutpoint estimation based on Mixtures of Gamma Distributions
#' @description This functions estimates cutpoint value to classify DMPs into
#'     two classes: 1) from treatment and 2) from control, based on Mixtures of
#'     Gamma Distributions.
#' @details After the estimation of potential DMPs, the pool of DMPs from
#'     control and treatment is assumed that follows mixtures of Gamma
#'     distributions corresponding to two populations. A posterior probability
#'     2d-vector is estimated for each DMP. By default DMPs with a posterior
#'     probability to belong to the treatment group greater than *post.cut =
#'     0.5* is classified as *DMP from treatment*. The post.cut can be modified.
#'     For all the cases \eqn{0 < post.cut < 1}. If parameter *find.cut = TRUE*,
#'     then a search for the best cutpoint in a predifined inteval
#'     (*cut.interval*) is performed calling function
#'     \code{\link{evaluateDIMPclass}}.
#' @param LR A "pDMP"or "InfDiv" object obtained with functions 
#'     \code{\link[MethylIT]{getPotentialDIMP}} or
#'     \code{\link[MethylIT]{estimateDivergence}}. These are list of GRanges
#'     objects, where each GRanges object from the list must have at least
#'     two columns: a column containing the total variation of methylation
#'     level (TV, difference of methylation levels) and a column containing a
#'     divergence of methylation levels (it could be TV or  Hellinger
#'     divergence).
#' @param post.cut Posterior probability to dicide whether a DMPs belong to 
#'     treatment group. Default *post.cut* = 0.5.
#' @param tv.col Column number where the total variation is located in the
#'     metadata from each GRanges object.
#' @param div.col Column number for divergence of methylation levels used in the
#'     estimation of the cutpoints. Default: 9L (hdiv column from an InfDiv
#'     object).
#' @param tv.cut A cutoff for the total variation distance to be applied to each
#'     site/range. Only sites/ranges *k* with \eqn{TVD_{k} > tv.cut} are 
#'     are used in the analysis. Its value must be a number 
#'     \eqn{0 < tv.cut < 1}. Default is \eqn{tv.cut = 0.25}.
#' @param control.names,treatment.names Optional. Names/IDs of the control and
#'     treatment samples, which must be include in the variable LR 
#'     (default, NULL). However, these are required if any of the parameters
#'     *find.cut* or *clas.perf* are set TRUE.
#' @param treatment.names Optional. Names/IDs of the treatment samples, which
#'     must be include in the variable LR (default, NULL).
#' @param column a logical vector for column names for the predictor variables
#'     to be used: Hellinger divergence "hdiv", total variation "TV",
#'     probability of potential DIMP "wprob", and the relative cytosine site
#'     position "pos" in respect to the chromosome where it is located. The
#'     relative position is estimated as (x - x.min)/(x.max - x), where x.min
#'     and x.max are the maximum and minimum for the corresponding chromosome,
#'     repectively. If "wprob = TRUE", then Logarithm base-10 of "wprob" will
#'     be used as predictor in place of "wprob" ( see
#'     \code{\link{evaluateDIMPclass}}).
#' @param find.cut Logic. Wether to search for an optimal cutoff value to
#'     classify DMPs based on given specifications.
#' @param classifier Classification model to use. Option "logistic" applies a
#'     logistic regression model; option "lda" applies a Linear Discriminant
#'     Analysis (LDA); "qda" applies a Quadratic Discriminant Analysis (QDA),
#'     "pca.logistic" applies logistic regression model using the Principal
#'     Component (PCs) estimated with Principal Component Analysis (PCA) as
#'     predictor variables. pca.lda" applies LDA using PCs as predictor
#'     variables, and the option "pca.qda" applies a Quadratic Discriminant
#'     Analysis (QDA) using PCs as predictor variables. 'SVM' applies Support
#'     Vector Machines classifier from R package e1071.
#' @param prop Proportion to split the dataset used in the logistic regression:
#'     group versus divergence (at DIMPs) into two subsets, training and
#'     testing.
#' @param clas.perf Logic. Whether to return the classificaiton performance for 
#'     the estimated cutpoint. Default, FALSE.
#' @param cut.interval 0 < *cut.interval* < 0.1. If *find.cut*= TRUE, the 
#'     interval of treatment group posterior probabilities where to search for a
#'     cutpoint. Deafult *cut.interval* = c(0.5, 0.8).
#' @param cut.incr 0 < *cut.incr* < 0.1. If *find.cut*= TRUE, the sucesive 
#'     increamental values runing on the interval *cut.interval*. Deafult, 
#'     *cut.incr* = 0.01.
#' @param num.cores,tasks Paramaters for parallele computation using package
#'     \code{\link[BiocParallel]{BiocParallel-package}}: the number of cores to
#'     use, i.e. at most how many child processes will be run simultaneously
#'     (see \code{\link[BiocParallel]{bplapply}} and the number of tasks per job
#'     (only for Linux OS).
#' @param stat An integer number indicating the statistic to be used in the
#'     testing when *find.cut* = TRUE. The mapping for statistic names are:
#'     0 = "Accuracy", 1 = "Sensitivity", 2 = "Specificity",
#'     3 = "Pos Pred Value", 4 = "Neg Pred Value", 5 = "Precision",
#'     6 = "Recall", 7 = "F1",  8 = "Prevalence", 9 = "Detection Rate",
#'     10 = "Detection Prevalence", 11 = "Balanced Accuracy", 12 = FDR.
#' @param maximize Whether to maximize the performance indicator given in
#'     parameter 'stat'. Default: TRUE.
#' @param ... Additional arameters to pass to functions
#'     \code{\link{evaluateDIMPclass}} and \code{\link[mixtools]{gammamixEM}}.
#' @importFrom mclust Mclust
#' @importFrom mixtools gammamixEM
#' @importFrom stats uniroot

gammaMixtCut <- function(LR, post.cut = 0.5, div.col=NULL, tv.col=NULL, 
                       tv.cut=0.25, find.cut=FALSE, 
                       control.names=NULL, treatment.names=NULL,
                       column=c(hdiv=FALSE, TV=FALSE, wprob=FALSE, pos=FALSE),
                       classifier=c("logistic", "pca.logistic", "lda",
                                    "svm", "qda","pca.lda", "pca.qda"),
                       prop=0.6, clas.perf = FALSE, cut.interval = c(0.5, 0.8),
                       cut.incr = 0.01, stat = 1, maximize = TRUE, num.cores=1L, 
                       tasks=0L, tol = .Machine$double.eps^0.5,
                       maxiter = 1000, ...) {
  if (!inherits(LR, "pDMP") && !inherits(LR, "InfDiv"))
       stop("* LR must an object from class 'pPDM' or 'InfDiv'")
  if ((find.cut || clas.perf) && is.null(div.col))    
       stop("If findcut or clas.perf is TRUE, a div.col value must be provided")
  
  divs = unlist(LR)
  # To remove divs == 0. The methylation signal only is given for divs > 0
  divs = divs[ abs(divs$hdiv) > 0 ]

  # =============  Obtain a prior classification model ===========
  fit <- Mclust(divs$hdiv, G=2, model="V", prior = priorControl())
  alpha <- fit$parameters$mean^2/fit$parameters$variance$sigmasq
  beta <- fit$parameters$variance$sigmasq/fit$parameters$mean
  pro <- fit$parameters$pro

  # ================= Fit Gamma mixture ==========================
  y1 <- gammamixEM(divs$hdiv, lambda = pro, alpha = alpha, beta = beta,
                   verb = FALSE)
  
  # Auxiliar function to find cutpoint/intersection point of the two gamma
  # distributions
  cutFun <- function(post.cut) {
       idx <- which(y1$posterior[, 2] > post.cut)
       TT <- divs[idx]
       CT <- divs[-idx]
       CT <- unlist(CT)
       TT <- unlist(TT)
       return(max(c(max(CT$hdiv),min(TT$hdiv))))
   }
   
   # === Optimal cutpoint from the intersection of gamma distributios ===
   idx <- y1$posterior[,1] > post.cut
   par1 <- y1$gamma.pars[,1]
   par2 <- y1$gamma.pars[,2]
   lower <- max(min(y1$x[idx]), min(y1$x[-idx]))
   upper <- min(max(y1$x[idx]), max(y1$x[-idx]))
   zerofun <- function(x) {
      dgamma(x, shape = par1[1], scale = par1[2]) - 
      dgamma(x, shape = par2[1], scale = par2[2])
   }
   zero <- unlist(uniroot(zerofun, interval = c(lower, upper), 
                  tol = .Machine$double.eps^0.5, maxiter = 1000))
   
   names(zero) <- c("cutpoint", "error", "iter", "init.it", "estim.prec")
  # -------------------------------------------------------------------- #
  
   if (find.cut) {
       cuts <- seq(cut.interval[1], cut.interval[2], cut.incr)
       k = 1; opt <- FALSE
       while (k < length(cuts) && !opt) {
           dmps <- selectDIMP(LR, div.col = div.col,
                               cutpoint = cutFun(cuts[k]),
                               tv.col=tv.col, tv.cut=tv.cut)
             
           conf.mat <- evaluateDIMPclass(dmps, column = column,
                                       control.names = control.names,
                                       treatment.names = treatment.names,
                                       classifier=classifier, prop=prop, 
                                       output = "conf.mat", num.cores=num.cores,
                                       tasks=tasks, verbose = FALSE, ...)
           if (stat == 0) {
               st <- conf.mat$Performance$overall[1]
               if (st == 1) opt <- TRUE
               k <- k + 1
           } 
           if (is.element(stat, 1:11)) {
               st <- conf.mat$Performance$byClass[stat]
               if (st == 1) opt <- TRUE
               k <- k + 1
           } 
           if (stat == 12) {
               st <- conf.mat$FDR
               if (st == 0) opt <- TRUE
               k <- k + 1
           }
       }
       conf.mat <- c(Cutpoint=cutFun(cuts[k]), PostProbCut = cuts[k], conf.mat)
   }
   # -------------------------------------------------------------------- #
   
   if (clas.perf && !find.cut) {
       dmps <- selectDIMP(LR, div.col = div.col, cutpoint = zero[1], 
                       tv.col=tv.col, tv.cut=tv.cut)
       conf.mat <- evaluateDIMPclass(dmps, column = column,
                                   control.names = control.names,
                                   treatment.names = treatment.names,
                                   classifier=classifier, prop=prop, 
                                   output = "conf.mat", num.cores=num.cores,
                                   tasks=tasks, verbose = FALSE, ...)                   
   }
   # -------------------------------------------------------------------- #
   if (clas.perf && !find.cut) {
       res <- list(gammaMixtureCut=zero, conf.mat = conf.mat, gammaMixture = y1)
       cat("\n")
       cat("Cutpoint estimation with Mixtures of Gamma Distributions \n")
       cat("\n")
       cat("Cutpoint =", zero[1], "\n")
       cat("\n")
       cat("The accessible objects in the output list are: \n")
       print(summary(res))
   }
   if (find.cut) {
       res <- list(gammaMixtureCut=zero, conf.mat = conf.mat, gammaMixture = y1)
       STAT <- c("Accuracy", "Sensitivity", "Specificity", "Pos Pred Value",
                 "Neg Pred Value","Precision", "Recall", "F1", "Prevalence", 
                 "Detection Rate", "Detection Prevalence", "Balanced Accuracy",
                 "FDR")
       cat("\n")
       cat("Cutpoint estimation with Mixtures of Gamma Distributions \n")
       cat("\n")
       cat("Cutpoint =", zero[1], "\n")
       cat("\n")
       cat("Cutpoint search performed using model posterior probabilities \n")
       cat("\n")
       cat("Optimized statistic:", STAT[stat + 1], "=", st, "\n")
       cat("Cutpoint =", cutFun(cuts[k]), "\n")
       cat("PostProbCut =", cuts[k], "\n")
       cat("\n")
       cat("Cytosine sites with treatment PostProbCut >=", cuts[k], "have a \n")
       cat("divergence value >=", cutFun(cuts[k]), "\n")
       cat("\n")
       cat("The accessible objects in the output list are: \n")
       print(summary(res))
   } 
   if (!clas.perf && !find.cut) {
       res <- list(gammaMixtureCut=zero, gammaMixture=y1)
       cat("\n")
       cat("Cutpoint estimation with Mixtures of Gamma Distributions \n")
       cat("\n")
       cat("Cutpoint =", zero[1], "\n")
       cat("\n")
       cat("The accessible objects in the output list are: \n")
       print(summary(res))
   }
   invisible(res)
}


