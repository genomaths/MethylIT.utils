#' @name predictDIMPclass
#' @rdname predictDIMPclass
#' @title Predict DIMP class
#' @description This function classify each DIMP as a control or a treatment
#'   DIMP
#' @details Predictions only makes sense if the query DIMPs belong to same
#'   methylation context and derive from an experiment accomplished under the
#'   same condition set for the DIMPs used to build the model.
#' @param LR A list of GRanges objects obtained through the through MethylIT
#'   downstream analysis. Basically, this object is a list of GRanges containing
#'   only differentially methylated position (DMPs). The metacolumn of each
#'   GRanges must contain the columna: Hellinger divergence "hdiv", total
#'   variation "TV", the probability of potential DMP "wprob", which naturally
#'   are added in the downstream analysis of MethylIT.
#' @param model A classifier model obtained with the function
#' 'evaluateDIMPclass'.
#' @param conf.matrix Optional. Logic, whether a confusion matrix should be
#'   returned (default, FALSE, see below).
#' @param control.names Optional. Names/IDs of the control samples, which must
#'   be include in thr variable LR (default, NULL).
#' @param treatment.names Optional. Names/IDs of the treatment samples, which
#'   must be include in the variable LR (default, NULL).
#' @return The same LR object with a column named "class" added to a GRanges
#'   object from LR (default). Based on the model prediction each DIMP is
#'   labeled as control "CT" or as treatment "TT". If "conf.matrix" is TRUE and
#'   the arguments control.names and treatment.names are provided, then the
#'   overall confusion matrix is returned
#' @importFrom MethylIT unlist
#' @importFrom caret confusionMatrix
#' @examples 
#' library(MethylIT)
#' 
#' data(cutpoint, PS, package = "MethylIT")
#' 
#' ## DIMPs are selected using the cupoints
#' DMPs <- selectDIMP(PS, div.col = 9L, cutpoint = cutpoint$cutpoint,
#'                    tv.cut = 0.92)
#' 
#' ## Classification of DIMPs into two clases: DIMPS from control and DIMPs from
#' ## treatment samples and evaluation of the classifier performance (for more
#' ## details see ?evaluateDIMPclass).
#' perf <- evaluateDIMPclass(LR = DMPs,
#'                           column = c(hdiv = TRUE, TV = TRUE,
#'                                      wprob = TRUE, pos = TRUE),
#'                           classifier = "lda", n.pc = 4L,
#'                           control.names =  c("C1", "C2", "C3"),
#'                           treatment.names = c("T1", "T2", "T3"),
#'                           center = TRUE, scale = TRUE,
#'                           prop = 0.6)
#' 
#' #' Now predictions of DIMP for control and treament can be obtained
#' pred = predictDIMPclass(LR = DMPs, model = perf$model,
#'                         conf.matrix = TRUE,
#'                         control.names = c("C1", "C2", "C3"),
#'                         treatment.names = c("T1", "T2", "T3"))
#' @export
predictDIMPclass <- function(LR, model, conf.matrix = FALSE,
                             control.names = NULL,
                             treatment.names = NULL) {
  if (conf.matrix && (is.null(control.names) || is.null(treatment.names))) {
       stop(paste0("* if conf.mat = TRUE, then the character vectors for ",
                   "control.names and treatment.names must be provided"))
  }
  
  if (inherits(LR, "GRangesList")) LR <- as(LR, "list")
  if (class(LR)[1] == "list") class(LR) <- "pDMP"
  r1 <- try(validateClass(LR), silent = TRUE)
  if (inherits(r1, "try-error"))
       stop("*** LR is not an object from 'pDMP' class or is not", 
            " coercible to 'pDMP' class" )
  
  if (conf.matrix && !is.null(control.names) && !is.null(treatment.names)) {
    sn = names(LR)
    idx.ct = match(control.names, sn)
    idx.tt = match(treatment.names, sn)
    CT = LR[ idx.ct ]
    CT = unlist(CT)
    TT = LR[ idx.tt ]
    TT = unlist(TT)
    classSet = list(CT = CT, TT = TT)
  } else classSet = LR

  position <- function(gr) {
    chrs = split(gr, seqnames(gr))
    gr = lapply(chrs, function(grc) {
      x = start(grc)
      x.min = min(x)
      x.max = max(x)
      delta =  max(c(x.max - x, 1))
      return((x - x.min)/(delta))
    })
    return(unlist(gr))
  }

  classifier <- function(GR) {
    cn = colnames(mcols(GR))
    vn = c("hdiv", "TV", "wprob")
    vn = cn[is.element(cn, vn)]
    m = mcols(GR[, vn])
    m$pos = position(GR)
    if (is.element("wprob", vn)) {
      vn = sub("wprob", "logP", vn)
      colnames(m) <- c(vn, "pos")
      m$logP = log10(m$logP + 2.2e-308)
    }

    if (inherits(model, "glm")) {
      model = structure(model, class = c("LogisticR", "glm"))
    }
    if (inherits(model, "lda") || inherits(model, "qda")) {
      pred = predict(model, newdata = m )$class
    } else pred = predict(model, newdata = m, type = "class")
    GR$class = pred
    return(GR)
  }
  LR = lapply(classSet, classifier)
  if (!conf.matrix) {
    return(LR)
  } else {
    TRUE_class <- factor(c(rep("CT", length(CT)), rep("TT", length(TT))), 
                       levels = c("CT", "TT"))
    PRED_class <- factor(c(as.character(LR$CT$class), 
                           as.character(LR$TT$class)),
                       levels = c("CT", "TT"))
    
    conf.mat <- confusionMatrix(data=PRED_class, reference=TRUE_class,
                               positive="TT")
    
    return(conf.mat = conf.mat)}
}

