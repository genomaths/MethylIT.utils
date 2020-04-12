#' @rdname shannonEntr
#' @title Compute Shannon Entropy
#' @description Compute Shannon Entropy of probability vector p.
#' @details By definition, if p_i = 0 for some i, then the entropy is 0. In 
#'     addition, since Shannon entropy is only defined for \eqn{0 <= p <= 1},
#'     the function will return 0 for any values of p out of the definition
#'     range of values.
#' @param p A probability vector, sum(p) = 1.
#' @param logbase A positive number: the base with respect to which logarithms
# are computed (default: logbase = 2).
#' @examples
#' counts = sample.int(10)
#' prob = counts/sum(counts)
#' shannonEntr(prob)
#' 
#' ## Out of the definition range of values for p, 'shannonEntr' will return 0
#' shannonEntr(c(0.5, 1.2, -0.3))
#' @export
#' @author Robersy Sanchez (\url{https://github.com/genomaths}).
shannonEntr <- function(p, logbase = 2) {
    logb <- function(p) {
        n <- length(p)
        logP <- integer(n)
        if (n > 1) {
            idx <- (p > 0 & p < 1)
            logP[idx] <- -log(p[idx], base = logbase)
        } else {
            if (p > 0 & p < 1) 
                logP <- -log(p, base = logbase) else logP <- 0
        }
        return(logP)
    }
    return(p * logb(p) + (1 - p) * logb(1 - p))
}


