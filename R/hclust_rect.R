#' @rdname hclust_rect
#' @title Draw Rectangles with Background Colors Around Hierarchical Clusters
#' @description Draws rectangles with background colors around the branches of a
#'     dendrogram highlighting the corresponding clusters. First the dendrogram 
#'     is cut at a certain level, then a rectangle is drawn around selected 
#'     branches.
#' @details This function is exactly function \code{\link[stats]{rect.hclust}} 
#'     with a nice feature added: "to draw the rectangles around hierarchical
#'     clusters with background colors".
#' @param tree The same as in \code{\link[stats]{rect.hclust}}
#' @param k,h The same as in \code{\link[stats]{rect.hclust}}
#' @param which,x The same as in \code{\link[stats]{rect.hclust}}
#' @param border The same as in \code{\link[stats]{rect.hclust}}
#' @param cluster The same as in \code{\link[stats]{rect.hclust}}
#' @param color Background color to use inside the rectangles around
#'     hierarchical clusters. Default is NULL.
#' @param cuts A numeric vector used to manually locate the rectangles around 
#'     hierarchical clusters in the rigth position. This is tricky since each
#'     experimental dataset yield different measurement scale and must be
#'     manually adjusted. Settings are cuts = c(xleft, ybottom, xright, ytop).
#'     Default is NULL. Use it only if need it.
#' @examples 
#' ### Violent crime rates by US state
#' hca <- hclust(dist(USArrests))
#' 
#' # Basic key RGB colors
#' # rgb(red, green, blue, alpha, names = NULL, maxColorValue = 1)
#' clusters.color = c(rgb(0, 0.7, 0, 0.1), rgb(0, 0, 1, 0.1),
#'                  rgb(1, 0.2, 0, 0.1))
#' 
#' plot(hca)
#' hclust_rect(hca, h = 150,  border = clusters.color, col = clusters.color )
#' @export

hclust_rect <- function (tree, k = NULL, which = NULL, x = NULL, h = NULL,
                       border = 2, cluster = NULL, cuts = NULL, 
                       color = NULL, ... ) {
  
   if (length(h) > 1L | length(k) > 1L) 
       stop("'k' and 'h' must be a scalar")
   if (!is.null(h)) {
       if (!is.null(k)) stop("specify exactly one of 'k' and 'h'")
       k <- min(which(rev(tree$height) < h))
       k <- max(k, 2)
   }
   else if (is.null(k)) stop("specify exactly one of 'k' and 'h'")
   if (k < 2 | k > length(tree$height)) 
       stop(gettextf("k must be between 2 and %d", length(tree$height)), 
           domain = NA)
   if (is.null(cluster)) cluster <- cutree(tree, k = k)
   clustab <- table(cluster)[unique(cluster[tree$order])]
   m <- c(0, cumsum(clustab))
   if (!is.null(x)) {
       if (!is.null(which)) stop("specify exactly one of 'which' and 'x'")
       which <- x
       for (n in seq_along(x)) which[n] <- max(which(m < x[n]))
   }
   else if (is.null(which))  which <- 1L:k
   if (any(which > k)) 
       stop(gettextf("all elements of 'which' must be between 1 and %d", k),
           domain = NA)
   border <- rep_len(border, length(which))
   retval <- list()
   l = seq_along(which)
   if (is.null(color) || length(color) < l) color = rep(NULL, l)
   if (is.null(cuts)) {
       for (n in seq_along(which)) {
           rect(m[which[n]] + 0.66, par("usr")[3L], m[which[n] +1] + 0.33,
                mean(rev(tree$height)[(k - 1):k]), border = border[n],
                col = color[n], ...)
           retval[[n]] <- which(cluster == as.integer(names(clustab)[which[n]]))
       }
   }
   for (n in l) {
       rect(m[which[n]] + cuts[1], cuts[2], m[which[n] +  1] + cuts[3], 
           mean(rev(tree$height)[(k - 1):k])-cuts[4], border = border[n],
           col = color[n], ...)
       retval[[n]] <- which(cluster == as.integer(names(clustab)[which[n]]))
   }
   invisible(retval)
}
