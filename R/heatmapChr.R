#' @rdname heatmapChr
#' @title Heatmap of GRanges Object
#' @description A function to create a heatmap or a graphical representation of
#'     data where the individual values contained in a matrix are represented as
#'     colors of the GRanges object.
#' @details This function creates a heatmap is a false color image with a color
#'     scale added to the right side and a chromosome scale to the bottom.
#' @param GR A Grange object with the variable of interest in its metacolumn.
#' @param filename This is the name of the image file in which we want to output
#'     the heatmap.
#' @param chr This is a required argument which corresponds to the chromosome of
#'     interest in the data.
#' @param sample.id This is the id or column name which is the sample of the
#'     GRange object.
#' @param factor.scale A number to scale chromosome position into "bp", "kbp",
#'     or "Mbp", depending on chromosome size
#' @param absolute If absolute == TRUE, all the values taken for the heatmap
#'     would be absolute values of the GRange object.
#' @param xtitle This is the x axis title for the heatmap which will be produced..
#' @param Barpalette This is a required argument which defines the barpalette
#'     for the heatmap which can be a colorRampPalette object.
#' @param format This is the format of the output file which will have the
#'     heatmap. Possible formats are "jpg", "png", "tiff" and "pdf". The default
#'     value for this is "tiff with compression = "lzw" and res = 600.
#' @param width This is the width of the heatmap image which will be produced.
#' @param height This is the height of the heatmap image which will be produced.
#' @param fontfamily This value defines the name of the font family which will
#'     be used in text or labels for the heatmap. Standard values are "serif",
#'     "sans" and "mono", and the Hershey font families are also available.
#' @param font value to pass in to parameter 'font' from function
#'     \code{\link[graphics]{par}}. An integer which specifies which font to use
#'     for text. If possible, device drivers arrange so that 1 corresponds to
#'     plain text (the default), 2 to bold face, 3 to italic and 4 to bold
#'     italic.
#' @param mar.scale A numeric vector of length 4, which sets the margin sizes in
#'     the following order: bottom, left, top, and right.
#' @param oma,oma.scale Same as 'oma' graphical parameter (see
#'     \code{\link[graphics]{par}}). 'oma' is used in the heatmap and oma.scale
#'     in the color bar.
#' @param mgp.scale A numeric vector of length 3, which sets the axis label
#'     locations relative to the edge of the inner plot window. The first value
#'     represents the location the labels (i.e. xlab and ylab in plot), the
#'     second the tick-mark labels, and third the tick marks.
#' @param compression This is the the type of compression to be used.
#'     The default compression type is "lzw".
#' @param res This is the nominal resolution in ppi which will be recorded in
#'     the bitmap file, if a positive integer. Also used for units other than
#'     the default, and to convert points to pixels. The defualt resolution is
#'     300ppi.
#' @param pointsize The pointsize of plotted text, interpreted as big points
#'     (1/72 inch) at res ppi. The default value for this is 1.
#' @param col.bar.lwd Line width grphical parameter for the color bar.
#' @param cex As in \code{\link[graphics]{par}}, it is numerical value giving
#'     the amount by which plotting text and symbols should be magnified
#'     relative to the default. This starts as 1 when a device is opened, and is
#'      reset when the layout is changed, e.g. by setting mfrow.
#' @param cex.xaxis,cex.yaxis,cex.xtitle,cex.bar.lab Cex values for x-axis,
#'      y-axis, x-title, and color-bar labels, respectively. As in
#'      \code{\link[graphics]{par}}, these are the magnification to be used for
#'      x-axis, y-axis, and label annotations relative to the current setting of
#'       cex.
#' @param lwd.ticks Line width for axis and ticks (heatmap only).
#' @param xaxis.labels.pos Y-coordinate for x-axis labels.
#' @param xaxis.adj Adjustment of the x-axis labels.
#' @param tick.breaks An integer number used to introduce the number breaks in
#'     the chromosome scale where the tick will be located.
#' @param line.xtitle specifying a value for xline.label overrides the default
#'     placement of x-axis title, and places them this many lines outwards from
#'     the plot edge
#' @param jpg.type Paramter 'type' from 'jpeg' functions (see ?jpeg).
#' @param ylas,bar.las numeric in {0,1,2,3}; the style of y-axis and colo-bar
#'     labels, as given for graphical parameters "las" (see
#'     \code{\link[graphics]{par}}).
#' @param ... Additional graphical parameters for 'par' R function used in the
#'     heatmap (not in the color bar).
#' @return A GRanges object with the new genomic regions and their corresponding
#'     summarized statistic.
#'
#' @examples
#' data(GRMatrix)
#' ## The pallette used in the bar color
#' bar.palette <- colorRampPalette(c(rep("cyan",4), "green",rep("yellow", 2),
#'                                   rep("red", 3), rep("darkblue", 2),
#'                                   rep("black",2)), bias = 1.1, space = "rgb")
#'
#' ## heatmap construction
#' k <- "Chr1"
#' file <- paste0(getwd(), "/heatmap_", k)
#' xlab <- paste0("CG. Chromosome ", sub("Chr", "", k), " (Mbp)")
#' heatmapChr(GR = GRMatrix, filename = file, format = "tiff", 
#'            chr = k, xtitle = xlab,
#'            Barpalette = bar.palette, mar = c(4, 6, 0, 0), res = 600,
#'            height = 350, width =2500, font = 2, fontfamily = "serif", ylas=1,
#'            lwd = 0.1, mar.scale = c(4, 5, 0.5, 3), cex.xtitle = 2.5,
#'            cex.bar.lab = 2, mgp = c(3,1,0), pointsize = 2, cex.xaxis = 2,
#'            cex.yaxis = 2, oma=c(1,1,1,0), xaxis.adj = c(0.5, 0.7),
#'            lwd.ticks = 0.1, line.xtitle = 3, mgp.scale = c(3, 1, 0),
#'            col.bar.lwd = 0.1, factor.scale = 10^3,
#'            sample.id = c("CT1", "CT2", "CT3", "TT1", "TT2", "TT3"))
#' ## To remove the file containing the heatmap
#' file.remove(paste0(file, ".tiff"))
#' 
#' @importFrom GenomicRanges GRanges
#' @importFrom grDevices dev.off heat.colors jpeg pdf png tiff
#' @importFrom graphics axis box image layout par plot polygon text title
#' @importFrom stats hclust
#' @export
#' @author Robersy Sanchez

heatmapChr <- function(GR, filename=NULL, chr, sample.id=NULL,
                       factor.scale=10^6, absolute=FALSE, xtitle=NULL,
                       Barpalette, format="tiff", width=4000, height=790,
                       fontfamily="sans", font=2,
                       mar.scale=c(2,2,2,2), mgp.scale=c(3, 1, 0),
                       mar.heatmap = c(2,2,2,2), mgp.heatmap = c(3.5,1,0),
                       compression="lzw", res=900, pointsize=1,
                       col.bar.lwd=1, cex=1, cex.xaxis=1.6,
                       cex.yaxis=2, cex.xtitle=2, cex.bar.lab=2, lwd.ticks=0.5,
                       xaxis.labels.pos=0.1, oma=c(2,2,2,2),
                       oma.scale = c(0, 0, 0, 0),
                       xaxis.adj=c(0.5, 1), tick.breaks=500, line.xtitle=NA,
                       jpg.type=c("cairo", "cairo-png", "Xlib", "quartz"),
                       ylas=1, bar.las=1, ...) {

   if(!inherits(GR, "GRanges")) stop("GR object must be a 'GRanges' object")
   seqlevels(GR, pruning.mode="coarse") <- chr
   if (is.null(sample.id)) {
       warning("'sample.id' is not provided","\n",
           "So, all columns from the GR meta-colums are used in the heatmap")
       sample.id <- try(colnames(mcols(GR)))
   }

   midpoints <- (start(GR) + end(GR)) / 2

   d <- mcols(GR[, sample.id])
   d <- as.matrix(d)
   if (absolute) d <- abs(d)

   breaks <- seq(min( d, na.rm=TRUE), max(d, na.rm=TRUE), length.out=tick.breaks)
   ind.ticks <- seq(1, nrow(d), tick.breaks)

   labels <- format(round(midpoints[ind.ticks] / factor.scale, digits=1),
                   digits=1)

   if (is.null(filename)) {
       filename <- paste0("heatmap_", chr, ".", format)
   } else filename <- paste0(filename, ".", format)

   if (is.null(xtitle)) paste("Chromosome ", chr, sep="")

   if (format == "jpeg") {
      jpeg(filename = filename, width = width, height = height, res = res,
           pointsize = pointsize, type = jpg.type[1])
   } else if (format == "png") {
       png(filename = filename, width = width, height = height,
           compression = compression, res = res, pointsize = pointsize)
   } else if (format == "pdf") {
       pdf(file = filename, width = width, height = height,
           pointsize = pointsize)
   } else if (format == "tiff") {
       tiff(filename = filename, width = width, height = height,
            compression = compression, res = res, pointsize = pointsize)
   }

   layout(matrix(c(1, 2), nrow = 1, ncol = 2), widths = c(4, 0.38), 
         heights = c(4 , 4))

   par(mar = mar.heatmap, mgp = mgp.heatmap, lwd = lwd.ticks, family = 
         fontfamily, oma = oma, cex = cex)
   image(1:nrow(d), 1:ncol(d), d, axes=F, ylab="", xlab="",
           col=Barpalette(length(breaks) - 1), lwd=lwd.ticks,
           family=fontfamily)
   axis(1, at = ind.ticks, labels = FALSE, font = font, cex.axis = cex.xaxis,
         lwd = lwd.ticks, family = fontfamily)
   axis(2, 1:ncol(d), labels = sample.id, font = font, cex.axis = cex.yaxis,
           lwd = lwd.ticks, las = ylas, family = fontfamily)
   par(family = fontfamily, cex = cex, cex.lab = cex.xtitle, font.lab = font)
   title(xlab = xtitle, line = line.xtitle)
   text(ind.ticks, y = xaxis.labels.pos, srt = 0, labels = labels, xpd = TRUE,
       font = font, cex = cex.xaxis, adj = xaxis.adj, lwd = lwd.ticks,
       family = fontfamily)
   box()
   par(mar = mar.scale, family = fontfamily, cex = cex, cex.axis = cex.bar.lab,
       mgp = mgp.scale, font = font, lwd = col.bar.lwd)
   .colorBar(d, col=Barpalette(length(breaks) - 1), breaks = breaks, 
            horiz = FALSE, lwd = col.bar.lwd, cex.axis = cex.bar.lab, 
            family = fontfamily, font = font, las = bar.las)
   box()
   dev.off()
   message("OK. The heatmap file is located at:")
   message(filename)
}


### ========================= Auxiliary function ============================

#' @name .colorBar
#' @title Color bar for heatmaps
#' @description This function is for internal use. It creates a color scaled
#'     bar for heatmap.
#' @param z matrix of values used to visualize the heatmap
#' @param zlim Limit for the numerical(color) scale, which must be consistent
#'     with parameter 'break'.
#' @param col Palette of color to use
#' @param breaks Numerical vector with the breaks used to produce the heatmap
#' @param horiz Whether to the color bar will be horizontal(= TRUE) or
#'     vertical(= FALSE)
#' @param ylim User defined limits for y-axis. Depending on the orientation, x-
#'     or y-limits may be defined that are different from the z-limits and will
#'     reduce the range of colors displayed.
#' @param xlim User defined limits for x-axis
#' @param lwd Line width.
#' @param cex.axis Cex values for color bar axis labels.
#' @param ... Additional parameter to pass to 'par' R function
#' @keywords internal
#' @return Image with color bar
#'
.colorBar <- function(z, zlim, col=heat.colors(12), breaks, horiz=TRUE,
                      ylim=NULL, xlim=NULL, lwd=0.5, cex.axis=1, ...){
   
   if (!missing(breaks)) {
      if (length(breaks) != (length(col) + 1)) {
         stop("must have one more break than colour")
      }
   }
   if (missing(breaks) & !missing(zlim)) {
      breaks <- seq(zlim[1], zlim[2], length.out=(length(col) + 1))
   }
   if (missing(breaks) & missing(zlim)) {
      zlim <- range(z, na.rm=TRUE)
      ## adds a bit to the range in both directions
      zlim[2] <- zlim[2] + c(zlim[2] - zlim[1]) * (1E-3)
      zlim[1] <- zlim[1] - c(zlim[2] - zlim[1]) * (1E-3)
      breaks <- seq(zlim[1], zlim[2], length.out=(length(col) + 1))
   }
   
   poly <- vector(mode="list", length(col))
   for (i in seq(poly)) {
      poly[[i]] <- c(breaks[i], breaks[i + 1], breaks[i + 1], breaks[i])
   }
   
   if (horiz) {YLIM <- c(0,1); XLIM <- range(breaks)}
   if (!horiz) {YLIM <- range(breaks); XLIM <- c(0, 1)}
   if (missing(xlim)) xlim <- XLIM
   if (missing(ylim)) ylim <- YLIM
   
   plot(x = 0.5, y = 0.5, type = "n", ylim = ylim, xlim = xlim, xaxt = "n",
        yaxt = "n", xaxs = "i", yaxs = "i", xlab = "", ylab = "",
        lwd = lwd, ...)
   if (!horiz) axis(side = 2, ylim = ylim, xlim = xlim, lwd = lwd,
                  cex.axis=cex.axis, ... )
   if (horiz) axis(side = 1, ylim = ylim, xlim = xlim, lwd = lwd,
                  cex.axis = cex.axis, ... )
   
   for (i in seq(poly)) {
      if (horiz) {
         polygon(poly[[i]], c(0,0,1,1), col = col[i], border = NA)
      }
      if (!horiz) {
         polygon(c(0,0,1,1), poly[[i]], col = col[i], border = NA)
      }
   }
}


