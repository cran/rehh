#'Plot EHHS around a focal marker
#'@description Plot curve of EHHS values around a focal marker.
#'@param x data (output of \code{\link{calc_ehhs}}).
#'@param nehhs logical. If \code{TRUE}, plot normalized EHHS.
#'@param ylim the y limits of the plot
#'@param type plot type (see \code{\link[graphics]{matplot}}).
#'@param main title for the plot (default \code{NA}, i.e. none).
#'@param xlab title for the x-axis.
#'@param ylab title for the y-axis.
#'@param bty box type around plot (see \code{\link[graphics]{par}}).
#'@param mrk.col color of the vertical line at the focal marker position.
#'@param ... further arguments to be passed to function \code{\link[graphics]{plot}}.
#'@seealso \code{\link{data2haplohh}}, \code{\link{plot.ehh}}, \code{\link{calc_ehhs}}, \code{\link{scan_hh}}.
#'@examples
#'#example haplohh object (280 haplotypes, 1424 SNPs)
#'#see ?haplohh_cgu_bta12 for details
#'data(haplohh_cgu_bta12)
#'#computing EHHS statisitics for the marker "F1205400"
#'#which displays a strong signal of selection
#'ehhs <- calc_ehhs(haplohh_cgu_bta12, mrk = "F1205400")
#'plot(ehhs)
#'@export
plot.ehhs <- function(x,
                      nehhs = FALSE,
                      ylim = c(0, 1),
                      type = "l",
                      main = paste0("EHHS around '", x$mrk.name, "'"),
                      xlab = "Position",
                      ylab = "Extended Haplotype Homozygosity per Site",
                      bty = "n",
                      mrk.col = "gray",
                      ...) {
  # check parameters
  if (is.null(x$ehhs)) {
    stop("Data does not contain a data frame with EHHS values.", call. = FALSE)
  }
  if (is.null(x$ehhs$POSITION)) {
    stop("Data frame seems not to contain a column with marker positions.",
         call. = FALSE)
  }
  if (nrow(x$ehhs) == 0) {
    stop("Empty data frame.", call. = FALSE)
  }
  if (is.null(x$mrk.name)) {
    stop("No marker id found.", call. = FALSE)
  }
  foc_pos <- x$ehhs[x$mrk.name, "POSITION"]
  if (is.na(foc_pos)) {
    stop(paste0(
      "No position information found for marker '",
      x$mrk.name,
      "'.",
      call. = FALSE
    ))
  }
  ehhs <- x$ehhs[[ifelse(nehhs, "NEHHS", "EHHS")]]
  if (length(ehhs) == 0) {
    stop("Data frame seems not to contain a column with EHHS values.",
         call. = FALSE)
  }
  
  # perform plot
  p <- floor(log(max(x$ehhs$POSITION), 1000))
  ## only shrink big scales, but never magnify small ones (p<0)
  scale <- 1000 ** max(0, p)
  ## no unit if p < 0
  unit <- c("", "(bp)", "(kb)", "(Mb)", "(Gb)")[max(-1, p)  + 2]

  dot_args <- list(...)
  if (!is.null(dot_args$xlim)) {
    dot_args$xlim <- dot_args$xlim / scale
  }
  
  do.call("plot", c(
    list(
      x = x$ehhs$POSITION / scale,
      y = ehhs,
      ylim = ylim,
      type = type,
      main = main,
      xlab = paste(xlab, unit),
      ylab = ylab,
      bty = "n"
    ),
    dot_args
  ))
  
  
  abline(v = foc_pos / scale,
         lty = 2,
         col = mrk.col)

}
