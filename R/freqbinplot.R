#'Plot of unstandardized iHS within frequency bins
#'@description Plot of unstandardized iHS within frequency bins.
#'@param x data (output of function \code{ihh2ihs})
#'@param spectrum logical. If \code{TRUE}, plot frequency spectrum instead of iHS.
#'@param main an overall title for the plot.
#'@param xlab a title for the x axis.
#'@param ylab a title for the y axis.
#'@param xlim the x coordinate range of the plot.
#'@param ylim the y coordinate range of the plot.
#'@param pch plotting 'character' see \code{points}.
#'@param ... further arguments to be passed to plot resp. points.
#'@details The plot shows the mean and the quantiles calculated by
#'function \code{\link{ihh2ihs}} for the unstandardized iHS in each frequency bin.
#'Note that the standardization of iHS is performed bin-wise in order
#'to reduce the frequency-dependence of
#'iHS values (expected under neutrality).
#'An implicit assumption of this procedure is that each bin is dominated
#'by neutral markers.
#'@seealso \code{\link{ihh2ihs}}
#'@examples library(rehh.data)
#'data(wgscan.cgu)
#'#results from a genome scan (44,057 SNPs)
#'#see ?wgscan.eut and ?wgscan.cgu for details
#'wgscan.cgu.ihs <- ihh2ihs(wgscan.cgu)
#'freqbinplot(wgscan.cgu.ihs)
#'@export
#'@importFrom graphics arrows points
#'@importFrom stats na.omit
freqbinplot <- function(x,
                        spectrum = FALSE,
                        main = NA,
                        xlab = "Derived allele frequency",
                        ylab = NA,
                        xlim = c(0, 1),
                        ylim = NULL,
                        pch = 20,
                        ...) {
  #pick table out of list, if necessary
  if (is.list(x) &
      !is.data.frame(x) & !is.null(x$frequency.class)) {
    x <- x$frequency.class
  }
  #omit NAs and standard deviation
  frequency.class <- na.omit(x[, -3])
  bin_names <- rownames(frequency.class)
  pos <- regexpr(',', bin_names)
  
  if (pos[1] == -1) {
    #discrete classes
    bin_center <- as.numeric(bin_names)
  } else{
    # continuous classes
    bin_a <- as.numeric(substr(bin_names, 2, pos - 1))
    bin_b <-
      as.numeric(substr(bin_names, pos + 1, nchar(bin_names) - 1))
    bin_center <- (bin_b + bin_a) / 2
  }
  
  lower_bound <- frequency.class[, "LOWER_QT"]
  upper_bound <- frequency.class[, "UPPER_QT"]
  
  if (spectrum) {
    if (is.na(main)) {
      main <- "Frequency spectrum in bins"
    }
    if (is.na(ylab)) {
      ylab <- "Number of markers"
    }
    if (is.null(ylim)) {
      ylim <- c(0, max(frequency.class[, "N_MRK"]))
    }
    
    plot(
      bin_center,
      frequency.class[, "N_MRK"],
      main = main,
      xlab = xlab,
      ylab = ylab,
      xlim = xlim,
      ylim = ylim,
      pch = pch,
      ...
    )
  } else{
    if (is.na(main)) {
      main <- "uniHS within frequency bins"
    }
    if (is.na(ylab)) {
      ylab <- "unstandardized iHS"
    }
    if (is.null(ylim)) {
      ylim <- c(-3.5, 3)
    }
    plot(
      NULL,
      main = main,
      xlab = xlab,
      ylab = ylab,
      xlim = xlim,
      ylim = ylim,
      ...
    )
    # error bars
    positive_length <-
      !is.na(upper_bound) &
      !is.na(lower_bound) & (upper_bound > lower_bound)
    
    arrows(
      bin_center[positive_length],
      lower_bound[positive_length],
      bin_center[positive_length],
      upper_bound[positive_length],
      code = 3,
      angle = 90,
      length = 0.2 / length(bin_center),
      # no subset!
      col = "darkgray"
    )
    
    # points on top
    points(bin_center,
           frequency.class[, 2],
           pch = pch,
           ...)
  }
}
