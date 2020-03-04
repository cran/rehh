#'Plot distribution of standardized iHS, Rsb or XP-EHH values
#'@description Plot the observed distribution of standardized iHS, Rsb or XP-EHH values together with
#'the standard Gaussian distribution.
#'@param data a vector of iHS, Rsb or XPEHH values.
#'@param col a vector describing the colors of the observed and Gaussian distribution, respectively.
#'@param lty line type.
#'@param lwd line width.
#'@param qqplot logical. If \code{TRUE} a qq-plot is drawn instead of the distribution density curve.
#'@param resolution affects only qqplot. Rasterize data points to a quadratic grid with the specified resolution and remove
#'duplicate points. Defaults to 0.01.
#'@param ... further arguments passed to \code{\link[graphics]{plot.default}}.
#'@return The function returns a plot.
#'@seealso \code{\link{ihh2ihs}}, \code{\link{ines2rsb}}, \code{\link{ies2xpehh}}, \code{\link{manhattanplot}}.
#'@examples library(rehh.data)
#'#results from a genome scan (44,057 SNPs) see ?wgscan.cgu for details
#'data(wgscan.cgu)
#'#extract vector with iHS values from data frame
#'IHS <- ihh2ihs(wgscan.cgu)$ihs[["IHS"]]
#'distribplot(IHS, main = "iHS (CGU population)")
#'distribplot(IHS, main = "iHS (CGU population)", qqplot = TRUE)
#'@export
#'@importFrom graphics abline curve legend par plot
#'@importFrom stats density dnorm pnorm qqnorm
distribplot <-
  function(data,
           lty = 1,
           lwd = 1.5,
           col = c("blue", "red"),
           qqplot = FALSE,
           resolution = 0.01,
           ...) {
    if (qqplot &
        (!is.numeric(resolution) |
         length(resolution) != 1 | resolution <= 0)) {
      stop("Resolution has to be specified by a positive number.",
           call. = FALSE)
    }
    
    dot.args <- list(...)
    
    if (!qqplot) {
      if (is.null(dot.args$main)) {
        dot.args$main <- "Genome-wide distribution"
      }
      if (is.null(dot.args$xlab)) {
        dot.args$xlab <- ""
      }
      
      do.call(plot,
              c(list(
                density(data, na.rm = TRUE),
                col = col[1],
                lwd = lwd
              ),
              dot.args))
      
      curve(dnorm,
            col = col[2],
            lty = lty,
            add = TRUE)
      
      legend(
        "topright",
        c("Observed", "Gaussian"),
        bty = "n",
        col = col,
        lty = lty,
        lwd = lwd
      )
    } else {
      coord <-
        unique(round(as.data.frame(qqnorm(data, plot.it = FALSE)) / resolution)) *
        resolution
      
      if (is.null(dot.args$main)) {
        dot.args$main <- "Normal Q-Q Plot"
      }
      if (is.null(dot.args$xlab)) {
        dot.args$xlab <- "Theoretical Quantiles"
      }
      if (is.null(dot.args$ylab)) {
        dot.args$ylab <- "Sample Quantiles"
      }
      
      do.call(plot,
              c(list(x = coord$x,
                     y = coord$y),
                dot.args))
      
      abline(a = 0, b = 1, lty = 2)
    }
  }
