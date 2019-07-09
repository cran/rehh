#'Plot distribution of standardized iHS, Rsb or XP-EHH values
#'@description Plot the observed distribution of standardized iHS, Rsb or XP-EHH values together with
#'the standard Gaussian distribution.
#'@param data a vector of iHS, Rsb or XPEHH values.
#'@param col a vector describing the colors of the observed and Gaussian distribution, respectively.
#'@param main an overall title for the plot.
#'@param xlab a title for the x axis.
#'@param lty line type.
#'@param lwd line width.
#'@param qqplot logical. If \code{TRUE} a qq-plot is drawn instead of the distribution density curve.
#'@param pch point character for qqplot (see \code{\link[graphics]{points}}).
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
           main = "Genome-wide distribution",
           xlab = "",
           qqplot = FALSE,
           pch = 20,
           ...) {
    if (!qqplot) {
      plot(
        density(data, na.rm = TRUE),
        main = main,
        xlab = xlab,
        col = col[1],
        lty = lty,
        lwd = lwd,
        ...
      )
      curve(dnorm, col = col[2], add = TRUE)
      legend(
        "topright",
        c("Observed", "Gaussian"),
        bty = "n",
        col = col,
        lty = lty,
        lwd = lwd
      )
    } else {
      qqnorm(data[!is.na(data)],
             pch = pch,
             ...)
      abline(a = 0, b = 1, lty = 2)
    }
  }
