% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/freqbinplot.R
\name{freqbinplot}
\alias{freqbinplot}
\title{Plot of unstandardized iHS within frequency bins}
\usage{
freqbinplot(
  x,
  spectrum = FALSE,
  main = NA,
  xlab = "Derived allele frequency",
  ylab = NA,
  xlim = c(0, 1),
  ylim = NULL,
  pch = 20,
  ...
)
}
\arguments{
\item{x}{data (output of function \code{ihh2ihs})}

\item{spectrum}{logical. If \code{TRUE}, plot frequency spectrum instead of iHS.}

\item{main}{an overall title for the plot.}

\item{xlab}{a title for the x axis.}

\item{ylab}{a title for the y axis.}

\item{xlim}{the x coordinate range of the plot.}

\item{ylim}{the y coordinate range of the plot.}

\item{pch}{plotting 'character' see \code{points}.}

\item{...}{further arguments to be passed to plot resp. points.}
}
\description{
Plot of unstandardized iHS within frequency bins.
}
\details{
The plot shows the mean and the quantiles calculated by
function \code{\link{ihh2ihs}} for the unstandardized iHS in each frequency bin.
Note that the standardization of iHS is performed bin-wise in order
to reduce the frequency-dependence of
iHS values (expected under neutrality).
An implicit assumption of this procedure is that each bin is dominated
by neutral markers.
}
\examples{
library(rehh.data)
data(wgscan.cgu)
#results from a genome scan (44,057 SNPs)
#see ?wgscan.eut and ?wgscan.cgu for details
wgscan.cgu.ihs <- ihh2ihs(wgscan.cgu)
freqbinplot(wgscan.cgu.ihs)
}
\seealso{
\code{\link{ihh2ihs}}
}
