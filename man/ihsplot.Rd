\name{ihsplot}
\alias{ihsplot}
\title{Plot iHS over a genome}
\description{Plot iHS over a genome.}
\usage{ihsplot(ihsdata, plot.pval = TRUE, ylim.scan = 2, pch = 16, cex = 0.5,
            cex.lab = 1.25, main = NA, cex.main = 1.5, cex.axis=1.)}
\arguments{
  \item{ihsdata}{A list obtained with the \code{ihh2ihs} function.}
  \item{plot.pval}{Either TRUE or FALSE if Pvalue should not be plotted}
  \item{ylim.scan}{An horizontal line is added at the corresponding coordinate, for
 instance to represent a significance threshold}
  \item{pch}{Type of the points representing SNPs in the plot(s)}
  \item{cex}{Size of the points representing SNPs in the plot(s)}
  \item{cex.lab}{Size of axis legends}
  \item{main}{Main Legend of the plot}
  \item{cex.main}{Size of the main legend}
  \item{cex.axis}{Size of the axis annotations}
}
\value{The function returns a plot}
\references{
   Gautier, M. and Naves, M. (2011). Footprints of selection in the ancestral admixture of a New World
   Creole cattle breed. \emph{Molecular Ecology}, \bold{20}, 3128--3143.

   Voight, B.F. and Kudaravalli, S. and Wen, X. and Pritchard, J.K. (2006).
   A map of recent positive selection in the human genome. \emph{Plos Biology}, \bold{4}, e72.
}
\seealso{ihh2ihs}
\examples{
data(wgscan.cgu)
## results from a genome scan (44,057 SNPs)
##see ?wgscan.eut and ?wgscan.cgu for details
res.ihs<-ihh2ihs(wgscan.cgu)
ihsplot(res.ihs)
}
