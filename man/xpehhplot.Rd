\name{xpehhplot}
\alias{xpehhplot}
\title{Plot XPEHH over a genome}
\description{Plot XPEHH over a genome.}
\usage{xpehhplot(data,plot.pval = TRUE,ylim.scan = 2,pch = 16,cex = 0.5,cex.lab = 1.25,
main = NA,cex.main = 1.5,cex.axis=1.)}

\arguments{
  \item{data}{A dataframe obtained using \code{ies2xpehh} function.}
  \item{plot.pval}{Either TRUE or FALSE if Pvalue should not be plotted}
  \item{ylim.scan}{An horizontal line is added at the corresponding coordinate, for instance to represent a significance threshold}
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

   Sabeti, P.C. et al. (2007). Genome-wide detection and characterization of positive selection in 
   human populations. \emph{Nature}, \bold{449}, 913--918.
}
\seealso{ies2xpehh}
\examples{
data(wgscan.cgu) ; data(wgscan.eut) 
## results from a genome scan (44,057 SNPs)
#see ?wgscan.eut and ?wgscan.cgu for details
res.xpehh<-ies2xpehh(wgscan.cgu,wgscan.eut,"CGU","EUT")
xpehhplot(res.xpehh)
}
