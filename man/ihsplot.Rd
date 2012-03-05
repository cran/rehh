\name{ihsplot}
\alias{ihsplot}
\title{Plot iHS over a genome}
\description{Plot iHS over a genome.}
\usage{ihsplot(data,plot.pval="TRUE",ylim.scan=2,pch=16,main="iHS")}
\arguments{
  \item{data}{A matrix with nsnps rows and 4 columns (Chromosome name, position of the SNP, iHS, Pvalue) (e.g. obtained using \code{ihh2ihs} function.}
  \item{plot.pval}{Either TRUE or FALSE if Pvalue should not be plotted}
  \item{ylim.scan}{An horizontal line is added at the corresponding coordinate, for instance to represent a significance threshold}
  \item{pch}{Type of the points representing SNPs in the plot(s)}
  \item{main}{Main Legend of the plot}
}
\value{The function returns a plot}
\references{
   Gautier, M. and Naves, M. (2011). Footprints of selection in the ancestral admixture of a New World
   Creole cattle breed. \emph{Molecular Ecology}, \bold{20}, 3128--3143.
}
\seealso{calc_ehh,scan_hh,distribplot,ihsplot}
\examples{
data(wgscan.cgu) 
## results from a genome scan (44,057 SNPs)
##see ?wgscan.eut and ?wgscan.cgu for details
res.ihs<-ihh2ihs(wgscan.cgu)
ihsplot(res.ihs$res.ihs)
}
