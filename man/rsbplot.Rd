\name{rsbplot}
\alias{rsbplot}
\title{Plot Rsb over a genome}
\description{Plot Rsb over a genome.}
\usage{rsbplot(data,plot.pval="TRUE",ylim.scan=2,pch=16,main=NA)}

\arguments{
  \item{data}{A matrix with nsnps rows and 4 columns (Chromosome name, position of the SNP, Rsb, Pvalue (either unilateral or bilateral)) (e.g. obtained using \code{ies2rsb} function.}
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
\seealso{calc_ehhs,scan_hh,distribplot,ihsplot}
\examples{
data(wgscan.cgu) ; data(wgscan.eut) 
## results from a genome scan (44,057 SNPs)
#see ?wgscan.eut and ?wgscan.cgu for details
res.rsb<-ies2rsb(wgscan.cgu,wgscan.eut,"CGU","EUT")
rsbplot(res.rsb$res.rsb)
}
