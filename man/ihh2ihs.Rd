\name{ihh2ihs}
\alias{ihh2ihs}
\title{Compute iHS (standardized iHH)}
\description{Compute iHS (standardized iHH).}
\usage{ihh2ihs(res_ihh,freqbin=0.025,minmaf=0.05)}
\arguments{
  \item{res_ihh}{A matrix with nsnps rows and six columns (Chromosome name, position of the SNP, Frequency of the ancestral allele, iHH for the ancestral allele, iHH for the derived allele and iES) obtained after performing a scan on the population of interest.}
  \item{freqbin}{Size of the bin to standardize log(iHH1/iHH2) according to the underlying Derived Allele frequency. Allele frequency bins vary from \code{minmaf} to 1-\code{minmaf} per step of size \code{freqbin}. If \code{freqbin} is set to 0 (e.g. in the case of a large number of SNPs and few haplotypes), standardization is performed
 considering each observed frequency as a frequency class.
}
  \item{minmaf}{SNPs with a MAF (Minor Allele Frequency) lower than \code{minmaf} are discarded from the analysis}
}
\details{iHS (standardized iHH) are standardized as described in Voight et al. (2006)}
\value{
  The returned value is a list containing two elements
  \item{res.ihs}{matrix with nsnps rows and four columns (Chromosome name, position of the SNP, Rsb and Pvalue)}
  \item{summary.class}{matrix with nclasses rows and four columns (Class, Number of SNPs belonging to this class,position of the SNP, mean iHH in this class, standard deviation of iHH in this class)}
}
\references{
   Gautier, M. and Naves, M. (2011). Footprints of selection in the ancestral admixture of a New World
   Creole cattle breed. \emph{Molecular Ecology}, \bold{20}, 3128--3143.

   Voight, B.F. and Kudaravalli, S. and Wen, X. and Pritchard, J.K. (2006). 
   A map of recent positive selection in the human genome. \emph{Plos Biology}, \bold{4}, e72.
}
\seealso{calc_ehh,scan_hh,distribplot,ihsplot}
\examples{
data(wgscan.cgu) 
## results from a genome scan (44,057 SNPs)
##see ?wgscan.eut and ?wgscan.cgu for details
res.ihs<-ihh2ihs(wgscan.cgu)
}
