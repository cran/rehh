\name{ies2xpehh}
\alias{ies2xpehh}
\title{Compute xpEHH (standardized ratio of iES from two populations) as described in Sabeti et al. (2007)}
\description{Compute xpEHH (standardized ratio of iES from two populations) as described in Sabeti et al. (2007).}
\usage{ies2xpehh(hh_pop1,hh_pop2,popname1=NA,popname2=NA,method="bilateral")}
\arguments{
  \item{hh_pop1}{A matrix with nsnps rows and six columns (Chromosome name, position of the SNP, Frequency of the ancestral allele, iHH for the ancestral allele, iHH for the derived allele and iES) obtained after performing a scan on the first population.}
  \item{hh_pop2}{A matrix with nsnps rows and six columns (Chromosome name, position of the SNP, Frequency of the ancestral allele, iHH for the ancestral allele, iHH for the derived allele and iES) obtained after performing a scan on the second population.}
  \item{popname1}{Name of the first population compared (character string).}
  \item{popname2}{Name of the second population compared (character string).}
  \item{method}{Either "bilateral" or "unilateral". If bilateral (resp. unilateral), pvalue (assuming Rsb follows a standard Gaussian distribution under neutrality) corresponds to a bilateral (resp. unilateral) tests}
}

\details{Ratio of iES (population 1 over population 2) computed and standardized as described in Sabeti et al. (2007)}
\value{
  The returned value is a matrix with nsnps rows and four columns (Chromosome name, position of the SNP, xpEHH and Pvalue)
}
\references{
   Gautier, M. and Naves, M. (2011). Footprints of selection in the ancestral admixture of a New World
   Creole cattle breed. \emph{Molecular Ecology}, \bold{20}, 3128--3143.

   Sabeti, P.C. et al. (2007). Genome-wide detection and characterization of positive selection in 
   human populations. \emph{Nature}, \bold{449}, 913--918.
}
\seealso{calc_ehhs,scan_hh,distribplot,rsbplot}
\examples{
data(wgscan.cgu) ; data(wgscan.eut) 
## results from a genome scan (44,057 SNPs)
##see ?wgscan.eut and ?wgscan.cgu for details
xpehh<-ies2xpehh(wgscan.cgu,wgscan.eut,"CGU","EUT")
}
