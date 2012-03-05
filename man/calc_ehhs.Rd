\name{calc_ehhs}
\alias{calc_ehhs}
\title{EHHS and iES computations at a given core SNP}
\description{Compute site Extended Haplotype Homozygosity (EHHS) and integrated EHH (iES) for a given focal SNPs.}
\usage{
calc_ehhs(haplohh,mrk,limhaplo=2,limehhs=0.05,
plotehhs=TRUE,main_leg="EHHS plot")
}
\arguments{
  \item{haplohh}{An object of class \code{haplohh} (see \code{data2haplohh}).}
  \item{mrk}{Integer representing the number of the focal marker}
  \item{limhaplo}{Minimal number of haplotypes to continue computing EHHS away from the core SNP. Useless, if no missing data. However, when some data are missing, 
                  haplotypes with missing data are removed from the computation. Hence as we compute EHH further from the core SNP, less haplotypes are expected}
  \item{limehhs}{Limit below which EHHS stops to be evaluated}
  \item{plotehhs}{If TRUE, EHHS estimates are plotted for each position}
  \item{main_leg}{Legend of the EHHS plot}
}

\details{
EHHS are computed at each position upstream and downstream the focal SNP. This allows
in turn the computation of the integrated EHHS relative to map distances (iES).
}
\value{
  The returned value is a list containing the following components:
  \item{ehhs}{A vector of nsnps colums containing EHHS estimates at each chromosome position relative to the focal SNP.}
  \item{nhaplo_eval}{A matrix of two rows and nsnps colums containing the number of evaluated haplotypes at each chromosome position relative to the focal SNP for the ancestral (first row) and derived (second row) alleles.}
  \item{freqall1}{Frequency of the ancestral allele matrix for the focal SNP.}
  \item{ies}{Integrated EHHS (IES) over the chromosome.}
  }
\references{
   Gautier, M. and Naves, M. (2011). Footprints of selection in the ancestral admixture of a New World
   Creole cattle breed. \emph{Molecular Ecology}, \bold{20}, 3128--3143.

   Tang, K. and Thornton, K.R. and Stoneking, M. (2007). A New Approach for Using Genome Scans
   to Detect Recent Positive Selection in the Human Genome. \emph{Plos Biology}, \bold{7}, e171.

}
\seealso{calc_ehh,data2haplohh,scan_hh}
\examples{
#example haplohh object (280 haplotypes, 1424 SNPs)
#see ?haplohh_cgu_bta12 for details
data(haplohh_cgu_bta12) 
#computing EHH statisitics for the focal SNP at position 456
#which displays a strong signal of selection
res.ehhs<-calc_ehhs(haplohh_cgu_bta12,mrk=456)
}
