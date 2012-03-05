\name{calc_ehh}
\alias{calc_ehh}
\title{EHH and iHH computations at a given core SNP}
\description{Compute Extended Haplotype Homozygosity (EHH) and integrated EHH (iHH) for a given focal SNPs.}
\usage{
calc_ehh(haplohh,mrk,limhaplo=2,limehh=0.05,
plotehh=TRUE,main_leg="EHH plot")
}
\arguments{
  \item{haplohh}{An object of class \code{haplohh} (see \code{data2haplohh}).}
  \item{mrk}{Integer representing the number of the focal marker}
  \item{limhaplo}{Minimal number of haplotypes to continue computing EHH away from the core SNP. Useless, if no missing data. However, when some data are missing, 
                  haplotypes with missing data are removed from the computation. Hence as we compute EHH further from the core SNP, less haplotypes are expected}
  \item{limehh}{Limit below which EHH stops to be evaluated}
  \item{plotehh}{If TRUE, EHH estimates for both the ancestral and derived allele are plotted for each position}
  \item{main_leg}{Legend of the EHH plot.}
}

\details{
EHH are computed at each position upstream and downstream the focal SNP for both the derived and ancestral allele. This allows
in turn the computation of the integrated EHH relative to map distances (iHH).
}
\value{
  The returned value is a list containing the
  following components:
  \item{ehh}{A matrix of two rows and nsnps colums containing EHH estimates at each chromosome position relative to the focal SNP for the ancestral (first row) and derived (second row) alleles.}
  \item{nhaplo_eval}{A matrix of two rows and nsnps colums containing the number of evaluated haplotypes at each chromosome position relative to the focal SNP for the ancestral (first row) and derived (second row) alleles.}
  \item{ihh}{A vector of two elements corresponding respectively to the iHH (integrated EHH) for the ancestral and derived allele.}
}
\references{
   Gautier, M. and Naves, M. (2011). Footprints of selection in the ancestral admixture of a New World
   Creole cattle breed. \emph{Molecular Ecology}, \bold{20}, 3128--3143.

   Sabeti, P.C. and Reich, D.E. and Higgins, J.M. and Levine, H.Z.P and Richter, D.J. and 
   Schaffner, S.F. and Gabriel, S.B. and Platko, J.V. and Patterson, N.J. and McDonald, G.J. 
   and Ackerman, H.C. and Campbell, S.J. and Altshuler, D. and Cooper, R. and Kwiatkowski, D.
   and Ward, R. and Lander, E.S. (2002). Detecting recent positive selection in the human genome 
   from haplotype structure. \emph{Nature}, \bold{419}, 832--837.

   Voight, B.F. and Kudaravalli, S. and Wen, X. and Pritchard, J.K. (2006). 
   A map of recent positive selection in the human genome. \emph{Plos Biology}, \bold{4}, e72.
}
\seealso{calc_ehhs,data2haplohh,scan_hh}
\examples{
#example haplohh object (280 haplotypes, 1424 SNPs)
#see ?haplohh_cgu_bta12 for details
data(haplohh_cgu_bta12) 

#computing EHH statisitics for the focal SNP at position 456
# which displays a strong signal of selection
res.ehh<-calc_ehh(haplohh_cgu_bta12,mrk=456)
}
