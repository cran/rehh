\name{scan_hh}
\alias{scan_hh}
\title{Computing EHH based statistics over a whole chromosome}
\description{Compute Extended Haplotype Homozygosity (EHH), site-specific EHH (EHHS), integrated EHH (iHH) and integrated EHHS (iES) for all SNPs of the chromosome (or linkage group).}
\usage{scan_hh(haplohh,limhaplo=2,limehh=0.05,limehhs=0.05)}
\arguments{
  \item{haplohh}{An object of class \code{haplohh} (see \code{data2haplohh}).}
  \item{limhaplo}{Minimal number of haplotypes to continue computing EHH away from the core SNP. Useless, if no missing data. However, when some data are missing, 
                  haplotypes with missing data are removed from the computation. Hence as we compute EHH further from the core SNP, less haplotypes are expected}
  \item{limehh}{Limit below which EHH stops to be evaluated}
  \item{limehhs}{Limit below which EHHS stops to be evaluated}
}

\details{
Extended Haplotype Homozygosity (EHH), site-specific EHH (EHHS), integrated EHH (iHH) and integrated EHHS (iES) are computed for all SNPs of the chromosome (or linkage group).
 This function is approximately twice as faster as a procedure calling in turn \code{calc_ehh} and \code{calc_ehhs} for all the SNP. To perform a whole genome-scan
 this function needs to be called for each chromosome and results concatenated. 
}
\value{
  The returned value is a matrix with nsnps rows and six columns (Chromosome name, position of the SNP, Frequency of the ancestral allele, iHH for the ancestral allele, iHH for the derived allele and iES)
}
\references{
   Gautier, M. and Naves, M. (2011). Footprints of selection in the ancestral admixture of a New World
   Creole cattle breed. \emph{Molecular Ecology}, \bold{20}, 3128--3143.

   Sabeti, P.C. and Reich, D.E. and Higgins, J.M. and Levine, H.Z.P and Richter, D.J. and 
   Schaffner, S.F. and Gabriel, S.B. and Platko, J.V. and Patterson, N.J. and McDonald, G.J. 
   and Ackerman, H.C. and Campbell, S.J. and Altshuler, D. and Cooper, R. and Kwiatkowski, D.
   and Ward, R. and Lander, E.S. (2002). Detecting recent positive selection in the human genome 
   from haplotype structure. \emph{Nature}, \bold{419}, 832--837.

   Tang, K. and Thornton, K.R. and Stoneking, M. (2007). A New Approach for Using Genome Scans
   to Detect Recent Positive Selection in the Human Genome. \emph{Plos Biology}, \bold{7}, e171.

   Voight, B.F. and Kudaravalli, S. and Wen, X. and Pritchard, J.K. (2006). 
   A map of recent positive selection in the human genome. \emph{Plos Biology}, \bold{4}, e72.
}
\seealso{calc_ehh,calc_ehhs,data2haplohh,ihh2ihs,ies2rsb}
\examples{
#example haplohh object (280 haplotypes, 1424 SNPs)
#see ?haplohh_cgu_bta12 for details
data(haplohh_cgu_bta12) 
res.scan<-scan_hh(haplohh_cgu_bta12)
}
