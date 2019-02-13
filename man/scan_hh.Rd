\name{scan_hh}
\alias{scan_hh}
\title{Computing EHH based statistics over a whole chromosome}
\description{Compute Extended Haplotype Homozygosity (EHH), site-specific EHH (EHHS), integrated EHH (iHH) and integrated EHHS (iES) for all SNPs of a chromosome (or linkage group).}
\usage{scan_hh(haplohh, limhaplo = 2, limehh = 0.05, limehhs = 0.05,
        scalegap = NA, maxgap = NA, 
        discard_integration_at_border = TRUE, threads = 1)}
\arguments{
  \item{haplohh}{An object of class \code{haplohh} (see \code{data2haplohh}).}
  \item{limhaplo}{Minimal number of haplotypes to continue computing EHH away from the core SNP. Useless, if no missing data. However, when some data are missing,
                  haplotypes with missing data are removed from the computation. Hence as we compute EHH further from the core SNP, less haplotypes are expected}
  \item{limehh}{Limit at which EHH stops to be evaluated}
  \item{limehhs}{Limit at which EHHS stops to be evaluated}
  \item{scalegap}{Scales gaps larger than the specified size to the specified size (default=NA, i.e. no scaling)}
  \item{maxgap}{Maximum allowed gap in bp between two SNPs below which EHH and EHHS stop to be evaluated (default=NA, i.e., no limitation)}
  \item{discard_integration_at_border}{If TRUE and if first or last marker or a gap (larger than \code{maxgap}) is reached and EHH(S) is greater than \code{limehh(s)}, then iHH/IES is set to \code{NA}}
  \item{threads}{Number of threads to parallelize compuation}
}

\details{
Extended Haplotype Homozygosity (EHH), site-specific EHH (EHHS), integrated EHH (iHH) and integrated EHHS (iES) are computed for all SNPs of the chromosome (or linkage group).
 This function is several times faster as a procedure calling in turn \code{calc_ehh} and \code{calc_ehhs} for all SNPs. To perform a whole genome-scan
 this function needs to be called for each chromosome and the results concatenated.
}
\value{
  The returned value is a dataframe with \code{haplohh@nsnps} rows and seven columns (Chromosome name, position of the SNP, Frequency of the ancestral allele, iHH for the ancestral allele, iHH for the derived allele, iES using the estimator by Sabeti et al. (2007) estimator and iES using the estimator by Tang et al. (2007))
}
\references{
   Gautier, M. and Naves, M. (2011). Footprints of selection in the ancestral admixture of a New World
   Creole cattle breed. \emph{Molecular Ecology}, \bold{20}, 3128--3143.

   Sabeti, P.C. et al. (2002). Detecting recent positive selection in the human genome
   from haplotype structure. \emph{Nature}, \bold{419}, 832--837.

   Sabeti, P.C. et al. (2007). Genome-wide detection and characterization of positive selection in
   human populations. \emph{Nature}, \bold{449}, 913--918.

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
