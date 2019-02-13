\name{calc_ehhs}
\alias{calc_ehhs}
\title{EHHS and iES computations at a given core SNP}
\description{Compute site-specific Extended Haplotype Homozygosity (EHHS) and integrated EHHS (iES) for a given focal SNP.}
\usage{
calc_ehhs(haplohh, mrk, limhaplo = 2, limehhs = 0.05,
    scalegap = NA, maxgap = NA,
    discard_integration_at_border = TRUE, plotehhs = TRUE,
    lty = 1, lwd = 1.5, col = c("blue", "red"), xlab = "Position",
    ylab = expression(Site ~ specific ~ italic(EHH) ~ (italic(EHHS))),
    cex.lab = 1.25, main = NA, cex.main = 1.5)
}
\arguments{
  \item{haplohh}{An object of class \code{haplohh} (see \code{data2haplohh}).}
  \item{mrk}{Either the number of the marker in the haplohh object (as integer) or its name (as string) to specify the focal marker}
  \item{limhaplo}{Minimal number of haplotypes to continue computing EHHS away from the core SNP. Useless, if no missing data. However, when some data are missing,
                  haplotypes with missing data are removed from the computation. Hence as we compute EHH further from the core SNP, less haplotypes are expected}
  \item{limehhs}{Limit at which EHHS stops to be evaluated}
  \item{scalegap}{Scales gaps larger than the specified size to the specified size (default=NA, i.e. no scaling)}
  \item{maxgap}{Maximum allowed gap in bp between two SNPs below which EHHS stops to be evaluated (default=NA, i.e., no limitation)}
  \item{discard_integration_at_border}{If TRUE and if first or last marker or a gap (larger than \code{maxgap}) is reached and EHHS is greater than \code{limehh}, then IES is set to \code{NA}}
  \item{plotehhs}{If TRUE, EHHS estimates are plotted for each position}
  \item{lty}{Line type for the EHHS_Sabeti_et_al_2007 and EHHS_Tang_et_al_2007 (respectively) curves}
  \item{lwd}{Line width for the EHHS_Sabeti_et_al_2007 and EHHS_Tang_et_al_2007 (respectively) curves}
  \item{col}{Color for the EHHS_Sabeti_et_al_2007 and EHHS_Tang_et_al_2007 (respectively) curves}
  \item{xlab}{Legend for the x--axis}
  \item{ylab}{Legend for the y--axis}
  \item{cex.lab}{Size of the axis legend}
  \item{main}{Main legend of the EHHS plot}
  \item{cex.main}{Size of the main legend}
}

\details{
EHHS are computed at each position upstream and downstream of the focal SNP. This allows
in turn the computation of the integrated EHHS relative to map distances (iES).
}
\value{
  The returned value is a list containing the following elements:
  \item{EHHS_Tang_et_al_2007}{A vector of nsnps columns containing EHHS estimates at each chromosome position relative to the focal SNP computed as described in the Tang et al. (2007).}
  \item{EHHS_Sabeti_et_al_2007}{A vector of nsnps columns containing EHHS estimates at each chromosome position relative to the focal SNP computed as described in the Sabeti et al. (2007).}
  \item{nhaplo_eval}{A matrix of two rows and nsnps columns containing the number of evaluated haplotypes at each chromosome position relative to the focal SNP for the ancestral (first row) and derived (second row) alleles.}
  \item{IES_Sabeti_et_al_2007}{Integrated EHHS (computed using the estimator by Sabeti et al. (2007)) over the chromosome.}
  \item{IES_Tang_et_al_2007}{Integrated EHHS (computed using the estimator by Tang et al. (2007)) over the chromosome.}
  }
\references{
   Gautier, M. and Naves, M. (2011). Footprints of selection in the ancestral admixture of a New World
   Creole cattle breed. \emph{Molecular Ecology}, \bold{20}, 3128--3143.

   Sabeti, P.C. et al. (2007). Genome-wide detection and characterization of positive selection in
   human populations. \emph{Nature}, \bold{449}, 913--918.

   Tang, K. and Thornton, K.R. and Stoneking, M. (2007). A New Approach for Using Genome Scans
   to Detect Recent Positive Selection in the Human Genome. \emph{Plos Biology}, \bold{7}, e171.

}
\seealso{calc_ehh,data2haplohh,scan_hh}
\examples{
#example haplohh object (280 haplotypes, 1424 SNPs)
#see ?haplohh_cgu_bta12 for details
data(haplohh_cgu_bta12)
#computing EHH statistics for the focal SNP at position 456
#which displays a strong signal of selection
res.ehhs<-calc_ehhs(haplohh_cgu_bta12,mrk=456)
}
