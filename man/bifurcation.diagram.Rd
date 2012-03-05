\name{bifurcation.diagram}
\alias{bifurcation.diagram}
\title{plot of an haplotype bifurcation diagram}
\description{Haplotype Bifurcation diagram visualizes the
breakdown of LD at increasing distances from the core allele at the selected focal SNPs.
}
\usage{
bifurcation.diagram(haplohh,mrk_foc,all_foc=1,nmrk_l=10,nmrk_r=10,
limhapcount=10,refsize=0.1,linecol="blue",
main_leg=NA,xlab_leg="Position")
}
\arguments{
  \item{haplohh}{An object of class \code{haplohh} (see \code{data2haplohh}).}
  \item{mrk_foc}{Integer representing the number of the focal marker}
  \item{all_foc}{either 1 or 2 depending on the chosen core allele (resp. ancestral or derived)}
  \item{nmrk_l}{Number of markers to be considered upstream of the focal SNP}
  \item{nmrk_r}{Number of markers to be considered downstream of the focal SNP}
  \item{limhapcount}{Minimal number of haplotypes containing the core allele at the focal SNP}
  \item{refsize}{Controls the relative width of the diagram lines on the plot}
  \item{linecol}{Color of the lines on the diagram}
  \item{main_leg}{Main legend of the diagram. By default, the name of the SNP together with the allele considered}
  \item{xlab_leg}{Legend on the xaxis of the diagram}
}
\details{
Haplotype Bifurcation diagram visualizes the
breakdown of LD at increasing distances from the core allele at the selected focal SNPs. The root (focal SNP) of each diagram is the core allele, 
identified by a vertical dashed line. The
diagram is bi-directional, portraying both centromere-proximal and centromere-distal LD.
Moving in one direction, each marker is an opportunity for a node; the diagram either
divides or not based on whether both or only one allele is present. Thus the breakdown of
LD on the core haplotype background is portrayed at progressively longer distances. The
thickness of the lines corresponds to the number of samples with the indicated long-
distance haplotype.
}
\value{The function returns a plot.}
\references{
   Sabeti, P.C. and Reich, D.E. and Higgins, J.M. and Levine, H.Z.P and Richter, D.J. and 
   Schaffner, S.F. and Gabriel, S.B. and Platko, J.V. and Patterson, N.J. and McDonald, G.J. 
   and Ackerman, H.C. and Campbell, S.J. and Altshuler, D. and Cooper, R. and Kwiatkowski, D.
   and Ward, R. and Lander, E.S. (2002). Detecting recent positive selection in the human genome 
   from haplotype structure. \emph{Nature}, \bold{419}, 832--837.
}
\examples{
#example haplohh object (280 haplotypes, 1424 SNPs)
#see ?haplohh_cgu_bta12 for details
data(haplohh_cgu_bta12) 
#plotting bifurcation diagram for both ancestral and derived allele
#from the focal SNP at position 456
#which display a strong signal of selection
layout(matrix(1:2,2,1))
#ancestral allele
bifurcation.diagram(haplohh_cgu_bta12,mrk_foc=456,all_foc=1,
nmrk_l=20,nmrk_r=20)
#derived allele
bifurcation.diagram(haplohh_cgu_bta12,mrk_foc=456,all_foc=2,
nmrk_l=20,nmrk_r=20)
## 
dev.off()
}
