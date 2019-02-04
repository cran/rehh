\name{bifurcation.diagram}
\alias{bifurcation.diagram}
\title{plot of an haplotype bifurcation diagram}
\description{A haplotype bifurcation diagram visualizes the
decay of EHH of a "core" allele of a focal SNP at increasing distances.
}
\usage{
bifurcation.diagram(haplohh, mrk_foc, all_foc=1, nmrk_l=10, nmrk_r=10,
                limhapcount = 10, refsize = 0.1, linecol = "blue",
                main_leg = NA, xlab_leg = "Position")
}
\arguments{
  \item{haplohh}{An object of class \code{haplohh} (see \code{data2haplohh}).}
  \item{mrk_foc}{Either the number of the marker in the haplohh object (as integer) or its name (as string) to specify the focal marker}
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
The function `bifurcation.diagram()` draws haplotype bifurcation diagrams [@Sabeti2002] that visualize the decay of \emph{EHH} around a focal SNP. A stark contrast of ancestral and derived bifurcation diagrams should correspond to outlier values of \emph{ihs}. In the plot the root (focal SNP) is identified by a vertical dashed line.
The diagram is bi-directional, portraying decay along both sides of the focal SNP. Moving in one direction, each marker is an opportunity for a bifurcation further differentiating (extended) haplotypes. The thickness of the lines corresponds to the number of chromosomes with the same extended haplotype.
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
