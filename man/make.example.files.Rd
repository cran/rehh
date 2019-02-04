\name{make.example.files}
\alias{make.example.files}
\title{Creating example input files}
\description{
This function copies the following example files to the working directory:
\itemize{
\item \code{bta12_cgu.hap} an haplotype input file in standard format
\item \code{bta12_cgu.thap} an haplotype input file in transposed format
\item \code{bta12_hapguess_switch.out} an haplotype input file in fastphase output format
\item \code{map.inp} a SNP information input file
}
These files contain data for 280 haplotypes (originating from 140 individuals belonging to the Creole
cattle breed from Guadeloupe) of 1,424 SNPs mapping to bovine chromosome 12 (BTA12) (see reference below).
}
\usage{
make.example.files()
}
%- maybe also 'usage' for other objects documented here.
\references{
   Gautier, M. and Naves, M. (2011). Footprints of selection in the ancestral admixture of a New World
   Creole cattle breed. \emph{Molecular Ecology}, \bold{20}, 3128--3143.
}
\seealso{data2haplohh}
\examples{make.example.files()}
