\name{rehh-package}
\alias{rehh-package}
\alias{rehh}
\docType{package}
\title{
Searching for footprints of selection using Haplotype Homozygosity based tests
}
\description{
Functions for the detection of footprints of selection on dense SNP data using Extended Homozygosity Haplotype (EHH) based tests. The package includes computation of EHH, iHS (within population) and Rsb (across pairs of populations) statistics. Various plotting functions are also included to facilate vizualisation and interpretation of the results.
}
\details{
\tabular{ll}{
Package: \tab rehh\cr
Version: \tab 1.0\cr
License: \tab GPL(>=2)\cr
Depends: \tab gplots , methods\cr
}

Index:
\preformatted{
bifurcation.diagram     plot of an haplotype bifurcation
                        diagram
calc_ehh                EHH and iHH computations at a 
                        given core SNP
calc_ehhs               EHHS and iES computations at a 
                        given core SNP
data2haplohh            Converting data into an object 
                        of class haplohh
distribplot             Distribution of standardized 
                        iHS or Rsb values
haplohh_cgu_bta12       Example of an 'haplohh' object
ies2rsb                 Compute Rsb (standardized ratio of 
                        iES from two populations)
ihh2ihs                 Compute iHS (standardized iHH)
ihsplot                 Plot iHS over a genome
make.example.files      Generate example input files
rsbplot                 Plot Rsb over a genome
scan_hh                 Computing EHH based statistics over 
                        a whole chromosome
wgscan.cgu              Whole genome scan results for the CGU
                        (Creole from Guadeloupe island)
wgscan.eut              Whole genome scan results for a pool of
                        European taurine cattle
}
}
\references{
   Gautier, M. and Vitalis, R. (2012). rehh: An R package to detect footprints of 
   selection in genome-wide SNP data from haplotype structure. \emph{Bioinformatics}, accepted.
% \emph{Bioinformatics}, \bold{20}, 3128--3143.
}
\keyword{ package }
