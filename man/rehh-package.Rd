\name{rehh-package}
\alias{rehh-package}
\alias{rehh}
\docType{package}
\title{
Searching for footprints of selection using Haplotype Homozygosity based tests
}
\description{
Functions for the detection of footprints of selection on dense SNP data using Extended Homozygosity Haplotype (EHH) based tests. The package includes computation of EHH, iHS (within population) and Rsb and XP-EHH (across pairs of populations) statistics. Various plotting functions are also included to facilitate visualization and interpretation of the results.
}
\details{
\tabular{ll}{
Package: \tab rehh\cr
Version: \tab 2.0.2\cr
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
                        iHS, Rsb or XP-EHH values
haplohh_cgu_bta12       Example of an 'haplohh' object
ies2rsb                 Compute Rsb (standardized ratio of 
                        iES from two populations)
ies2xpehh               Compute XP-EHH (standardized ratio of 
                        iES from two populations)
ihh2ihs                 Compute iHS (standardized iHH)
ihsplot                 Plot iHS over a genome
make.example.files      Generate example input files
rsbplot                 Plot Rsb over a genome
xpehhplot               Plot XP-EHH over a genome
scan_hh                 Computing EHH based statistics over 
                        a whole chromosome
}
}
\references{
   Gautier M., Klassmann A., and Vitalis R. (2016). rehh: An R package to detect footprints of selection in genome-wide SNP data from haplotype structure.
   \emph{Molecular Ecology Resources}, accepted

   Gautier M. and Vitalis R. (2012). rehh: An R package to detect footprints of 
   selection in genome-wide SNP data from haplotype structure. \emph{Bioinformatics}, \bold{28}(8), 1176--1177.
   
   Gautier M. and Naves M. (2011). Footprints of selection in the ancestral admixture of a New World Cattle breed. \emph{Molecular Ecology}, (20), 3128--3143.
}
\keyword{ package }
