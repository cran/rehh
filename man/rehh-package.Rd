\name{rehh-package}
\alias{rehh-package}
\alias{rehh}
\docType{package}
\title{
Searching for footprints of selection using "Extended Haplotype Homozygosity" based statistics
}
\description{
Population genetic data in form of "Single Nucleotide Polymorphisms" (SNPs) is often used to identify genomic regions that have been under recent natural or artificial selection and might provide clues about the molecular mechanisms of adaptation.
The concept of an "Extended Haplotype Homozygosity" (EHH), introduced by (Sabeti 2002), has given rise to several derived statistics designed for whole genome scans. The package provides functions to compute three of these, namely: "iHS" (Voight 2006) for detecting selection within a single population as well as "Rsb" (Tang 2007) and "XP-EHH" (Sabeti 2007) to detect (differential) selection between two populations. Various plotting functions are also included to facilitate visualization and interpretation of these statistics.
}
\details{
\tabular{ll}{
Package: \tab rehh\cr
Version: \tab 2.0.4\cr
License: \tab GPL(>=2)\cr
Depends: \tab gplots , methods\cr
}

Index:
\preformatted{
bifurcation.diagram     Plot of an haplotype bifurcation
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
   Gautier M., Klassmann A., and Vitalis R. (2017). rehh 2.0: a reimplementation of the R package rehh to detect positive selection from haplotype structure.
   \emph{Molecular Ecology Resources}, \bold{17}, 78--90.

   Gautier M. and Vitalis R. (2012). rehh: An R package to detect footprints of
   selection in genome-wide SNP data from haplotype structure. \emph{Bioinformatics}, \bold{28}(8), 1176--1177.

   Gautier M. and Naves M. (2011). Footprints of selection in the ancestral admixture of a New World Cattle breed. \emph{Molecular Ecology}, (20), 3128--3143.

   Sabeti P. C., Reich D. E., Higgins J. M., Levine H. Z. P., Richter D. J., others, (2002). Detecting recent positive selection in the human genome from haplotype structure. \emph{Nature} \bold{419}, 832-837.

   Sabeti P. C., Varilly P., Fry B., Lohmueller J., Hostetter E., others, (2007). Genome-wide detection
and characterization of positive selection in human populations. \emph{Nature} \bold{449}, 913-918.

Tang K., Thornton K. R., Stoneking M., (2007). A new approach for using genome scans to detect recent positive selection in the human genome. \emph{PLoS Biol} \bold{5}, e171.

Voight B. F., Kudaravalli S., Wen X., Pritchard J. K., (2006). A map of recent positive selection in the
human genome. \emph{PLoS Biol} \bold{4}, e72.
}
\keyword{ package }
