\name{haplohh-class}
\Rdversion{1.1}
\docType{class}
\alias{haplohh-class}

\title{Class \code{"haplohh"}}
\description{
An object of class haplohh contains all relevant haplotype information (see below).
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("haplohh", ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{haplo}}{Object of class \code{"matrix"}: haplotypes with alleles coded as 0 (missing data), 1 (ancestral allele) or 2 (derived allele)}
    \item{\code{position}}{Object of class \code{"numeric"}:position of the SNPs in the chromosome}
    \item{\code{snp.name}}{Object of class \code{"character"}:names of the SNP}
    \item{\code{chr.name}}{Object of class \code{"numeric"}:name of the chromosome SNPs are mapping to}
    \item{\code{nhap}}{Object of class \code{"numeric"}:number of haplotypes}
    \item{\code{nsnp}}{Object of class \code{"numeric"}:number of SNPs in the haplotypes}
  }
}
\seealso{data2haplohh}
\examples{
showClass("haplohh")
}
\keyword{classes}
