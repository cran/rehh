#'Translate object of \code{\link{haplohh-class}} into SweepFinder format
#'@description Extract allele frequencies of an object of class \code{\link{haplohh-class}}
#'and returns a table in SweepFinder input format.
#'@param haplohh object of class \code{\link{haplohh-class}}.
#'@param polarized logical. If \code{TRUE} (default), flag "folded" is set to 0, otherwise to 1.
#'@param verbose logical. If \code{TRUE} (default), prints filter statements.
#'@details SweepFinder and SweeD are two stand-alone programs which
#'implement the same method to detect selective sweeps using the
#'allele frequency at each site. This function calculates these frequencies
#'from a \code{\link{haplohh-class}} and returns a table which
#'can be saved into a file (with tabs as separators, without row names and quotes) that can
#'be used as input for the two programs.
#'
#'Sites with less than two haplotypes genotyped or with more than two alleles are removed.
#'If \code{polarized}, sites monomorphic for the ancestral allele are removed, too.
#'
#'@return A dataframe with four columns:
#'\itemize{
#'\item \strong{position} marker position
#'\item \strong{x} (absolute) frequency of the alternative (derived) variant
#'\item \strong{n} number of non-missing genotypes
#'\item \strong{folded} a flag marking polarization
#'}
#'@seealso \code{\link{haplohh-class}}, \code{\link{data2haplohh}}
#'@references DeGiorgio, M., and, Huber, CD and Hubisz, MJ and, Hellmann, I. and Nielsen, R. (2016)
#'SweepFinder2: increased robustness and flexibility. \emph{Bioinformatics} \strong{32}:1895-1897
#'
#'Pavlidis, P., D. Zivkovic, A. Stamatakis, and N. Alachiotis, (2013)
#'SweeD: likelihood-based detection of selective sweeps in thousands of genomes.
#'\emph{Molecular Biology and Evolution} \strong{30}: 2224-34.
#'@examples #example
#'# sweepfinder example from vignette
#'make.example.files()
#'hh <- data2haplohh("example_sweep_with_recombination.vcf")
#'haplohh2sweepfinder(hh)
#'remove.example.files()
#'@export
haplohh2sweepfinder <-
  function(haplohh,
           polarized = TRUE,
           verbose = TRUE) {
    # remove multi-allelic sites (sites with more than two alleles)
    # and sites with less than two sequences genotyped
    haplohh <-
      subset(
        haplohh,
        max_alleles = 2,
        min_perc_geno.mrk = floor(2 / nhap(haplohh) * 100),
        verbose = verbose
      )
    
    n <- colSums(!is.na(haplohh@haplo))
    x <- colSums(!is.na(haplohh@haplo) & haplohh@haplo != 0L)
    df <- data.frame(
      position = haplohh@positions,
      x = x,
      n = n,
      folded = (!polarized) * 1L
    )
    # sites that are monomorphic for the ancestral allele are not allowed
    if (polarized) {
      df <- df[df$x > 0, ]
    }
    return(df)
  }
