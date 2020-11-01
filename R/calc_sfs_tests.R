#'Calculate site frequency spectrum test statistics
#'@description Calculate site frequency spectrum (SFS) tests Tajima's D, Fay & Wu's H and Zeng's E.
#'@param haplohh an object of class \code{haplohh} (see \code{\link{data2haplohh}})
#'@param polarized logical. \code{TRUE} by default. If \code{FALSE}, use major and minor allele instead of ancestral and derived. If there
#'are more than two alleles then the minor allele refers to the second-most frequent allele.
#'Note that Tajima's D remains unchanged, Fay & Wu's H is always zero for folded spectra and Zeng's E becomes equal to Tajima's D.
#'@param window_size size of sliding windows. If \code{NA} (default), there will be only
#'one window covering the whole length of the chromosome.
#'@param overlap size of window overlap (default 0, i.e. no overlap).
#'@param right logical, indicating if the windows should be closed on the right and open on the left (default) or vice versa.
#'@param min_n_mrk minimum number of (polymorphic) markers per window.
#'@param verbose logical. \code{TRUE} by default; reports if multi-allelic sites are removed.
#'@details Neutrality tests based on the site frequency spectrum (SFS) are
#'largely unrelated to EHH-based methods. The tests provided here are implemented
#'elsewhere, too (e.g. in package \href{https://cran.r-project.org/package=PopGenome}{PopGenome}).
#'
#'Each test compares two estimations of the \emph{scaled mutation rate} theta,
#'which all have the same expected value under neutrality. Deviations from zero indicate
#'violations of the neutral null model, typically population size changes, population subdivision or selection.
#'Tajima's D and Fay & Wu's H become negative in presence of an almost completed sweep, Zeng's E becomes
#'positive for some time after it. Significance can typically be assigned only by
#'simulations.
#'
#'The standard definition of the tests cannot cope with missing values and typically markers
#'with missing genotypes must be discarded. Ferretti (2012) provides an extension
#'that can handle missing values (without discarding any non-missing values). In this package, 
#'only the first moments (the theta-estimators themselves) are adapted accordingly, 
#'but not the second moments (their variances), because the latter is computationally demanding
#'and the resulting bias relatively small. It is recommended, though, to discard markers or haplotypes 
#'with more than 20\% missing values.
#'
#'Multi-allelic markers are always removed since the tests rely on the "infinite sites model" which
#'implies that all polymorphic markers are bi-allelic. 
#'Monomorphic markers can be present, but are irrelevant for the tests.
#'
#'@return A data frame with window coordinates, the number of contained (polymorphic) markers, Watterson's, Tajima's and Zeng's
#'estimators of theta and the test statistics of Tajima's D, Fay & Wu's H and Zeng's E.
#'@examples
#'make.example.files()
#'# neutral evolution
#'hh <- data2haplohh("example_neutral.vcf", verbose = FALSE)
#'calc_sfs_tests(hh)
#'# strong selective sweep
#'hh <- data2haplohh("example_sweep.vcf", verbose = FALSE)
#'calc_sfs_tests(hh)
#'remove.example.files()
#'@references Watterson, G.A. (1975). On the number of segregating sites in genetical models without recombination.
#'\emph{Theoretical Population Biology} \strong{7}(2) 256-276.
#'
#'Tajima, F. (1983). Evolutionary relationship of DNA sequences in finite populations.
#'\emph{Genetics} \strong{105}(2) 437-60.
#'
#'Tajima, F. (1989). Statistical method for testing the neutral mutation hypothesis by DNA polymorphism.
#'\emph{Genetics} \strong{123}(3) 585-95.
#'
#'Fay, J. and Wu, C. (2000). Hitchhiking under positive Darwinian selection. \emph{Genetics}
#'\strong{155}(3) 1405-13.
#'
#'Zeng, E. et al. (2006). Statistical tests for detecting positive selection by utilizing high-frequency variants.
#'\emph{Genetics} \strong{174}(3) 1431-9.
#'
#'Ferretti, L. and Raineri, E. and Ramos-Onsins, S. (2012). Neutrality tests for sequences with missing data.
#'\emph{Genetics} \strong{191}(4) 1397-401.
#'@export
calc_sfs_tests <-
  function(haplohh,
           polarized = TRUE,
           window_size = NA,
           overlap = 0,
           right = TRUE,
           min_n_mrk = 1,
           verbose = TRUE) {
    # need integer numbers, otherwise "%%" causes trouble with i386
    overlap <- as.integer(overlap)
    min_n_mrk <- as.integer(min_n_mrk)
    
    if (is.na(window_size)) {
      window_size <- ceiling(max(haplohh@positions))
    }
    window_size <- as.integer(window_size)
    
    if (window_size < 1) {
      stop("Window size has to be a positive integer number.", call. = FALSE)
    }
    if (is.na(overlap) |
        overlap < 0 |
        overlap >= window_size |
        overlap %% 1L != 0L |
        (overlap != 0L & window_size %% overlap != 0L)) {
      stop("'overlap' has to be zero or an integer factor of 'window_size'.",
           call. = FALSE)
    }
    if (min_n_mrk < 1 | min_n_mrk %% 1 != 0L) {
      stop("'min_n_mrk' has to be a positive integer number.",
           call. = FALSE)
    }
    
    # if present, remove multi-allelic sites
    if (max(apply(haplohh@haplo, 2, function(x) {
      length(na.omit(unique(x)))
    })) > 2) {
      haplohh <-
        subset(
          haplohh,
          max_alleles = 2,
          min_perc_geno.mrk = floor(2 / nhap(haplohh) * 100),
          verbose = verbose
        )
    }
    
    step <- ifelse(overlap != 0, overlap, window_size)
    window_left <-
      seq(
        floor(min(positions(haplohh) - right * 1) / window_size) * window_size,
        ceiling(max(positions(haplohh))) - window_size + 1 - right * 1,
        step
      )
    window_right <- window_left + window_size
    windows <- cbind(window_left, window_right)
    
    n_mrk <- vector(mode = "integer", length = nrow(windows))
    results <- matrix(0, nrow = nrow(windows), ncol = 6)
    
    .Call(
      "CALL_SFS_TESTS",
      haplo(haplohh),
      nhap(haplohh),
      nmrk(haplohh),
      positions(haplohh),
      polarized,
      windows,
      nrow(windows),
      right,
      min_n_mrk,
      n_mrk,
      results
    )
    
    df <- data.frame(chr.name(haplohh), windows, n_mrk, results)
    colnames(df) <-
      c(
        "CHR",
        "START",
        "END",
        "N_MRK",
        "THETA_S",
        "THETA_PI",
        "THETA_L",
        "TAJIMA_D",
        "FAY_WU_H",
        "ZENG_E"
      )
    
    return(df)
  }
