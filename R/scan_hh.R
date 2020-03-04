#'Compute iHH, iES and inES over a whole chromosome
#'@description Compute integrated EHH (iHH), integrated EHHS (iES) and integrated normalized EHHS (inES)
#'for all markers of a chromosome (or linkage group).
#'
#'@param haplohh an object of class \code{haplohh} (see \code{\link{data2haplohh}})
#'@param limhaplo if there are less than \code{limhaplo} chromosomes that can be used for
#'the calculation of EHH(S), the calculation is stopped. The option is intended for the case of missing data,
#'which leads to the successive exclusion of haplotypes: the further away from the focal marker
#'the less haplotypes contribute to EHH(S).
#'@param limhomohaplo if there are less than \code{limhomohaplo} homozygous chromosomes, the
#'calculation is stopped. This option is intended for unphased data and should be invoked only
#'if relatively low frequency variants are not filtered subsequently (see main vignette and Klassmann et al. 2020). 
#'@param limehh limit at which EHH stops to be evaluated.
#'@param limehhs limit at which EHHS stops to be evaluated.
#'@param phased logical. If \code{TRUE} (default) chromosomes are expected to be phased. If \code{FALSE}, the haplotype data is assumed to
#'consist of pairwise ordered chromosomes belonging to diploid individuals.
#'EHH(S) is then estimated over individuals which are homozygous at the focal marker.
#'@param polarized logical. \code{TRUE} by default. If \code{FALSE}, use major and minor allele instead of ancestral and derived. If there
#'are more than two alleles then the minor allele refers to the second-most frequent allele.
#'@param scalegap scale or cap gaps larger than the specified size to the specified size (default=\code{NA}, i.e. no scaling).
#'@param maxgap maximum allowed gap in bp between two markers. If exceeded, further calculation of EHH(S) is stopped at the gap
#'(default=\code{NA}, i.e no limitation).
#'@param discard_integration_at_border logical. If \code{TRUE} (default) and computation reaches first or last marker or a gap larger than \code{maxgap},
#'iHH, iES and inES are set to \code{NA}.
#'@param lower_ehh_y_bound lower y boundary of the area to be integrated over (default: \code{limehh}). Can be set
#'to zero for compatibility with the program hapbin.
#'@param lower_ehhs_y_bound lower y boundary of the area to be integrated (default: \code{limehhs}). Can be set
#'to zero for compatibility with the program hapbin.
#'@param threads number of threads to parallelize computation
#'
#'@details Integrated EHH (iHH), integrated EHHS (iES) and integrated normalized EHHS (inES)
#'are computed for all markers of the chromosome (or linkage group). This function is several
#'times faster as a procedure calling in turn \code{calc_ehh} and \code{calc_ehhs}
#'for all markers. To perform a whole genome-scan this function needs
#'to be called for each chromosome and results concatenated.
#'@return The returned value is a dataframe with markers in rows and the following columns
#'\enumerate{
#'\item chromosome name
#'\item position in the chromosome
#'\item sample frequency of the ancestral / major allele
#'\item sample frequency of the second-most frequent remaining allele
#'\item number of evaluated haplotypes at the focal marker for the ancestral / major allele
#'\item number of evaluated haplotypes at the focal marker for the second-most frequent remaining allele  
#'\item iHH of the ancestral / major allele
#'\item iHH of the second-most frequent remaining allele
#'\item iES (used by Sabeti et al 2007)
#'\item inES (used by Tang et al 2007)}
#'Note that in case of unphased data the evaluation is restricted to
#'haplotypes of homozygous individuals which reduces the power
#'to detect selection, particularly for iHS (for appropriate parameter setting 
#'see the main vignette and Klassmann et al (2020)).
#'
#
#'@references Gautier, M. and Naves, M. (2011). Footprints of selection in the ancestral admixture of a New World Creole cattle breed. \emph{Molecular Ecology}, \strong{20}, 3128-3143.
#'
#'Klassmann, A. et al. (2020). Detecting selection using Extended Haplotype
#'Homozygosity (EHH)-based statistics on unphased or unpolarized data. (submitted).
#'
#'Sabeti, P.C. et al. (2002). Detecting recent positive selection in the human genome from haplotype structure. \emph{Nature}, \strong{419}, 832-837.
#'
#'Sabeti, P.C. et al. (2007). Genome-wide detection and characterization of positive selection in human populations. \emph{Nature}, \strong{449}, 913-918.
#'
#'Tang, K. and Thornton, K.R. and Stoneking, M. (2007). A New Approach for Using Genome Scans to Detect Recent Positive Selection in the Human Genome. \emph{Plos Biology}, \strong{7}, e171.
#'
#'Voight, B.F. and Kudaravalli, S. and Wen, X. and Pritchard, J.K. (2006). A map of recent positive selection in the human genome. \emph{Plos Biology}, \strong{4}, e72.
#'@seealso \code{\link{data2haplohh}}, \code{\link{calc_ehh}}, \code{\link{calc_ehhs}}
#'\code{\link{ihh2ihs}},\code{\link{ines2rsb}}, \code{\link{ies2xpehh}}
#'@examples
#'#example haplohh object (280 haplotypes, 1424 SNPs)
#'#see ?haplohh_cgu_bta12 for details
#'data(haplohh_cgu_bta12)
#'scan <- scan_hh(haplohh_cgu_bta12)
#'@export
scan_hh <-
  function(haplohh,
           limhaplo = 2,
           limhomohaplo = 2,
           limehh = 0.05,
           limehhs = 0.05,
           phased = TRUE,
           polarized = TRUE,
           scalegap = NA,
           maxgap = NA,
           discard_integration_at_border = TRUE,
           lower_ehh_y_bound = limehh,
           lower_ehhs_y_bound = limehhs,
           threads = 1) {
    ##check parameters
    if (!(is.haplohh(haplohh))) {
      stop("Data is not a valid haplohh object.", call. = FALSE)
    }
    if (limhaplo < 2) {
      stop("limhaplo must be larger than 1.", call. = FALSE)
    }
    if (limhomohaplo < 2) {
      stop("limhomohaplo must be larger than 1.", call. = FALSE)
    }
    if (limehh < 0 |
        limehh > 1) {
      stop("limehh must lie between 0 and 1.", call. = FALSE)
    }
    if (limehhs < 0 |
        limehhs > 1) {
      stop("limehhs must lie between 0 and 1.", call. = FALSE)
    }
    if (is.na(maxgap)) {
      maxgap <- (max(positions(haplohh)) + 1)
    }
    
    if (is.na(scalegap)) {
      scalegap <- (max(positions(haplohh)) + 1)
    } else if (scalegap > maxgap) {
      stop("scalegap has to be smaller than maxgap in order to have an effect.",
           call. = FALSE)
    }
    
    ##perform calculation
    
    #calculate for each column (=marker) the alleles to be used.
    #if polarized = TRUE the first allele is the ancestral allele
    #and the second allele the derived allele with highest frequency
    #if polarized = FALSE the first allele is the allele with highest frequency
    #and the second allele the allele with second-highest frequency
    #Additionally, their respective frequencies and the number of all alleles
    #is calculated (= number of chromosomes, if no missing values).
    #Use integers to save memory space.
    allele_stats <- apply(haplo(haplohh), 2, function(x) {
      t <- tabulate(x + 1L)
      #only one allele: return it as first and "-1" for the second
      if (length(t) == 1) {
        return(c(0L, NA, t[1], 0L, t[1]))
      } else{
        if (polarized) {
          #first allele is '0' (ancestral), second allele is most frequent derived
          alleles <- c(0L, which.max(t[-1]))
        } else{
          #first allele is major, second allele is minor (second most frequent)
          alleles <- order(t, decreasing = TRUE)[1:2] - 1L
        }
        return(c(alleles,
                 t[alleles + 1L],
                 sum(t)))
      }
    })
    
    res.list <- .Call(
      "CALL_SCAN_HH",
      haplo(haplohh),
      nhap(haplohh),
      nmrk(haplohh),
      #ancestral or major allele
      allele_stats[1, ],
      #derived or minor allele
      allele_stats[2, ],
      as.integer(limhaplo),
      as.integer(limhomohaplo),
      as.double(limehh),
      as.double(limehhs),
      as.integer(scalegap),
      as.integer(maxgap),
      positions(haplohh),
      as.integer(phased),
      as.integer(discard_integration_at_border),
      as.double(lower_ehh_y_bound),
      as.double(lower_ehhs_y_bound),
      as.integer(threads)
    )
    
    #there is only one chromosome in a haplohh
    CHR <- rep(chr.name(haplohh), nmrk(haplohh))
    #position of marker within chromosome
    POSITION <- positions(haplohh)
    ##output
    
    res <-
      data.frame(
        CHR,
        POSITION,
        FREQ_A = allele_stats[3, ] / allele_stats[5, ],
        FREQ_D = allele_stats[4, ] / allele_stats[5, ],
        NHAPLO_A = res.list[[1]],
        NHAPLO_D = res.list[[2]],
        IHH_A = res.list[[3]],
        IHH_D = res.list[[4]],
        IES = res.list[[5]],
        INES = res.list[[6]],
        row.names = mrk.names(haplohh)
      )
    
    #if not polarized, change column names
    if (!polarized) {
      colnames(res)[3:8] <-
        c("FREQ_MAJ",
          "FREQ_MIN",
          "NHAPLO_MAJ",
          "NHAPLO_MIN",
          "IHH_MAJ",
          "IHH_MIN")
    }
    
    return(res)
  }
