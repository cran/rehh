#'Compute iHH, iES and inES over a whole chromosome without cut-offs
#'@description Compute integrated EHH (iHH), integrated EHHS (iES) and integrated normalized EHHS (inES) for all markers of a chromosome (or linkage group).
#'This function computes the statistics by a slightly different algorithm than \code{\link{scan_hh}}: it sidesteps the calculation of EHH and EHHS values and their subsequent integration and
#'consequently no cut-offs relying on these values can be specified. Instead,
#'it computes the (full) lengths of pairwise shared haplotypes and averages them afterwords.
#'
#'This function is primarily intended for the study of general properties of these statistics using simulated data.
#'@param haplohh an object of class \code{haplohh} (see \code{\link{data2haplohh}})
#'@param phased logical. If \code{TRUE} (default) chromosomes are expected to be phased. If \code{FALSE}, the haplotype data is assumed to
#'consist of pairwise ordered chromosomes belonging to diploid individuals.
#'EHH(S) is then estimated over individuals which are homozygous at the focal marker.
#'@param polarized logical. \code{TRUE} by default. If \code{FALSE}, use major and minor allele instead of ancestral and derived. If there
#'are more than two alleles then the minor allele refers to the second-most frequent allele.
#'@param maxgap maximum allowed gap in bp between two markers. If exceeded, further calculation of EHH(S) is stopped at the gap
#'(default=\code{NA}, i.e. no limitation).
#'@param max_extend maximum distance in bp to extend shared haplotypes away from the focal marker.
#'(default \code{NA}, i.e. no limitation).
#'@param discard_integration_at_border logical. If \code{TRUE} (default) and computation of any of the statistics reaches first or last
#'marker or a gap larger than \code{maxgap}, iHH, iES and inES are set to \code{NA}.
#'@param geometric.mean logical. If \code{FALSE} (default), the standard arithmetic mean is used to average
#'shared haplotype lengths. If \code{TRUE}
#'the geometric mean is used instead (IES values are undefined in this case). Note that usage of the geometric mean has not
#'yet been studied formally and should be considered experimental!
#'@param threads number of threads to parallelize computation
#'
#'@details Integrated EHH (iHH), integrated EHHS (iES) and integrated normalized EHHS (inES)
#'are computed for all markers of the chromosome (or linkage group). This function sidesteps
#'the computation of EHH and EHHS values and their stepwise integration. Instead, the length of all shared haplotypes
#'is computed and afterwords averaged.  In the absence of missing values the
#'statistics are identical to those calculated by \code{\link{scan_hh}} with settings
#'\code{limehh = 0}, \code{limehhs = 0}, \code{lower_ehh_y_bound = 0} and \code{interpolate = FALSE}, yet this function is faster.
#'
#'Application of a cut-off is necessary for reducing the spurious signals
#'of selection caused by single shared haplotypes of extreme length. Hence, e.g. for human experimental data
#'it might be reasonable to set \code{max_extend} to 1 or 2 Mb.
#'
#'\code{\link{scan_hh}} computes the statistics iHH_A, ihh_D and iES/inES separately,
#'while this function calculates them simultaneously. Hence, 
#'if \code{discard_integration_at_border} is set to \code{TRUE} and the extension of shared haplotypes
#'reaches a border (i.e. chromosomal boundaries or a gap larger than \code{maxgap}), this function discards
#'all statistics.
#'
#'The handling of missing values is different, too: \code{\link{scan_hh}} "removes" chromosomes with missing values from further calculations.
#'EHH and EHHS are then calculated for the remaining chromosomes which can accidentally yield an increase in EHH or EHHS.
#'This can not happen with \code{scan_hh_full()} which treats each missing value of a marker
#'as if it were a new allele - terminating any shared haplotype, but does changing the
#'set of considered chromosomes. Thus, missing values
#'cause a faster decay of EHH(S) with function \code{scan_hh_full()}.
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
#'Klassmann, A. and Gautier, M. (2020). Detecting selection using Extended Haplotype
#'Homozygosity-based statistics on unphased or unpolarized data (preprint). 
#'https://doi.org/10.22541/au.160405572.29972398/v1
#'
#'Sabeti, P.C. et al. (2002). Detecting recent positive selection in the human genome from haplotype structure. \emph{Nature}, \strong{419}, 832-837.
#'
#'Sabeti, P.C. et al. (2007). Genome-wide detection and characterization of positive selection in human populations. \emph{Nature}, \strong{449}, 913-918.
#'
#'Tang, K. and Thornton, K.R. and Stoneking, M. (2007). A New Approach for Using Genome Scans to Detect Recent Positive Selection in the Human Genome. \emph{Plos Biology}, \strong{7}, e171.
#'
#'Voight, B.F. and Kudaravalli, S. and Wen, X. and Pritchard, J.K. (2006). A map of recent positive selection in the human genome. \emph{Plos Biology}, \strong{4}, e72.
#'@seealso \code{\link{data2haplohh}}, \code{\link{scan_hh}},
#'\code{\link{ihh2ihs}}, \code{\link{ines2rsb}}, \code{\link{ies2xpehh}}
#'@examples
#'#example haplohh object (280 haplotypes, 1424 SNPs)
#'#see ?haplohh_cgu_bta12 for details
#'data(haplohh_cgu_bta12)
#'#using function scan_hh() with no cut-offs
#'scan <- scan_hh(haplohh_cgu_bta12, discard_integration_at_border = FALSE,
#'limehh = 0, limehhs = 0, lower_ehh_y_bound = 0, interpolate = FALSE)
#'#using function scan_hh_full()
#'scan_full <- scan_hh_full(haplohh_cgu_bta12, discard_integration_at_border = FALSE)
#'#both yield identical results within numerical precision
#'all.equal(scan, scan_full)
#'@export
scan_hh_full <-
  function(haplohh,
           phased = TRUE,
           polarized = TRUE,
           maxgap = NA,
           max_extend = NA,
           discard_integration_at_border = TRUE,
           geometric.mean = FALSE,
           threads = 1) {
    ##check parameters
    if (!(is.haplohh(haplohh))) {
      stop("Data is not a valid haplohh object.", call. = FALSE)
    }
    if (!is.na(maxgap) & (!is.numeric(maxgap) | maxgap < 1)) {
      stop("maxgap must be a positive integer number.", call. = FALSE)
    }
    if (!is.na(max_extend) &
        (!is.numeric(max_extend) | max_extend < 1)) {
      stop("max_extend must be a positive integer number.", call. = FALSE)
    }
    
    if (is.na(maxgap)) {
      maxgap <- diff(range(positions(haplohh))) + 1
    }
    if (is.na(max_extend)) {
      max_extend <- diff(range(positions(haplohh))) + 1
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
      "CALL_SCAN_HH_FULL",
      haplo(haplohh),
      nhap(haplohh),
      nmrk(haplohh),
      #ancestral or major allele
      allele_stats[1, ],
      #derived or minor allele
      allele_stats[2, ],
      positions(haplohh),
      as.integer(maxgap),
      as.integer(max_extend),
      as.integer(phased),
      as.integer(discard_integration_at_border),
      as.integer(geometric.mean),
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
      colnames(res)[3:10] <-
        c(
          "FREQ_MAJ",
          "FREQ_MIN",
          "NHAPLO_MAJ",
          "NHAPLO_MIN",
          "IHH_MAJ",
          "IHH_MIN",
          "IES",
          "INES"
        )
    }
    
    return(res)
  }
