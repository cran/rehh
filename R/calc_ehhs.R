#'EHHS and iES computation for a given focal marker
#'@description Compute site-specific Extended Haplotype Homozygosity (EHHS) and integrated EHHS (iES) for a given focal marker.
#'@param haplohh an object of class \code{haplohh} (see \code{\link{data2haplohh}}).
#'@param mrk integer representing the number of the focal marker within the haplohh object
#'or string representing its ID/name.
#'@param limhaplo if there are less than \code{limhaplo} chromosomes that can be used for
#'the calculation of EHH, the calculation is stopped. The option is intended for the case of missing data,
#'which leads to the successive exclusion of haplotypes: the further away from the focal marker
#'the less haplotypes contribute to EHH.
#'@param limhomohaplo if there are less than \code{limhomohaplo} homozygous chromosomes, the
#'calculation is stopped. This option is intended for unphased data and should be invoked only
#'if relatively low frequency variants are not filtered subsequently (see main vignette and Klassmann et al. 2020). 
#'@param limehhs limit at which EHHS stops to be evaluated.
#'@param include_zero_values logical. If \code{FALSE}, return values only for those positions where the calculation is
#'actually performed, i.e. until stopped by reaching either \code{limehh} or \code{limhaplo}. If \code{TRUE}, report EHH values for
#'all markers, the additional ones being zero.
#'@param include_nhaplo logical. If \code{TRUE}, report the number of evaluated haplotypes at each marker
#'(only informative, if missing data leads to a decrease of evaluated haplotypes).
#'@param phased logical. If \code{TRUE} (default) chromosomes are expected to be phased. If \code{FALSE}, the haplotype data is assumed to
#'consist of pairwise ordered chromosomes belonging to diploid individuals.
#'EHHS is then estimated over individuals which are homozygous at the focal marker.
#'@param scalegap scale or cap gaps larger than the specified size to the specified size (default=\code{NA}, i.e. no scaling).
#'@param maxgap maximum allowed gap in bp between two markers. If exceeded, further calculation of EHHS is stopped at the gap
#'(default=\code{NA}, i.e no limitation).
#'@param discard_integration_at_border logical. If \code{TRUE} (default) and computation reaches first or last marker or a gap larger than \code{maxgap},
#'iHH is set to \code{NA}.
#'@param lower_y_bound lower y boundary of the area to be integrated over (default: \code{limehhs}). Can be set
#'to zero for compatibility with the program hapbin.
#'@param interpolate logical. Affects only IES and INES values. If \code{TRUE} (default), integration
#'is performed over a continuous EHHS curve (values are interpolated linearly between consecutive markers),
#'otherwise the EHHS curve decreases stepwise at markers.
#'@details Values for site-specific Extended Haplotype Homozygosity (EHHS) are computed at each position upstream and downstream
#'of the focal marker. These values are integrated with respect to their
#'genomic position to yield an 'integrated EHHS' (iES) value.
#'@return The returned value is a list containing the following elements:
#'\describe{
#'\item{mrk.name}{The name/identifier of the focal marker.}
#'\item{ehhs}{A table containing EHHS values as used by Sabeti et al. (2007),
#'resp. the same values normalized to 1 at the focal marker (nEHHS) as used by Tang et al. (2007).}
#'\item{IES}{Integrated EHHS.}
#'\item{INES}{Integrated  normalized EHHS.}
#'}
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
#'@seealso \code{\link{data2haplohh}}, \code{\link{plot.ehhs}}, \code{\link{calc_ehh}}, \code{\link{scan_hh}}.
#'@examples
#'#example haplohh object (280 haplotypes, 1424 SNPs)
#'#see ?haplohh_cgu_bta12 for details
#'data(haplohh_cgu_bta12)
#'#computing EHHS statistics for the marker "F1205400"
#'#which displays a strong signal of selection
#'ehhs <- calc_ehhs(haplohh_cgu_bta12, mrk = "F1205400")
#'@export
calc_ehhs <-
  function(haplohh,
           mrk,
           limhaplo = 2,
           limhomohaplo = 2,
           limehhs = 0.05,
           include_zero_values = FALSE,
           include_nhaplo = FALSE,
           phased = TRUE,
           scalegap = NA,
           maxgap = NA,
           discard_integration_at_border = TRUE,
           lower_y_bound = limehhs,
           interpolate = TRUE) {
    ##check parameters
    if (!(is.haplohh(haplohh))) {
      stop("Data is not a valid haplohh object.", call. = FALSE)
    }
    
    if (is.numeric(mrk)) {
      mrk <- as.integer(mrk)
      if (mrk < 1) {
        stop(paste0("No marker numbers smaller than 1 allowed."), call. = FALSE)
      }
      if (mrk > nmrk(haplohh)) {
        stop(
          paste0(
            "The marker number ",
            mrk,
            " is bigger than the number of markers in the data set (",
            nmrk(haplohh),
            ")"
          ),
          call. = FALSE
        )
      }
    } else{
      mrk <- as.character(mrk)
      if (!(mrk %in% mrk.names(haplohh))) {
        stop(paste0("Marker '", mrk, "' not found."), call. = FALSE)
      }
      mrk <- which(mrk.names(haplohh) == mrk)
    }
    
    if (limhaplo < 2) {
      stop("limhaplo must be larger than 1.", call. = FALSE)
    }
    if (limhomohaplo < 2) {
      stop("limhomohaplo must be larger than 1.", call. = FALSE)
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
    
    ##perfrom calculations
    res.list <- .Call(
      "CALL_EHHS",
      haplo(haplohh),
      nhap(haplohh),
      nmrk(haplohh),
      mrk,
      limhaplo,
      limhomohaplo,
      limehhs,
      phased
    )
    
    ehhs <-
      data.frame(positions(haplohh), res.list[[2]], res.list[[3]])
    
    colnames(ehhs) <- c("POSITION",
                        "EHHS",
                        "NEHHS")
    
    if (include_nhaplo) {
      ehhs$NHAPLO <- res.list[[1]]
    }
    
    row.names(ehhs) <- mrk.names(haplohh)
    
    ines <- .Call(
      "CALL_INTEGRAL",
      positions(haplohh),
      ehhs$NEHHS,
      mrk,
      limehhs,
      scalegap,
      maxgap,
      discard_integration_at_border,
      lower_y_bound,
      interpolate
    )
    
    ies <- .Call(
      "CALL_INTEGRAL",
      positions(haplohh),
      ehhs$EHHS,
      mrk,
      limehhs,
      scalegap,
      maxgap,
      discard_integration_at_border,
      lower_y_bound,
      interpolate
    )
    
    if (!include_zero_values) {
      nonzeros <- which(rowSums(ehhs[2:3]) > 0)
      
      if (length(nonzeros > 0)) {
        first <- min(nonzeros)
        last <- max(nonzeros)
        
        #if limehh is zero report only the positions next to non-zero positions
        if (limehhs == 0) {
          first <- max(first - 1, 1)
          last <- min(last + 1, nmrk(haplohh))
        }
        
        ehhs <- ehhs[first:last, ]
      } else{
        ehhs <- ehhs[mrk, ]
      }
    }
    
    ##output
    
    l <- list(
      mrk.name = ifelse(
        is.null(mrk.names(haplohh)),
        as.character(mrk),
        mrk.names(haplohh)[mrk]
      ),
      ehhs = ehhs,
      IES = ies,
      INES = ines
    )
    
    return(new("ehhs", l))
  }

#'@rdname calc_ehhs
#'@aliases calc_ehhs
#'@importFrom methods setClass
ehhs <- setClass("ehhs",
                 contains = "list")

#'@importFrom methods is validObject
is.ehhs <- function(x) {
  res <- (is(x, "ehhs") & validObject(x))
  return(res)
}
