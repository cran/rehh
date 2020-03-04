#'EHH and iHH computation for a given focal marker
#'@description Compute Extended Haplotype Homozygosity (EHH) and integrated EHH (iHH) for a given focal marker.
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
#'@param limehh limit at which EHH stops to be evaluated
#'@param include_zero_values logical. If \code{FALSE}, return values only for those positions where the calculation is
#'actually performed, i.e. until stopped by reaching either \code{limehh} or \code{limhaplo}. If \code{TRUE}, report EHH values for
#'all markers, the additional ones being zero.
#'@param include_nhaplo logical. If \code{TRUE}, report the number of evaluated haplotypes at each marker
#'(only informative, if missing data leads to a decrease of evaluated haplotypes).
#'@param phased logical. If \code{TRUE} (default) chromosomes are expected to be phased. If \code{FALSE}, the haplotype data is assumed to
#'consist of pairwise ordered chromosomes belonging to diploid individuals.
#'EHH is then estimated over individuals which are homozygous at the focal marker.
#'@param polarized logical. \code{TRUE} by default. If \code{FALSE}, use major and minor allele instead of ancestral and derived. If there
#'are more than two alleles then the minor allele refers to the second-most frequent allele.
#'@param scalegap scale or cap gaps larger than the specified size to the specified size (default=\code{NA}, i.e. no scaling).
#'@param maxgap maximum allowed gap in bp between two markers. If exceeded, further calculation of EHH is stopped at the gap
#'(default=\code{NA}, i.e no limitation).
#'@param discard_integration_at_border logical. If \code{TRUE} (default) and computation reaches first or last marker or a gap larger than \code{maxgap},
#'iHH is set to \code{NA}.
#'@param lower_y_bound lower y boundary of the area to be integrated over (default: \code{limehh}). Can be set
#'to zero for compatibility with the program hapbin.
#'@details Values for allele-specific Extended Haplotype Homozygosity (EHH) are computed
#'upstream and downstream of the focal marker for each of its alleles.
#'These values are integrated with respect to their genomic positions to yield
#'an 'integrated EHH' (iHH) value for each allele.
#'@return The returned value is a list containing the following elements:
#'\describe{
#'\item{mrk.name}{The name/identifier of the focal marker.}
#'\item{freq}{A vector with the frequencies of the alleles of the focal marker.}
#'\item{ehh}{A data frame with EHH values for each allele of the focal marker.}
#'\item{ihh}{A vector with iHH (integrated EHH) values for each allele of the focal marker.}
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
#'@seealso \code{\link{data2haplohh}}, \code{\link{plot.ehh}}, \code{\link{calc_ehhs}}, \code{\link{scan_hh}}.
#'@examples
#'#example haplohh object (280 haplotypes, 1424 SNPs)
#'#see ?haplohh_cgu_bta12 for details
#'data(haplohh_cgu_bta12)
#'#computing EHH statistics for the marker "F1205400"
#'#which displays a strong signal of selection
#'ehh <- calc_ehh(haplohh_cgu_bta12, mrk = "F1205400")
#'@export
calc_ehh <-
  function(haplohh,
           mrk,
           limhaplo = 2,
           limhomohaplo = 2,
           limehh = 0.05,
           include_zero_values = FALSE,
           include_nhaplo = FALSE,
           phased = TRUE,
           polarized = TRUE,
           scalegap = NA,
           maxgap = NA,
           discard_integration_at_border = TRUE,
           lower_y_bound = limehh) {
    # check parameters
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
    if (limehh < 0 |
        limehh > 1) {
      stop("limehh must lie between 0 and 1.", call. = FALSE)
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
    
    # calculation
    t <- tabulate(haplo(haplohh)[, mrk] + 1L)
    
    ## order alleles by their internal coding
    mrk_alleles <- sort(unique(haplo(haplohh)[, mrk]))
    
    ## first allele is either ancestral or major allele
    first_allele <- ifelse(polarized, 0L, mrk_alleles[1L])
    other_alleles <- mrk_alleles[mrk_alleles != first_allele]
    
    ## create result data frame with one column
    ehh <- data.frame(positions(haplohh))
    
    # if there is no ancestral allele, prepend zeros
    if (!(first_allele %in% mrk_alleles)) {
      freq <- 0
      ehh <- cbind(ehh, 0)
      nhaplo <- data.frame(rep(0L, nmrk(haplohh)))
      ihh <- 0
    } else{
      freq <- t[1] / sum(t)
      res.list <- .Call(
        "CALL_EHH",
        haplo(haplohh),
        nhap(haplohh),
        nmrk(haplohh),
        mrk,
        first_allele,
        limhaplo,
        limhomohaplo,
        limehh,
        phased
      )
      
      ehh <- cbind(ehh, res.list[[2]])
      nhaplo <- data.frame(res.list[[1]])
      
      ihh <- .Call(
        "CALL_INTEGRAL",
        positions(haplohh),
        res.list[[2]],
        mrk,
        limehh,
        scalegap,
        maxgap,
        discard_integration_at_border,
        lower_y_bound
      )
    }
    
    # if there are no derived alleles, append zeros
    if (length(other_alleles) == 0) {
      freq <- c(freq, 0)
      ehh <- cbind(ehh, 0)
      nhaplo <- cbind(nhaplo, 0L)
      ihh <- c(ihh, 0)
    } else{
      for (allele in other_alleles) {
        freq <- c(freq, t[allele + 1] / sum(t))
        res.list <- .Call(
          "CALL_EHH",
          haplo(haplohh),
          nhap(haplohh),
          nmrk(haplohh),
          mrk,
          allele,
          limhaplo,
          limhomohaplo,
          limehh,
          phased
        )
        
        ehh <- cbind(ehh, res.list[[2]])
        nhaplo <- cbind(nhaplo, res.list[[1]])
        
        ihh_allele <- .Call(
          "CALL_INTEGRAL",
          positions(haplohh),
          res.list[[2]],
          mrk,
          limehh,
          scalegap,
          maxgap,
          discard_integration_at_border,
          lower_y_bound
        )
        
        ihh <- c(ihh, ihh_allele)
      }
    }
    
    first_name <- ifelse(polarized, "A", "MAJ")
    other_name <- ifelse(polarized, "D", "MIN")
    
    if (length(other_alleles) > 1) {
      index_other <- seq_along(other_alleles)
    } else{
      index_other <- ""
    }
    
    colnames(ehh) <- c("POSITION",
                       paste0("EHH_", first_name),
                       paste0("EHH_", other_name, index_other))
    colnames(nhaplo) <- c(paste0("NHAPLO_", first_name),
                          paste0("NHAPLO_", other_name, index_other))
    
    if (include_nhaplo) {
      ehh <- cbind(ehh, nhaplo)
    }
    
    row.names(ehh) <- mrk.names(haplohh)
    
    if (!include_zero_values) {
      nonzeros <- which(rowSums(ehh[2:(1 + length(freq))]) > 0)
      
      if (length(nonzeros) > 0) {
        first <- min(nonzeros)
        last <- max(nonzeros)
        
        #if limehh is zero report only the positions next to non-zero positions
        if (limehh == 0) {
          first <- max(first - 1, 1)
          last <- min(last + 1, nmrk(haplohh))
        }
        
        ehh <- ehh[first:last, ]
      } else{
        ehh <- ehh[mrk, ]
      }
    }
    
    
    names(ihh) <- c(paste0("IHH_", first_name),
                    paste0("IHH_", other_name, index_other))
    
    names(freq) <- c(paste0("FREQ_", first_name),
                     paste0("FREQ_", other_name, index_other))
    
    
    #output
    l <- list(
      mrk.name = ifelse(
        is.null(mrk.names(haplohh)),
        as.character(mrk),
        mrk.names(haplohh)[mrk]
      ),
      freq = freq,
      ehh = ehh,
      ihh = ihh
    )
    
    return(new("ehh", l))
  }

#'@rdname calc_ehh
#'@aliases calc_ehh
#'@importFrom methods setClass
ehh <- setClass("ehh",
                contains = "list")

#'@importFrom methods is validObject
is.ehh <- function(x) {
  res <- (is(x, "ehh") & validObject(x))
  return(res)
}
