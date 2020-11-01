#'Calculate pairwise shared haplotype length between all chromosomes
#'@description Calculate pairwise shared haplotype length between all chromosomes at a focal marker.
#'@param haplohh an object of class \code{haplohh} (see \code{\link{data2haplohh}}).
#'@param mrk integer representing the number of the focal marker within the haplohh object
#'or string representing its ID/name.
#'@param phased logical. If \code{TRUE} (default) chromosomes are expected to be phased. If \code{FALSE}, the haplotype data is assumed to
#'consist of pairwise ordered chromosomes belonging to diploid individuals and only the two chromosomes of
#'each individual are compared.
#'@param maxgap maximum allowed gap in bp between two markers. If exceeded, further calculation is stopped at the gap
#'(default=\code{NA}, i.e. no limitation).
#'@param max_extend maximum distance in bp to extend shared haplotypes away from the focal marker.
#'(default \code{NA}, i.e. no limitation).
#'@param side side to consider, either "left" (positions lower than focal position), "right" (positions higher than focal position)
#'or "both" (default).
#'@details The function computes the length of shared haplotypes (stretches of identical sequence) around
#'the focal marker.
#'
#'Note that the function \code{\link{calc_haplen}} calculates for each chromosome
#'the boundaries of its longest shared haplotype; separately upstream and downstream of
#'the focal marker.
#'@return The returned value is a matrix with pairwise shared haplotype lengths.
#'@seealso \code{\link{data2haplohh}}, \code{\link{scan_hh_full}}.
#'@examples
#'#example haplohh object (280 haplotypes, 1424 SNPs)
#'#see ?haplohh_cgu_bta12 for details
#'data(haplohh_cgu_bta12)
#'#computing shared haplotype lengths around the marker "F1205400"
#'#which displays a strong signal of selection
#'m <- calc_pairwise_haplen(haplohh_cgu_bta12, mrk = "F1205400")
#'@export
calc_pairwise_haplen <-
  function(haplohh,
           mrk,
           phased = TRUE,
           maxgap = NA,
           max_extend = NA,
           side = "both") {
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
    
    if (!(side %in% c("both", "left", "right"))) {
      stop("side must be either \"both\", \"left\" or \"right\".", call. = FALSE)
    }
    
    ##perform calculations
    pairwise_haplen <-
      matrix(
        0,
        nrow = nhap(haplohh),
        ncol = nhap(haplohh),
        dimnames = list(row.names(haplo(haplohh)), row.names(haplo(haplohh)))
      )
    
    .Call(
      "CALL_PAIRWISE_HAPLEN",
      haplo(haplohh),
      nhap(haplohh),
      nmrk(haplohh),
      positions(haplohh),
      mrk,
      as.integer(maxgap),
      as.integer(max_extend),
      ifelse(side == "both", 0L, ifelse(side == "left", 1L, 2L)),
      as.integer(phased),
      pairwise_haplen
    )
    
    ##output
    return(pairwise_haplen)
  }
