#'Update object of class haplohh
#'@description Update object of class \code{\link{haplohh-class}}
#'constructed by rehh versions up to version 2.0.4.
#'@param haplohh an object of an old version of \code{\link{haplohh-class}}.
#'@details This function is intended to update \code{haplohh} objects
#'that have been built by rehh versions up to 2.0.4. These objects cannot
#'be used in functions of the current version. 
#'The following changes have been made to the class definition:
#'The internal representation of the haplotype matrix followed the encoding
#'\itemize{
#'\item 0 missing value
#'\item 1 ancestral allele
#'\item 2 derived allele}
#'and has been replaced by a vcf-like encoding:
#'\itemize{
#'\item \code{NA} missing value
#'\item 0 ancestral allele
#'\item 1 derived allele.}
#'Furthermore the slots \code{nsnp}, \code{snp.name} and \code{nhap} have been removed
#'and slot \code{position} renamed to \code{positions}.
#'An update of an old \code{haplohh} object is done as follows:
#'
#'\code{new_haplohh = update_haplohh(old_haplohh)}.
#'@seealso \code{\link{haplohh-class}}, \code{\link{data2haplohh}}.
#'@export
update_haplohh <- function(haplohh) {
  if (any(!(haplohh@haplo %in% 0:2))) {
    cat("Alleles are not coded in the appropriate format:\n")
    cat("0 (missing data), 1 (ancestral allele) or 2 (derived allele).\n")
    stop("Update failed.", call. = FALSE)
  }
  
  haplo <- matrix(
    as.integer(haplohh@haplo) - 1L,
    nrow = nrow(haplohh@haplo),
    ncol = ncol(haplohh@haplo)
  )
  haplo[haplo == -1] <- NA
  colnames(haplo) <- haplohh@snp.name
  
  return(
    new(
      "haplohh",
      chr.name = haplohh@chr.name,
      positions = haplohh@position,
      haplo = haplo
    )
  )
}
