#'class for haplotype length
#'@export
haplen <- setClass("haplen",
                   contains = "list")

is.haplen <- function(x){
  res <- (is(x, "haplen") & validObject(x))
  return(res)
}

#'Calculate length of longest shared haplotypes around a focal marker
#'@description Calculate for each chromosome the maximum length of its extended haplotype homozygosity.
#'@param furcation an object of class \code{furcation} calculated by \code{\link{calc_furcation}}.
#'@details Extended haplotype homozygosity is defined as the region
#'around a focal marker in which a particular chromosome shares
#'a haplotype with (its sequence is identical to) another chromosome.
#'The function calculates for each chromosome the boundaries of its longest
#'shared haplotype. These correspond to the last furcations of a chromsome 
#'in a furcation diagram. Note that the calculation is performed independently
#'upstream and downstream of the focal marker and hence upper and lower 
#'boundaries do not necessarily arise from the same chromosomal pair.
#'@return The functions returns a list containing four elements:
#'\describe{
#'\item{mrk.name}{name/identifier of the focal marker.}
#'\item{position}{position of the focal marker.}
#'\item{xlim}{positions of left- and rightmost markers covered by extended haplotypes.}
#'\item{haplen}{a data frame with the coordinates of extended haplotypes around the focal marker.}
#'}
#'@examples #example haplohh object (280 haplotypes, 1424 SNPs)
#'#see ?haplohh_cgu_bta12 for details
#'data(haplohh_cgu_bta12)
#'#plotting haplotype lengths for both ancestral and derived allele
#'#of the marker "F1205400"
#'#which displays a strong signal of selection
#'f <- calc_furcation(haplohh_cgu_bta12, mrk = "F1205400")
#'h <- calc_haplen(f)
#'plot(h)
#'@export
calc_haplen <- function(furcation) {
  ##check parameters
  if (!is.furcation(furcation)) {
    stop("The data is not a valid object of a furcation.", call. = FALSE)
  }

  ##perform calculation
  hapl_lim_l <- rep(NA, furcation@nhap)
  hapl_lim_r <- rep(NA, furcation@nhap)
  allele <- rep(NA, furcation@nhap)
  description <- rep(NA, furcation@nhap)

  #number of chromosomes without missing values at focal marker
  nhap_without_missing <- 0
  for (f in furcation) {
    nhap_without_missing <- nhap_without_missing + f@count
  }

  for (i in seq_along(furcation)) {
    ftree <- furcation[[i]]@left
    #x coordinates of all inner nodes
    x <- ftree@node_pos
    label_parent <- ftree@label_parent
    node_parent <- ftree@node_parent
    node_size <- calc_node_size(node_parent, label_parent)

    for (j in seq_along(label_parent)) {
      if (!is.na(label_parent[j])) {
        if (node_size[label_parent[j]] == 1) {
          hapl_lim_l[j] <- x[label_parent[j]]
        } else{
          hapl_lim_l[j] <- furcation@xlim[1]
        }
        allele[j] <- furcation[[i]]@allele
        description[j] <- furcation[[i]]@description
      }
    }

    ftree <- furcation[[i]]@right
    #x coordinates of all inner nodes
    x <- ftree@node_pos
    label_parent <- ftree@label_parent
    node_parent <- ftree@node_parent
    node_size <- calc_node_size(node_parent, label_parent)

    for (j in seq_along(label_parent)) {
      if (!is.na(label_parent[j])) {
        if (node_size[label_parent[j]] == 1) {
          hapl_lim_r[j] <- x[label_parent[j]]
        } else{
          hapl_lim_r[j] <- furcation@xlim[2]
        }
        allele[j] <- furcation[[i]]@allele
        description[j] <- furcation[[i]]@description
      }
    }
  }

  ## output
  l <- list(
    mrk.name = furcation@mrk.name,
    position = furcation@position,
    xlim = furcation@xlim,
    haplen =  data.frame(
      ALLELE = allele,
      DESCRIPTION = description,
      MIN = hapl_lim_l,
      MAX = hapl_lim_r
    )
  )

  h <- new("haplen", l)

  return(h)
}
