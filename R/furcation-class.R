#'An S4 class to represent a furcation tree on one side of one allele of a focal marker
#'@details A furcation structure consists of two trees ("left" and "right") for each allele of
#'a focal marker. If there are only bi-allelic markers
#'and no missing values, the trees are bifurcating.
#'
#'Missing values are treated similarly to an extra allele
#'in so far as they cause a furcation. However, the resulting daughter node is marked
#'accordingly and the chromosomes excluded from further calculations.
#'If all chromosomes of a parent node have missing values, the "furcation" is
#'degenerated and yields a single daughter node.
#'
#'Note that a tree with n leaves can have at most 2n-1 nodes.
#'
#'In a furcation tree, the leaves do not necessarily represent
#'single chromosomes, either due to multiple missing data or
#'because the first/last marker was reached before all
#'extended haplotypes were distinct.
#'
#'@slot node_parent a vector, representing the tree structure. 
#'Each node (number) is assigned its parent node (number).
#'@slot node_pos a vector, assigning to each node (number) its position in the chromosome, i.e.
#'at which marker position the furcation occurred.
#'@slot node_with_missing_data a vector of type \code{logical}. 
#'Pseudo-furcations arise due to missing data at a marker.
#'The daughter node (number) is marked accordingly.
#'@slot label_parent a vector, that attaches an "extra leave", representing
#'the haplotype number (defined by the order in the haplotype data file) to leaves
#'of the tree.
#'This is necessary because in general not all leaves of the original tree represent
#'a single haplotype/chromosome.
ftree <- setClass(
  "ftree",
  slots = c(
    node_parent = "integer",
    node_pos = "numeric",
    node_with_missing_data = "logical",
    label_parent = "integer"
  )
)

#'An S4 class containing furcation trees for one allele
#'of a focal marker
#'@description An S4 class containing the furcation trees
#'for both sides of a focal marker for one allele.
#'@slot allele the allele of the focal marker.
#'@slot description "ancestral", "derived", "major", "minor", etc.
#'@slot count the number of chromosomes with that allele.
#'@slot left furcation tree to the left of the marker.
#'@slot right furcation tree to the right of the marker.
#'@seealso \code{\link{ftree}}, \code{\link{furcation}}
allelefurcation <- setClass(
  "allelefurcation",
  slots = c(
    allele = "integer",
    description = "character",
    count = "integer",
    left = "ftree",
    right = "ftree"
  )
)

#'An S4 class representing the complete furcation pattern around a focal marker.
#'@slot .Data a list containing for each allele an object of \code{allelefurcation-class}.
#'@slot mrk.name the name/identifier of the focal marker.
#'@slot position the chromosomal position of the focal marker.
#'@slot xlim the range of marker positions.
#'@slot nhap the number of haplotypes in the sample.
#'@seealso \code{\link{calc_furcation}}
#'@examples # copy example files into working directory
#'make.example.files()
#'# read first example file
#'hh <- data2haplohh("example1.vcf")
#'# remove example files
#'remove.example.files()
#'# calculate furcation structure around marker "rs6"
#'f <- calc_furcation(hh, mrk = "rs6")
#'# extract left side tree of ancestral allele (which is coded by '0')
#'f[['0']]@@left
#'# the tree consists of seven nodes, '1' being the root node
#'# nodes 2 and 3 have the root node as parent, etc.
#'# the first chromosome is attached as a label node to node 7, etc.
#'# For comparison, a plot of the complete furcation structure:
#'plot(f)
#'@export
furcation <- setClass(
  "furcation",
  contains = "list",
  slots = c(
    mrk.name = "character",
    position = "numeric",
    xlim = "numeric",
    nhap = "integer"
  )
)

#'@importFrom methods is validObject
is.ftree <- function(x) {
  res <- (is(x, "ftree") & validObject(x))
  return(res)
}

#'@importFrom methods is validObject
is.allelefurcation <- function(x) {
  res <- (is(x, "allelefurcation") & validObject(x))
  return(res)
}

#'@importFrom methods is validObject
is.furcation <- function(x) {
  res <- (is(x, "furcation") & validObject(x))
  return(res)
}
