#'Class "haplohh"
#'@description An object of this class contains the information needed for computation of EHH based statistics.
#'@slot chr.name name of the chromosome/scaffold to which the markers belong.
#'@slot positions vector of type \code{numeric} containing the marker positions within the chromosome.
#'@slot haplo matrix of type \code{integer} containing haplotypes in rows and markers in columns.
#'@details This class is the basis for all calculations done by this package.
#'Note that the matrix in slot \code{haplo} has to be of type \code{integer}, not \code{numeric}.
#'Objects built by versions of rehh up to 2.0.4 coded this matrix as \code{numeric} and used
#'a different coding scheme. They can be converted e.g. by
#'\code{haplohh <- update_haplohh(old_haplohh)} in order be used
#'with the present version.
#'@seealso \code{\link{data2haplohh}}, \code{\link{update_haplohh}}
#'@examples showClass("haplohh")
#'@export
#'@importFrom methods setClass
#'@aliases chr.name positions haplo nmrk nhap mrk.names hap.names
setClass(
  Class = "haplohh",
  slots = c(
    chr.name = "character",
    positions = "numeric",
    haplo = "matrix"
  )
)

setGeneric("chr.name", function(x)
  standardGeneric("chr.name"))

#'@rdname haplohh-class
#'@aliases chr.name, haplohh-class
#'@param x an object of this class.
#'@export
setMethod("chr.name", "haplohh", function(x)
  x@chr.name)

setGeneric("positions", function(x)
  standardGeneric("positions"))

#'@rdname haplohh-class
#'@aliases positions, haplohh-class
#'@export
setMethod("positions", "haplohh", function(x)
  x@positions)

setGeneric("haplo", function(x)
  standardGeneric("haplo"))

#'@rdname haplohh-class
#'@aliases haplo, haplohh-class
#'@export
setMethod("haplo", "haplohh", function(x)
  x@haplo)

setGeneric("nmrk", function(x)
  standardGeneric("nmrk"))

#'@rdname haplohh-class
#'@aliases nmrk, haplohh-class
#'@export
setMethod("nmrk", "haplohh", function(x)
  ncol(x@haplo))

setGeneric("mrk.names", function(x)
  standardGeneric("mrk.names"))

#'@rdname haplohh-class
#'@aliases mrk.names, haplohh-class
#'@export
setMethod("mrk.names", "haplohh", function(x)
  colnames(x@haplo))

setGeneric("nhap", function(x)
  standardGeneric("nhap"))

#'@rdname haplohh-class
#'@aliases nhap, haplohh-class
#'@export
setMethod("nhap", "haplohh", function(x)
  nrow(x@haplo))

setGeneric("hap.names", function(x)
  standardGeneric("hap.names"))

#'@rdname haplohh-class
#'@aliases hap.names, haplohh-class
#'@export
setMethod("hap.names", "haplohh", function(x)
  rownames(x@haplo))

#'@importFrom methods is validObject
is.haplohh <- function(x) {
  check <- (is(x, "haplohh") & validObject(x))
  ## explicitly check if haplo matrix is integer
  if (check) {
    check <- is.integer(x@haplo)
    if (!check) {
      warning("Matrix 'haplo' of class 'haplohh' must be of type 'integer'.")
    }
  }
  ## check whether number of markers is equal to number of positions
  if (check) {
    check <- nmrk(x) == length(positions(x))
    if (!check) {
      warning("Number of positions must be equal to number of markers.")
    }
  }
  return(check)
}
