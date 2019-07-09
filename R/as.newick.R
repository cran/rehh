#'Convert a furcation tree into Newick format
#'@description Convert a furcation tree into Newick format.
#'@param furcation an object of \code{\link{furcation-class}}.
#'@param allele the allele to be considered (default 0).
#'@param side side (either \code{"left"} or \code{"right"}).
#'@param hap.names names/labels of chromosomes in haplotype data file. Per default
#'haplotypes are numbered by their order in the input file.
#'@seealso \code{\link{ftree-class}}, \code{\link{calc_furcation}}, \code{\link{plot.furcation}}
#'@examples #example haplohh object (280 haplotypes, 1424 SNPs)
#'#see ?haplohh_cgu_bta12 for details
#'data(haplohh_cgu_bta12)
#'#calculate furcation for the marker "F1205400"
#'#which displays a strong signal of selection
#'f <- calc_furcation(haplohh_cgu_bta12, mrk = "F1205400")
#'#get left tree of ancestral allele (coded as '0')
#'as.newick(f, 0, "left")
#'@export
as.newick <-
  function(furcation,
           allele = 0,
           side,
           hap.names = seq_len(furcation@nhap)) {
    ##check parameters
    if (!is.furcation(furcation)) {
      stop("The data is not a valid furcation object.", call. = FALSE)
    }
    if (!is.numeric(allele) | round(allele) != allele) {
      stop("Allele has to be specified as integer number.", call. = FALSE)
    }
    allele <- as.character(allele)
    if (!(allele %in% names(furcation))) {
      stop(paste("Could not find allele", allele, "."), call. = FALSE)
    }
    if (length(allele) != 1 & length(furcation) > 1) {
      stop("Exactly one allele has to be specified.", call. = FALSE)
    }
    if (side != "left" & side != "right") {
      stop("Side has to be either 'left' or 'right'", call. = FALSE)
    }
    if (!is.null(hap.names)) {
      if (length(hap.names) != furcation@nhap) {
        stop(
          "Number of specified haplotype names has to be equal to number of haplotypes.",
          call. = FALSE
        )
      }
    } else{
      hap.names <- seq_len(furcation@nhap)
    }
    
    ##calculations
    if (side == "left") {
      ftree <- furcation[[allele]]@left
    } else{
      ftree <- furcation[[allele]]@right
    }
    
    ## usage of temporary file is a work-around for
    ## C string streams are missing under Windows 
    tmp_file_name <- tempfile()
    
    #calculation and output done by C
    if (.Call(
      "CALL_ASNEWICK",
      tmp_file_name,
      as.integer(ftree@node_parent - 1),
      as.integer(ftree@label_parent - 1),
      as.double(ftree@node_pos),
      as.double(ifelse(
        side == "left", furcation@xlim[1], furcation@xlim[2]
      )),
      as.character(hap.names)
    )) {
      newick <- readLines(tmp_file_name, 1, warn = FALSE)
    } else{
      stop("Could not write to a temporary file.", call. = FALSE)
    }
    
    unlink(tmp_file_name)
    return(newick)
  }
