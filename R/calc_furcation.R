#'calculate furcation trees around a focal marker
#'@description  Calculate furcation trees around a focal marker. A furcation tree captures
#'in greater detail than EHH values the decrease of extended haplotype homozygosity at
#'increasing distances from the selected focal marker.
#'@param haplohh an object of class haplohh (see \code{\link{data2haplohh}}).
#'@param mrk integer representing the number of the focal marker within the haplohh object
#'or string representing its ID/name.
#'@param allele a vector of alleles as coded internally, i.e. in case of polarized alleles,
#'0 represents the ancestral, 1 or higher the derived alleles.
#'If \code{NULL}, all alleles of the focal marker are considered.
#'@param limhaplo if there are less than \code{limhaplo} chromosomes that can be used for
#'the calculation, it is stopped. This is useful in case of missing data,
#'which lead to a successive exclusion of haplotypes: the further away from the focal marker
#'the less haplotypes are evaluated.
#'@param phased logical. If \code{TRUE} (default), chromosomes are expected to be phased.
#'If \code{FALSE}, consecutive chromosomes are assumed to
#'belong to diploid individuals and furcation trees are limited to within individuals which
#'are homozygous at the focal marker.
#'@param polarized logical. Affects only the order of furcations. If \code{TRUE} (default), the ancestral allele
#'becomes the first furcation and derived alleles are sorted by their internal coding. Otherwise all alleles
#'are sorted by their internal coding.
#'@details A haplotype furcation tree visualizes the breakdown
#'of LD at increasing distances from the focal marker.
#'The root of each tree is an allele of the focal marker, which in turn is identified
#'by a vertical dashed line.
#'Moving either to the "left" or to the "right" of the focal marker, each further marker is an opportunity for a node;
#'the tree either divides or does not, based on whether alleles at that marker
#'distinguish between hitherto identical extended haplotypes.
#'The thickness of the lines corresponds to the number of chromosomes sharing an extended haplotype.
#'@return An object of class furcation, containing the furcation structure of the specified alleles at the focal marker.
#'@seealso \code{\link{plot.furcation}}, \code{\link{calc_haplen}}.
#'@references Sabeti, P.C. and Reich, D.E. and Higgins, J.M. and Levine, H.Z.P and Richter, D.J. and Schaffner, S.F. and Gabriel, S.B. and Platko, J.V. and Patterson, N.J. and McDonald, G.J. and Ackerman, H.C. and Campbell, S.J. and Altshuler, D. and Cooper, R. and Kwiatkowski, D. and Ward, R. and Lander, E.S. (2002). Detecting recent positive selection in the human genome from haplotype structure. Nature, 419, 832-837.
#'@examples #example haplohh object (280 haplotypes, 1424 SNPs)
#'#see ?haplohh_cgu_bta12 for details
#'data(haplohh_cgu_bta12)
#'#plotting a furcation diagram for both ancestral and derived allele
#'#from the marker "F1205400"
#'#which display a strong signal of selection
#'f <- calc_furcation(haplohh_cgu_bta12, mrk = "F1205400")
#'plot(f)
#'@export
#'@importFrom methods new
calc_furcation <-
  function(haplohh,
           mrk,
           allele = NA,
           limhaplo = 2,
           phased = TRUE,
           polarized = TRUE) {
    ##check parameters
    if (!(is.haplohh(haplohh))) {
      stop("Data is not a valid haplohh object.", call. = FALSE)
    }
    if (limhaplo < 2) {
      stop("limhapcount must be an integer greater than 1.", call. = FALSE)
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
            " is greater than the number of markers in the data set (",
            nmrk(haplohh),
            ")."
          ),
          call. = FALSE
        )
      }
    } else{
      mrk <- as.character(mrk)
      if (!(mrk %in% mrk.names(haplohh))) {
        stop(paste0("Marker '",
                    mrk,
                    "' not found."), call. = FALSE)
      }
      mrk <- which(mrk.names(haplohh) == mrk)
    }
    
    ## order alleles by their internal coding
    mrk_allele <- sort(unique(haplo(haplohh)[, mrk]))
    
    ## create description of alleles ("Ancestral", "Major", etc.)
    if (polarized) {
      if (sum(mrk_allele != 0L) > 1) {
        index_other <- seq_len(sum(mrk_allele != 0L))
      } else{
        index_other <- ""
      }
      description <- paste0("Derived", index_other)
      
      # if ancestral allele is present, set it first
      if (0L %in% mrk_allele) {
        mrk_allele <- c(0L, mrk_allele[mrk_allele != 0L])
        description <- c("Ancestral", description)
      }
      
    } else{
      if (length(mrk_allele) > 2) {
        index_other <- seq_along(mrk_allele[-1])
      } else{
        index_other <- ""
      }
      description <- c("Major", paste0("Minor", index_other))
    }
    
    if (anyNA(allele) | length(allele) == 0) {
      allele <- mrk_allele
    } else{
      if (!is.numeric(allele)) {
        stop("Allele has to be specified by an integer.")
      }
      allele <- as.integer(allele)
      for (i in allele) {
        if (!(i %in% mrk_allele)) {
          stop(paste0("Marker has no allele '", i, "'."))
        }
      }
    }
    
    ##perform calculation
    
    f <- furcation(
      mrk.name = ifelse(
        is.null(mrk.names(haplohh)),
        as.character(mrk),
        mrk.names(haplohh)[mrk]
      ),
      position = positions(haplohh)[mrk],
      xlim = range(positions(haplohh)),
      nhap = nhap(haplohh)
    )
    
    #calculate allele furcations
    for (a in allele) {
      # calculate ftree for left and right side
      allelefurcation <- new("allelefurcation",
                             allele = a,
                             description = description[which(mrk_allele == a)])
      
      for (side in 1:2) {
        # first resp. last marker defines side
        allelefurcation_list <- .Call(
          "CALL_FURCATION",
          haplohh@haplo,
          nhap(haplohh),
          mrk,
          ifelse(side == 1, 1L, nmrk(haplohh)),
          a,
          limhaplo,
          phased
        )
        
        ftree <- new("ftree")
        
        #change indexing from C to R
        ftree@node_parent <- allelefurcation_list[[2]] + 1L
        #replace marker index by its chromosomal position
        ftree@node_pos <-
          positions(haplohh)[allelefurcation_list[[1]] + 1L]
        #transform to logical
        ftree@node_with_missing_data <-
          as.logical(allelefurcation_list[[3]])
        #change indexing from C to R
        ftree@label_parent <- allelefurcation_list[[4]] + 1L
        
        #add node "names" for human readability;
        #this is NOT needed for computation
        names(ftree@node_pos) <- seq_along(ftree@node_pos)
        names(ftree@node_parent) <- seq_along(ftree@node_parent)
        names(ftree@node_with_missing_data) <-
          seq_along(ftree@node_with_missing_data)
        names(ftree@label_parent) <- seq_along(ftree@label_parent)
        
        if (side == 1) {
          allelefurcation@left <- ftree
        } else{
          allelefurcation@right <- ftree
        }
      }
      
      allelefurcation@count <- as.integer(sum(!(
        is.na(allelefurcation@left@label_parent) &
          is.na(allelefurcation@right@label_parent)
      )))
      
      if (allelefurcation@count > 0) {
        f[[as.character(a)]] <- allelefurcation
      }
    }
    
    return(f)
  }
