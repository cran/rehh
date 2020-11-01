#'Subsets object of \code{\link{haplohh-class}}
#'@description Subsets the data of an object of class \code{\link{haplohh-class}},
#'meeting certain conditions.
#'@param x object of class \code{\link{haplohh-class}} to be subset.
#'@param select.hap expression, indicating haplotypes to select.
#'@param select.mrk expression, indicating markers to select.
#'@param min_perc_geno.hap threshold on percentage of missing data for haplotypes
#'(haplotypes with less than \code{min_perc_geno.hap} percent of markers genotyped are discarded). Default is \code{NA},
#'hence no constraint.
#'@param min_perc_geno.mrk threshold on percentage of missing data for markers (markers genotyped on less than
#'\code{min_perc_geno.mrk} percent of haplotypes are discarded). By default, \code{min_perc_geno.mrk=100},
#'hence only fully genotyped markers are retained.
#'This value cannot be set to \code{NA} or zero.
#'@param min_maf threshold on the Minor Allele Frequency. Markers having a MAF lower than or equal to minmaf are discarded.
#'In case of multi-allelic markers the second-most frequent allele is referred to as minor allele.
#'Setting this value to zero eliminates monomorphic sites. Default is \code{NA},
#'hence no constraint.
#'@param max_alleles threshold for the maximum number of different alleles at a site. Default is \code{NA},
#'hence no restriction. In order to retain only bi-allelic markers, set this parameter to 2.
#'@param verbose logical. If \code{TRUE} (default), report verbose progress.
#'@param ... further arguments are ignored.
#'@seealso \code{\link{haplohh-class}}, \code{\link{data2haplohh}}
#'@examples #example haplohh object (280 haplotypes, 1424 SNPs)
#'#see ?haplohh_cgu_bta12 for details
#'data(haplohh_cgu_bta12)
#'#select subset of first 10 hyplotypes and first 5 markers
#'subset(haplohh_cgu_bta12, select.hap = 1:10, select.mrk = 1:5)
#'@export
subset.haplohh <-
  function(x,
           select.hap = NULL,
           select.mrk = NULL,
           min_perc_geno.hap = NA,
           min_perc_geno.mrk = 100,
           min_maf = NA,
           max_alleles = NA,
           verbose = TRUE,
           ...) {
    # check parameters
    ### empty haplotype makes no sense, but is of no harm
    if (!is.na(min_perc_geno.hap) & (min_perc_geno.hap < 0 |
                                     min_perc_geno.hap > 100)) {
      stop("min_perc_geno.hap should lie in the interval [0,100].",
           call. = FALSE)
    }
    ### empty marker will cause trouble -> forbid NA or zero !!
    if (is.na(min_perc_geno.mrk) | min_perc_geno.mrk <= 0 |
        min_perc_geno.mrk > 100) {
      stop("min_perc_geno.mrk should lie in the interval (0,100].",
           call. = FALSE)
    }
    ### minor frequency can be maximal 0.5
    if (!is.na(min_maf)) {
      if (!is.numeric(min_maf) | min_maf < 0 |
          min_maf > 0.5) {
        stop("min_maf should lie in the interval [0,0.5].", call. = FALSE)
      }
    }
    ### max_alleles must be at least 2
    if (!is.na(max_alleles)) {
      if (!is.numeric(max_alleles) | max_alleles < 2) {
        stop("max_alleles should be at least 2.", call. = FALSE)
      }
    }
    ### check if object is valid haplohh
    if (!is.haplohh(x)) {
      stop("Data is not a valid object of class haplohh.", call. = FALSE)
    }
    ### check if object is empty
    if (min(dim(haplo(x))) == 0) {
      stop("Empty data object.", call. = FALSE)
    }
    
    # filtering
    
    if (!is.null(select.hap)) {
      if (verbose)
        cat("* Subset haplotypes *\n")
      # drop = FALSE ensures that haplo is a matrix
      haplo <- x@haplo[select.hap, , drop = FALSE]
      
      if (nrow(haplo) == nhap(x)) {
        if (verbose)
          cat("No haplotype discarded.\n")
      } else{
        if (verbose)
          cat(nhap(x) - nrow(haplo), "haplotypes discarded.\n")
        x@haplo <- haplo
        if (verbose)
          cat(nhap(x), "haplotypes remaining.\n")
      }
      if (nhap(x) == 0) {
        warning(
          "No haplotype left after filtering on 'select.hap'.\n",
          call. = FALSE,
          immediate. = TRUE
        )
        return(x)
      }
    }
    
    if (!is.null(select.mrk)) {
      if (verbose)
        cat("* Subset markers *\n")
      # drop = FALSE ensures that haplo is a matrix
      haplo <- x@haplo[, select.mrk, drop = FALSE]
      
      if (ncol(haplo) == nmrk(x)) {
        if (verbose)
          cat("No marker discarded.\n")
      } else{
        if (verbose)
          cat(nmrk(x) - ncol(haplo), "markers discarded.\n")
        # positions are unnamed in haplohh
        positions <- positions(x)
        # name them to subset by marker name
        names(positions) <- mrk.names(x)
        # update slots in object
        x@haplo <- haplo
        x@positions <- unname(positions[select.mrk])
        if (verbose)
          cat(nmrk(x), "markers remaining.\n")
        if (nmrk(x) == 0) {
          warning(
            "No marker left after filtering on 'select.mrk'.\n",
            call. = FALSE,
            immediate. = TRUE
          )
          return(x)
        }
      }
    }
    
    if (verbose)
      cat("* Filtering data *\n")
    
    if (!is.na(min_perc_geno.hap)) {
      if (verbose)
        cat("Discard haplotypes with less than",
            min_perc_geno.hap,
            "% of genotyped markers.\n")
      
      hap_sel <-
        (100 * rowMeans(!is.na(x@haplo))) >= min_perc_geno.hap
      
      if (sum(hap_sel) == nhap(x)) {
        if (verbose)
          cat("No haplotype discarded.\n")
      } else{
        if (verbose)
          cat(nhap(x) - sum(hap_sel), "haplotypes discarded.\n")
        x@haplo <- x@haplo[hap_sel, , drop = FALSE]
        if (verbose)
          cat(nhap(x), "haplotypes remaining.\n")
        
        if (nhap(x) == 0) {
          warning(
            "No haplotype left after filtering of missing data.\n",
            "If applicable, reduce min_perc_geno.hap to allow for more missing data.\n",
            call. = FALSE,
            immediate. = TRUE
          )
          return(x)
        }
      }
    }
    #selection des markers d'apres les donnees manquantes
    if (verbose)
      cat("Discard markers genotyped on less than",
          min_perc_geno.mrk,
          "% of haplotypes.\n")
    mrk_sel <-
      (100 * colMeans(!is.na(x@haplo))) >= min_perc_geno.mrk
    if (sum(mrk_sel) == nmrk(x)) {
      if (verbose)
        cat("No marker discarded.\n")
    } else{
      if (verbose)
        cat(nmrk(x) - sum(mrk_sel), "markers discarded.\n")
      x@haplo <- x@haplo[, mrk_sel, drop = FALSE]
      x@positions <- x@positions[mrk_sel]
      if (verbose)
        cat(nmrk(x), "markers remaining.\n")
      
      if (nmrk(x) == 0) {
        warning(
          "No marker left after filtering of missing data.\n",
          "If applicable, reduce min_perc_geno.mrk to allow for more missing data.\n",
          call. = FALSE,
          immediate. = TRUE
        )
        return(x)
      }
    }
    #selection des mrks sur MAF
    if (!is.na(min_maf)) {
      if (verbose)
        cat("Discard markers with Minor Allele Frequency equal to or below",
            min_maf,
            ".\n")
      mrk_sel <- apply(x@haplo, 2, function(x) {
        t <- tabulate(x + 1L)
        if (length(t) == 1) {
          ## only one (ancestral) allele
          return(FALSE)
        } else{
          ## at least one (derived) allele
          alleles <- order(t, decreasing = TRUE)[1:2]
          #check frequency of most and second-most common allele
          return(min(t[alleles]) / sum(t) > min_maf)
        }
      })
      if (sum(mrk_sel) == nmrk(x)) {
        if (verbose)
          cat("No marker discarded.\n")
      } else{
        if (verbose)
          cat(nmrk(x) - sum(mrk_sel), "markers discarded.\n")
        x@haplo <- x@haplo[, mrk_sel, drop = FALSE]
        x@positions <- x@positions[mrk_sel]
        if (verbose)
          cat(nmrk(x), "markers remaining.\n")
        
        if (nmrk(x) == 0) {
          warning(
            "No marker left after filtering on Minor Allele Frequency.\n",
            "If applicable, reduce min_maf to allow for less frequent minor alleles.",
            call. = FALSE,
            immediate. = TRUE
          )
          return(x)
        }
      }
    }
    
    if (!is.na(max_alleles)) {
      if (verbose)
        cat("Discard markers with more than",
            max_alleles,
            "different alleles.\n")
      mrk_sel <-
        apply(x@haplo, 2, function(x) {
          length(unique(na.omit(x))) <= max_alleles
        })
      if (sum(mrk_sel) == nmrk(x)) {
        if (verbose)
          cat("No marker discarded.\n")
      } else{
        if (verbose)
          cat(nmrk(x) - sum(mrk_sel), "markers discarded.\n")
        x@haplo <- x@haplo[, mrk_sel, drop = FALSE]
        x@positions <- x@positions[mrk_sel]
        if (verbose)
          cat(nmrk(x), "markers remaining.\n")
        
        if (nmrk(x) == 0) {
          warning(
            "No marker left after filtering on maximal number of different alleles.\n",
            call. = FALSE,
            immediate. = TRUE
          )
          return(x)
        }
      }
    }
    
    if (verbose)
      cat("Data consists of",
          nhap(x),
          "haplotypes and",
          nmrk(x),
          "markers.\n")
    
    multicity <- tabulate(apply(x@haplo, 2, function(x) {
      sum(tabulate(x + 1) != 0)
    }))
    names(multicity) <- seq_along(multicity)
    
    if (verbose) {
      cat("Number of mono-, bi-, multi-allelic markers:\n")
      cat(names(multicity), "\n")
      cat(multicity, "\n")
    }
    
    return(x)
  }
