#'Compute iHS
#'@description Compute iHS (standardized ratio of iHH values of two alleles).
#'@param scan a data frame with chromosome name,
#'marker position, frequency of ancestral (resp. major) allele, frequency of derived (resp. minor)
#'allele, and iHH for both alleles, as obtained from function \code{\link{scan_hh}}.
#'@param freqbin size of the bins to standardize log(iHH_A/iHH_D). Markers are binned with
#'respect to the derived allele frequency at the focal marker. The bins are built from
#'\code{min_maf} to \code{1-min_maf} in steps of size \code{freqbin}. If set to 0, standardization
#'is performed considering each observed frequency as a discrete frequency
#'class (useful in case of a large number of markers and few different haplotypes).
#'If set to an integer of 1 or greater, a corresponding number of equally sized bins are created.
#'@param min_maf focal markers with a MAF (Minor Allele Frequency) lower than or equal to \code{min_maf}
#'are discarded from the analysis (default 0.05).
#'@param min_nhaplo focal markers with least one of the two compared alleles carried by fewer
#'than \code{min_nhaplo} haplotypes, are discarded (default \code{NA}).
#'@param standardize logical. If \code{TRUE} (default), then standardize iHS, else report unstandardized iHS.
#'@param include_freq logical. If \code{TRUE} include columns with allele frequencies into result.
#'@param right logical. If \code{TRUE} the bin intervals are closed on the right (and open on the left).
#'@param alpha calculate quantiles \code{alpha/2} and \code{(1-alpha/2)} for unstandardized binned iHS.
#'@param p.side side to which refers the p-value. Default \code{NA}, meaning two-sided. Can be set
#'to \code{"left"} or \code{"right"}.
#'@param p.adjust.method method passed to function \code{\link[stats]{p.adjust}} to correct the p-value for
#' multiple testing. Default \code{"none"}.
#'@param verbose logical. If \code{TRUE} (default), report number of markers of the source data frame and result data frame.
#'@details Computes log ratio of iHH of two focal alleles as described in Voight et al. (2006). The standardization
#'is performed within each bins separately because of the frequency-dependence
#'of expected iHS values under neutrality. An implicit assumption of this
#'approach is that each bin is dominated by neutral markers.
#'
#'Since the standardized iHS values follow, if markers evolve predominantly neutrally, approximately
#'a standard Gaussian distribution, it is practical to assign to the values a p-value relative
#'to the null-hypothesis of neutral evolution. The parameter \code{p.side} determines
#'if the p-value is assigned to both sides of the distribution or to one side of interest.
#'@return The returned value is a list containing two elements
#'\describe{
#'\item{ihs}{a data frame with markers in rows and the columns for chromosome name, marker position,
#'iHS and, if standardized, p-value in a negative log10 scale. Optionally, allele frequencies are included.}
#'\item{frequency.class}{a data frame with bins in rows and columns for
#'the number of markers, mean uniHS, standard deviation uniHS, lower quantile uniHS, upper quantile uniHS.}
#'}
#'@references Gautier, M. and Naves, M. (2011). Footprints of selection in the ancestral admixture of a New World Creole cattle breed. \emph{Molecular Ecology}, \strong{20}, 3128-3143.
#'
#'Voight, B.F. and Kudaravalli, S. and Wen, X. and Pritchard, J.K. (2006). A map of recent positive selection in the human genome. \emph{Plos Biology}, \strong{4}, e72.
#'@seealso \code{\link{scan_hh}}, \code{\link{distribplot}}, \code{\link{freqbinplot}}, \code{\link{manhattanplot}}
#'@examples library(rehh.data)
#'data(wgscan.cgu)
#'#results from a genome scan (44,057 SNPs)
#'#see ?wgscan.eut and ?wgscan.cgu for details
#'wgscan.cgu.ihs <- ihh2ihs(wgscan.cgu)
#'@export
#'@importFrom stats sd quantile p.adjust p.adjust.methods
ihh2ihs <- function(scan,
                    freqbin = 0.025,
                    min_maf = 0.05,
                    min_nhaplo = NA,
                    standardize = TRUE,
                    include_freq = FALSE,
                    right = FALSE,
                    alpha = 0.05,
                    p.side = NA,
                    p.adjust.method = "none",
                    verbose = TRUE) {
  final_score_colname <- "IHS"
  #check parameters
  if (!("CHR" %in% colnames(scan))) {
    stop("Data does not contain a column named 'CHR'.", call. = FALSE)
  }
  if (!("POSITION" %in% colnames(scan))) {
    stop("Data does not contain a column named 'POSITION'.", call. = FALSE)
  }
  if (anyDuplicated(scan[c("CHR", "POSITION")])) {
    stop("Data contains duplicated chromosomal positions.", call. = FALSE)
  }
  if (is.na(freqbin) | !is.numeric(freqbin) | freqbin < 0) {
    stop("Specify a number or size of frequency bins.", call. = FALSE)
  }
  if (!is.na(min_nhaplo) &
      (round(min_nhaplo) != min_nhaplo | min_nhaplo < 2)) {
    stop("min_nhaplo must be an integer number greater than 1.")
  }
  # check if interval can be covered
  ## note: %% cannot be used, because i386 behaves differently
  if (freqbin > 0 &
      freqbin < 1) {
    n_bin <- (1 - 2 * min_maf) / freqbin
    if (abs(n_bin - round(n_bin)) > 1E-6) {
      stop(
        paste(
          "Cannot cover interval [",
          min_maf,
          ",",
          1 - min_maf,
          "] with bins of size",
          freqbin,
          "."
        ),
        call. = FALSE
      )
    }
  }
  if (standardize) {
    if (!is.na(p.side) & !(p.side %in% c("left", "right"))) {
      stop("'p.side' must be either NA, 'left' or 'right'.",
           call. = FALSE)
    }
    if (!(p.adjust.method %in% p.adjust.methods)) {
      stop("Unknown adjust method.", call. = FALSE)
    }
  }
  
  # ignore case of colnames
  uppercolnames <- toupper(colnames(scan))
  #in case of results from older package versions (bi-allelic sites),
  #rename ancestral column, add derived column
  freq1_col_nr <- which(colnames(scan) == "freq_A")
  if (length(freq1_col_nr) == 1) {
    colnames(scan)[freq1_col_nr] <- "FREQ_A"
    scan$FREQ_D <- 1 - scan$FREQ_A
    freq2_col_nr <- which(colnames(scan) == "FREQ_D")
  }
  else{
    freq1_col_nr <- which(uppercolnames == "FREQ_A")
    freq2_col_nr <- which(uppercolnames == "FREQ_D")
    if (length(freq1_col_nr) != 1 | length(freq2_col_nr) != 1) {
      freq1_col_nr <- which(uppercolnames == "FREQ_MAJ")
      freq2_col_nr <- which(uppercolnames == "FREQ_MIN")
      if (length(freq1_col_nr) != 1 | length(freq2_col_nr) != 1) {
        stop("Columns(s) with allele frequencies not found.", call. = FALSE)
      }
    }
  }
  
  if (!is.na(min_nhaplo)) {
    nhaplo1_col_nr <- which(uppercolnames == "NHAPLO_A")
    nhaplo2_col_nr <- which(uppercolnames == "NHAPLO_D")
    if (length(nhaplo1_col_nr) != 1 | length(nhaplo2_col_nr) != 1) {
      nhaplo1_col_nr <- which(uppercolnames == "NHAPLO_MAJ")
      nhaplo2_col_nr <- which(uppercolnames == "NHAPLO_MIN")
      if (length(nhaplo1_col_nr) != 1 |
          length(nhaplo2_col_nr) != 1) {
        stop("Column(s) with number of haplotypes not found.", call. = FALSE)
      }
    }
  }
  
  # try polarized iHH
  ihh1_col_nr <- which(uppercolnames == "IHH_A")
  ihh2_col_nr <- which(uppercolnames == "IHH_D")
  if (length(ihh1_col_nr) != 1 | length(ihh2_col_nr) != 1) {
    # try unpolarized iHH
    ihh1_col_nr <- which(uppercolnames == "IHH_MAJ")
    ihh2_col_nr <- which(uppercolnames == "IHH_MIN")
    # if not found ...
    if (length(ihh1_col_nr) != 1 | length(ihh2_col_nr) != 1) {
      stop("Could not find columns with iHH values.", call. = FALSE)
    }
    # unpolarized markers should not be binned
    if (freqbin != 1 & standardize) {
      warning(paste(
        "If alleles are unpolarized, 'freqbin' should be set to 1 (one bin)."
      ),
      call. = FALSE)
    }
  }
  
  # perform calculation
  ## subset to markers with high enough minor allele frequency
  if (!is.na(min_maf)) {
    if (verbose)
      cat("Discard focal markers with Minor Allele Frequency equal to or below",
          min_maf,
          ".\n")
    
    mrk_sel <-
      scan[[freq1_col_nr]] > min_maf &
      scan[[freq2_col_nr]] > min_maf
    
    if (sum(mrk_sel) == nrow(scan)) {
      if (verbose)
        cat("No marker discarded.\n")
    } else{
      if (verbose)
        cat(nrow(scan) - sum(mrk_sel), "markers discarded.\n")
      
      scan <- scan[mrk_sel,]
      
      if (verbose)
        cat(nrow(scan), "markers remaining.\n")
      
    }
  }
  
  ## subset to markers with enough haplotypes evaluated
  if (!is.na(min_nhaplo)) {
    if (verbose)
      cat(
        "Discard focal markers with less than",
        min_nhaplo,
        "evaluated haplotypes for an allele.\n"
      )
    
    mrk_sel <-
      scan[[nhaplo1_col_nr]] >= min_nhaplo &
      scan[[nhaplo2_col_nr]] >= min_nhaplo
    
    if (sum(mrk_sel) == nrow(scan)) {
      if (verbose)
        cat("No marker discarded.\n")
    } else{
      if (verbose)
        cat(nrow(scan) - sum(mrk_sel), "markers discarded.\n")
      
      scan <- scan[mrk_sel,]
      
      if (verbose)
        cat(nrow(scan), "markers remaining.\n")
      
    }
  }
  
  if (nrow(scan) == 0) {
    stop("Scan contains no marker.",
         call. = FALSE)
  }
  
  ## un-standardized log ratio
  un_log_ratio <- log(scan[[ihh1_col_nr]] / scan[[ihh2_col_nr]])
  un_log_ratio[un_log_ratio == "Inf" | un_log_ratio == "-Inf"] <- NA
  
  if (freqbin > 0) {
    if (freqbin >= 1) {
      freqbin <- (1 - 2 * min_maf) / round(freqbin)
    }
    bins <- cut(scan[[freq2_col_nr]],
                breaks = seq(min_maf, 1 - min_maf, freqbin),
                right = right)
  } else {
    bins <- as.factor(scan[[freq2_col_nr]])
  }
  
  frequency.class <- as.data.frame(cbind(
    tapply(un_log_ratio, bins, length),
    tapply(un_log_ratio, bins, mean, na.rm = TRUE),
    tapply(un_log_ratio, bins, sd, na.rm = TRUE),
    tapply(
      un_log_ratio,
      bins,
      quantile,
      probs = alpha / 2,
      na.rm = TRUE
    ),
    tapply(
      un_log_ratio,
      bins,
      quantile,
      probs = 1 - alpha / 2,
      na.rm = TRUE
    )
  ))
  
  #empty classes get NA as length; replace by 0
  frequency.class[is.na(frequency.class[, 1]), 1] <- 0
  
  #mean of empty class returns NaN; replace by NA
  frequency.class[is.nan(frequency.class[, 2]), 2] <- NA
  
  if (standardize) {
    log_ratio <-
      (un_log_ratio - frequency.class[bins, 2]) / frequency.class[bins, 3]
    
    if (length(bins) > 1) {
      for (i in which(frequency.class[, 1] < 10)) {
        warning(
          paste(
            "The number of markers with allele frequencies in bin",
            rownames(frequency.class)[i],
            "is less than 10: you should probably increase bin width."
          )
        )
      }
    }
    
    if (is.na(p.side)) {
      pval <- 2 * pnorm(-abs(log_ratio))
    } else if (p.side == "left") {
      pval <- pnorm(log_ratio)
    } else if (p.side == "right") {
      pval <- 1 - pnorm(log_ratio)
    }
    
    log_pval <- -log10(p.adjust(pval, method = p.adjust.method))
    
    max_log_pval <- max(log_pval[log_pval != "Inf"])
    log_pval[log_pval == "Inf"] <- max_log_pval + 1
  }
  
  if (include_freq) {
    res <-
      scan[c("CHR", "POSITION", colnames(scan)[c(freq1_col_nr, freq2_col_nr)])]
  } else{
    res <- scan[c("CHR", "POSITION")]
  }
  
  if (standardize) {
    res[[final_score_colname]] <- log_ratio
    log_column_name <- "LOGPVALUE"
    if (!is.na(p.side)) {
      log_column_name <- paste0(log_column_name, "_", toupper(p.side))
    }
    res[[log_column_name]] <- log_pval
  } else{
    res[[paste0("UN", final_score_colname)]] <- un_log_ratio
  }
  
  colnames(frequency.class) <- c("N_MRK",
                                 "MEAN_UNIHS",
                                 "SD_UNIHS",
                                 "LOWER_QT",
                                 "UPPER_QT")
  
  return(list(ihs = res, frequency.class = frequency.class))
}
