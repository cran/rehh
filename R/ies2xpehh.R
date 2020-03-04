#'Compute XP-EHH
#'@description Compute XP-EHH (standardized ratio of iES of two populations).
#'@param scan_pop1 a data frame with markers in rows and columns with chromosome name, position of the
#'marker, frequency of the ancestral allele and iES as obtained by \code{\link{scan_hh}} on the first population.
#'@param scan_pop2 a data frame with markers in rows and columns with chromosome name, position of the
#'marker, frequency of the ancestral allele and iES as obtained by \code{\link{scan_hh}} on the second population.
#'@param popname1 short ID/name of the first population; to be added to an output column name.
#'@param popname2 short ID/name of the second population; to be added to an output column name.
#'@param min_nhaplo discard positions where in at least one of the populations fewer than \code{min_nhaplo} haplotypes
#'have been evaluated (default \code{NA}).
#'@param standardize logical. If \code{TRUE} (default), then standardize XP-EHH, else report unstandardized XP-EHH.
#'@param include_freq logical. If \code{TRUE} include columns with allele frequencies into result.
#'@param p.side side to which refers the p-value. Default \code{NA}, meaning two-sided. Can be set
#'to \code{"left"} or \code{"right"}.
#'@param p.adjust.method method passed to function \code{\link[stats]{p.adjust}} to correct the p-value for
#' multiple testing. Default \code{"none"}.
#'@param verbose logical. If \code{TRUE} (default), report number of markers of the two source data frames and result data frame.
#'@details Log ratio of iES (population 1 over population 2) computed as described in Sabeti et al. (2007).

#'Note that the two data frames are merged on the basis of chromosome and position. Marker names
#'are kept, if they are identical and unique in both data frames.
#'
#'Since the standardized XP-EHH values follow, if markers evolve predominantly neutrally, approximately
#'a standard Gaussian distribution, it is practical to assign to the values a p-value relative
#'to the null-hypothesis of neutral evolution. The parameter \code{p.side} determines
#'if the p-value is assigned to both sides of the distribution or to one side of interest.
#'@return The returned value is a data frame with markers in rows and columns for chromosome name, marker position,
#'XP-EHH and, if standardized, p-value in a negative log10 scale. Optionally, allele frequencies are included.
#'@references Gautier, M. and Naves, M. (2011). Footprints of selection in the ancestral admixture of a New World Creole cattle breed. \emph{Molecular Ecology}, \strong{20}, 3128-3143.
#'
#'Sabeti, P.C. et al. (2007). Genome-wide detection and characterization of positive selection in human populations. \emph{Nature}, \strong{449}, 913-918.
#'@seealso \code{\link{scan_hh}}, \code{\link{distribplot}}, \code{\link{manhattanplot}}
#'@examples library(rehh.data)
#'data(wgscan.cgu) ; data(wgscan.eut)
#'## results from a genome scan (44,057 SNPs)
#'##see ?wgscan.eut and ?wgscan.cgu for details
#'wgscan.xpehh <- ies2xpehh(wgscan.cgu, wgscan.eut, "CGU", "EUT")
#'@export
#'@importFrom stats sd p.adjust p.adjust.methods
ies2xpehh <-
  function(scan_pop1,
           scan_pop2,
           popname1 = NA,
           popname2 = NA,
           min_nhaplo = NA,
           standardize = TRUE,
           include_freq = FALSE,
           p.side = NA,
           p.adjust.method = "none",
           verbose = TRUE) {
    ## check parameters
    initial_score_colname <- "IES"
    old_initial_score_colname <- "iES_Sabeti_et_al_2007"
    final_score_colname <- "XPEHH"
    
    if (is.null(scan_pop1$CHR)) {
      stop("Data for population 1 does not contain a column named 'CHR'.",
           call. = FALSE)
    }
    if (is.null(scan_pop2$CHR)) {
      stop("Data for population 2 does not contain a column named 'CHR'.",
           call. = FALSE)
    }
    if (is.null(scan_pop1$POSITION)) {
      stop("Data for population 1 does not contain a column named 'POSITION'.",
           call. = FALSE)
    }
    if (is.null(scan_pop2$POSITION)) {
      stop("Data for population 2 does not contain a column named 'POSITION'.",
           call. = FALSE)
    }
    if (anyDuplicated(scan_pop1[c("CHR", "POSITION")])) {
      stop("Data for population 1 contains duplicated chromosomal positions.",
           call. = FALSE)
    }
    if (anyDuplicated(scan_pop2[c("CHR", "POSITION")])) {
      stop("Data for population 2 contains duplicated chromosomal positions.",
           call. = FALSE)
    }
    if (!is.na(min_nhaplo) &
        (round(min_nhaplo) != min_nhaplo | min_nhaplo < 2)) {
      stop("min_nhaplo must be an integer number greater than 1.")
    }
    
    if (is.null(scan_pop1[[initial_score_colname]])) {
      #if column with old name exists, rename to new name
      if (!is.null(scan_pop1[[old_initial_score_colname]])) {
        colnames(scan_pop1)[which(colnames(scan_pop1) == old_initial_score_colname)] <-
          initial_score_colname
      } else{
        stop("Data for population 1 does not contain a column with iES values.",
             call. = FALSE)
      }
    }
    if (is.null(scan_pop2[[initial_score_colname]])) {
      #if column with old name exists, rename to new name
      if (!is.null(scan_pop2[[old_initial_score_colname]])) {
        colnames(scan_pop2)[which(colnames(scan_pop2) == old_initial_score_colname)] <-
          initial_score_colname
      } else{
        stop("Data for population 2 does not contain a column with iES values.",
             call. = FALSE)
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
    
    ## perform calculation
    if (is.na(popname1) |
        is.na(popname2)) {
      suffixes <- c(".POP1", ".POP2")
    } else{
      final_score_colname <-
        paste0(final_score_colname, "_", popname1, "_", popname2)
      suffixes <- c(paste0(".", popname1), paste0(".", popname2))
    }
    
    if (include_freq) {
      freq_colnames1 <- intersect(colnames(scan_pop1),
                                  c("FREQ_A", "FREQ_D", "FREQ_MAJ", "FREQ_MIN"))
      colnames1 <-
        c("CHR", "POSITION", freq_colnames1, initial_score_colname)
      freq_colnames2 <- intersect(colnames(scan_pop2),
                                  c("FREQ_A", "FREQ_D", "FREQ_MAJ", "FREQ_MIN"))
      colnames2 <-
        c("CHR", "POSITION", freq_colnames2, initial_score_colname)
    } else{
      colnames1 <-
        colnames2 <- c("CHR", "POSITION", initial_score_colname)
    }
    
    if (!is.na(min_nhaplo)) {
      nhaplo1_col_nr <- which(colnames(scan_pop1) == "NHAPLO_A")
      nhaplo2_col_nr <- which(colnames(scan_pop1) == "NHAPLO_D")
      if (length(nhaplo1_col_nr) != 1 |
          length(nhaplo2_col_nr) != 1) {
        nhaplo1_col_nr <- which(colnames(scan_pop1) == "NHAPLO_MAJ")
        nhaplo2_col_nr <- which(colnames(scan_pop1) == "NHAPLO_MIN")
        if (length(nhaplo1_col_nr) != 1 |
            length(nhaplo2_col_nr) != 1) {
          stop("Column(s) with number of haplotypes not found in population 1.",
               call. = FALSE)
        }
      }
      
      if (verbose)
        cat(
          "Discard focal markers with less than",
          min_nhaplo,
          "evaluated haplotypes in population 1.\n"
        )
      
      mrk_sel <-
        scan_pop1[[nhaplo1_col_nr]] + scan_pop1[[nhaplo2_col_nr]] >= min_nhaplo
      
      if (sum(mrk_sel) == nrow(scan_pop1)) {
        if (verbose)
          cat("No marker discarded.\n")
      } else{
        if (verbose)
          cat(nrow(scan_pop1) - sum(mrk_sel), "markers discarded.\n")
        
        scan_pop1 <- scan_pop1[mrk_sel, ]
        
        if (verbose)
          cat(nrow(scan_pop1), "markers remaining.\n")
        
      }
      
      nhaplo1_col_nr <- which(colnames(scan_pop2) == "NHAPLO_A")
      nhaplo2_col_nr <- which(colnames(scan_pop2) == "NHAPLO_D")
      if (length(nhaplo1_col_nr) != 1 |
          length(nhaplo2_col_nr) != 1) {
        nhaplo1_col_nr <- which(colnames(scan_pop2) == "NHAPLO_MAJ")
        nhaplo2_col_nr <- which(colnames(scan_pop2) == "NHAPLO_MIN")
        if (length(nhaplo1_col_nr) != 1 |
            length(nhaplo2_col_nr) != 1) {
          stop("Column(s) with number of haplotypes not found in population 2.",
               call. = FALSE)
        }
      }
      if (verbose)
        cat(
          "Discard focal markers with less than",
          min_nhaplo,
          "evaluated haplotypes in population 2.\n"
        )
      
      mrk_sel <-
        scan_pop2[[nhaplo1_col_nr]] + scan_pop2[[nhaplo2_col_nr]] >= min_nhaplo
      
      if (sum(mrk_sel) == nrow(scan_pop2)) {
        if (verbose)
          cat("No marker discarded.\n")
      } else{
        if (verbose)
          cat(nrow(scan_pop2) - sum(mrk_sel), "markers discarded.\n")
        
        scan_pop2 <- scan_pop2[mrk_sel, ]
        
        if (verbose)
          cat(nrow(scan_pop2), "markers remaining.\n")
        
      }
      
    }
    
    if (verbose) {
      cat("Scan of pop1 contains", nrow(scan_pop1), "markers.\n")
      cat("Scan of pop2 contains", nrow(scan_pop2), "markers.\n")
    }
    
    # work-around to keep row names: set them as columns
    scan_pop1$rownames <- row.names(scan_pop1)
    colnames1 <- c("rownames", colnames1)
    scan_pop2$rownames <- row.names(scan_pop2)
    colnames2 <- c("rownames", colnames2)
    
    res <- merge(
      scan_pop1[, colnames1],
      scan_pop2[, colnames2],
      by = c("CHR", "POSITION"),
      suffixes = suffixes,
      sort = FALSE
    )
    
    if (nrow(res) == 0) {
      stop("Merged data contains no marker.",
           call. = FALSE)
    }
    
    # remove big objects
    rm(scan_pop1)
    rm(scan_pop2)
    
    # get row names right
    col_rn1 <- paste0("rownames", suffixes[1])
    col_rn2 <- paste0("rownames", suffixes[2])
    
    # if merged rownames are identical and unique ...
    if (identical(res[[col_rn1]], res[[col_rn2]]) &
        !anyDuplicated(res[[col_rn1]])) {
      row.names(res) <- res[[col_rn1]]
    }
    res[[col_rn1]] <- res[[col_rn2]] <- NULL
    
    if (verbose)
      cat("Merged data contains",
          nrow(res),
          "markers.\n")
    
    un_log_ratio <-
      log(res[[paste0(initial_score_colname, suffixes[1])]] / res[[paste0(initial_score_colname, suffixes[2])]])
    
    # replace infinity by NA
    un_log_ratio[is.infinite(un_log_ratio)] <- NA
    
    # remove iES values
    res[paste0(initial_score_colname, suffixes)] <- NULL
    
    if (standardize) {
      mean <- mean(un_log_ratio, na.rm = TRUE)
      sd <- sd(un_log_ratio, na.rm = TRUE)
      log_ratio <- (un_log_ratio - mean) / sd
      
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
    
    # add final score and, if standardized, p-values
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
    
    return(res)
  }
