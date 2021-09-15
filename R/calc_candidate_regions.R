#'Determine candidate regions of selection
#'@description Determine candidate regions of selection.
#'@param scan a data frame containing scores (output of \code{\link{ihh2ihs}}, \code{\link{ines2rsb}} or \code{\link{ies2xpehh}}).
#'@param threshold boundary score above which markers are defined as "extreme".
#'@param pval logical. If \code{TRUE} use the (negative log-) p-value instead of the score.
#'@param ignore_sign logical. If \code{TRUE} (default), take absolute values of score.
#'@param window_size size of sliding windows. If set to 1, no windows
#'are constructed and only the individual extremal markers are reported.
#'@param overlap size of window overlap (default 0, i.e. no overlap).
#'@param right logical, indicating if the windows should be closed on the right (and open on the left) or vice versa.
#'@param min_n_mrk minimum number of markers per window.
#'@param min_n_extr_mrk minimum number of markers with extreme
#'value in a window.
#'@param min_perc_extr_mrk minimum percentage of extremal markers among all markers.
#'@param join_neighbors logical. If \code{TRUE} (default), merge neighboring windows with
#'extreme values.
#'@details There is no generally agreed method how to determine genomic
#'regions which might have been under recent selection. Since selection tends
#'to yield clusters of markers with outlier values, a common approach is
#'to search for regions with an elevated number or fraction
#'of outlier or extremal markers.
#'This function allows to set three conditions a window must fulfill
#'in order to classify as candidate region:
#'\itemize{
#'\item \code{min_n_mrk} a minimum number of (any) markers.
#'\item \code{min_n_extr_mrk} a minimum number of markers with outlier / extreme value.
#'\item \code{min_perc_extr_mrk} a minimum percentage of extremal markers among all markers.
#'}
#'"Extreme" markers are defined by having a score above the specified \code{threshold}.
#'@return A data frame with chromosomal regions, i.e. windows that fulfill
#'the necessary conditions to qualify as candidate regions under selection.
#'For each region the overall number of markers, their mean and maximum,
#'the number of markers with extremal values, their percentage of all markers
#'and their average are reported.
#'@seealso \code{\link{calc_region_stats}}
#'@export
calc_candidate_regions <- function(scan,
                                   threshold = NA,
                                   pval = FALSE,
                                   ignore_sign = FALSE,
                                   window_size = 1000000,
                                   overlap = 0,
                                   right = TRUE,
                                   min_n_mrk = 1,
                                   min_n_extr_mrk = 1,
                                   min_perc_extr_mrk = 0,
                                   join_neighbors = TRUE) {
  ## check parameters
  
  if (!is.null(scan$iHS)) {
    scan <- scan$iHS
  } else if (!is.null(scan$ihs)) {
    scan <- scan$ihs
  }
  
  if (!("CHR" %in% colnames(scan))) {
    stop("Scan does not contain a column named 'CHR'.", call. = FALSE)
  }
  if (!("POSITION" %in% colnames(scan))) {
    stop("Scan does not contain a column named 'POSITION'.", call. = FALSE)
  }
  if (is.na(threshold) | !is.numeric(threshold)) {
    stop("Specify a numeric value for 'threshold'.",
         call. = FALSE)
  }
  # need integer numbers, otherwise "%%" causes trouble with i386
  window_size <- as.integer(window_size)
  overlap <- as.integer(overlap)
  min_n_mrk <- as.integer(min_n_mrk)
  min_n_extr_mrk <- as.integer(min_n_extr_mrk)
  
  if (window_size < 1 | window_size %% 1L != 0L) {
    stop("Window size has to be a positive integer number.", call. = FALSE)
  }
  if (is.na(overlap) |
      overlap < 0 |
      overlap >= window_size |
      overlap %% 1L != 0L |
      (overlap != 0L & window_size %% overlap != 0L)) {
    stop("'overlap' has to be zero or an integer factor of 'window_size'.",
         call. = FALSE)
  }
  if (min_n_mrk < 1 | min_n_mrk %% 1 != 0L) {
    stop("'min_n_mrk' has to be a positive integer number.", call. = FALSE)
  }
  if (min_n_extr_mrk < 0 | min_n_extr_mrk %% 1L != 0L) {
    stop("'min_n_extr_mrk' has to be a positive integer number.",
         call. = FALSE)
  }
  if (min_perc_extr_mrk < 0 | min_perc_extr_mrk > 100) {
    stop("'min_perc_extr_mrk' has to be a number between 0 and 100.",
         call. = FALSE)
  }
  
  if (pval) {
    matched_columns <-
      grep("PVALUE", colnames(scan), ignore.case = TRUE)
  } else{
    matched_columns <-
      grep("IHS|RSB|XPEHH", colnames(scan), ignore.case = TRUE)
  }
  if (length(matched_columns) == 1) {
    score_column_nr <- matched_columns[1]
  } else{
    # guess that it is the third column
    warning("Could not identify data column. Trying third column.")
    score_column_nr <- 3
  }
  
  # perform calculation
  ## remove NA
  scan <- scan[!is.na(scan[[score_column_nr]]), ]
  
  if (ignore_sign) {
    score <- abs(scan[[score_column_nr]])
  } else{
    score <- scan[[score_column_nr]]
  }
  
  if (window_size > 1) {
    colnames <- c(
      "CHR",
      "START",
      "END",
      "N_MRK",
      "MEAN_MRK",
      "MAX_MRK",
      "N_EXTR_MRK",
      "PERC_EXTR_MRK",
      "MEAN_EXTR_MRK"
    )
    
    #create empty data frame
    cand_reg <- data.frame(matrix(ncol = 9, nrow = 0))
    
    chromosomes <- unique(scan$CHR)
    
    for (chr in chromosomes) {
      chr_pos <- scan[scan$CHR == chr, "POSITION"]
      chr_score <- score[scan$CHR == chr]
      
      chr_cand_reg <- data.frame(matrix(ncol = 6, nrow = 0))
      
      if (overlap != 0) {
        offsets <- seq(0, window_size - overlap, overlap)
      } else{
        offsets <- 0
      }
      
      for (offset in offsets) {
        breaks <- seq(
          floor(min(chr_pos) / window_size) * window_size + offset,
          max(chr_pos) + window_size - 1,
          window_size
        )
        # empty window (may arise if only a single value in chr_pos)
        if (length(breaks) < 2)
          next
        
        windows <- cut(chr_pos,
                       breaks = breaks,
                       right = right,
                       labels = FALSE)
        
        n_mrk <- tapply(chr_score, windows, function(x) {
          length(x)
        })
        
        mean_mrk <- tapply(chr_score, windows, function(x) {
          mean(x)
        })
        
        max_mrk <- tapply(chr_score, windows, function(x) {
          max(x)
        })
        
        n_extr_mrk <- tapply(chr_score, windows, function(x) {
          sum(x >= threshold)
        })
        
        mean_extr_mrk <- tapply(chr_score, windows, function(x) {
          mean(x[x >= threshold])
        })
        
        perc_extr_mrk <- n_extr_mrk / n_mrk * 100
        
        condition <- n_mrk >= min_n_mrk &
          n_extr_mrk >= min_n_extr_mrk &
          perc_extr_mrk >= min_perc_extr_mrk
        
        # subset all vectors by condition
        n_mrk <- n_mrk[condition]
        mean_mrk <- mean_mrk[condition]
        max_mrk <- max_mrk[condition]
        n_extr_mrk <- n_extr_mrk[condition]
        perc_extr_mrk <- perc_extr_mrk[condition]
        mean_extr_mrk <- mean_extr_mrk[condition]
        
        cand_window_nr <- as.integer(names(n_extr_mrk))
        
        if (length(cand_window_nr) > 0) {
          # append chromosomal candidates for each window offset
          chr_cand_reg <- rbind(
            chr_cand_reg,
            data.frame(
              chr,
              breaks[cand_window_nr],
              breaks[cand_window_nr + 1],
              n_mrk,
              mean_mrk,
              max_mrk,
              n_extr_mrk,
              round(perc_extr_mrk, 2),
              mean_extr_mrk
            )
          )
        }
      }
      
      if (nrow(chr_cand_reg) > 1) {
        # sort different offset windows by position
        chr_cand_reg <- chr_cand_reg[order(chr_cand_reg[[2]]), ]
        
        # join neighboring windows
        if (join_neighbors) {
          i <- 1 # index of row
          # no for-loop since nrows may decrease
          while (i < nrow(chr_cand_reg)) {
            j <- 0 # number of overlapping windows to the right
            # increase j until there is no overlap between
            # consecutive windows any more
            while (i + j < nrow(chr_cand_reg) &
                   chr_cand_reg[i + j, 3] >= chr_cand_reg[i + j + 1, 2]) {
              j <- j + 1
            }
            # if there are overlapping windows ...
            if (j > 0) {
              #update leftmost neighbor
              chr_cand_reg[i, 3] <- chr_cand_reg[i + j, 3]
              index <- which(chr_pos > chr_cand_reg[i, 2] &
                               chr_pos <= chr_cand_reg[i + j, 3])
              chr_cand_reg[i, 4] <- length(index)
              chr_cand_reg[i, 5] <- mean(chr_score[index])
              chr_cand_reg[i, 6] <-
                max(chr_score[index])
              chr_cand_reg[i, 7] <-
                sum(chr_score[index] >= threshold)
              chr_cand_reg[i, 8] <- round(chr_cand_reg[i, 7] /
                                            chr_cand_reg[i, 4] * 100, 2)
              chr_cand_reg[i, 9] <-
                mean(chr_score[index][chr_score[index] >=
                                        threshold])
              #eliminate all right neighbors
              chr_cand_reg <- chr_cand_reg[-((i + 1):(i +
                                                        j)), ]
            }
            i <- i + 1
          }
        }
      }
      cand_reg <- rbind(cand_reg, chr_cand_reg)
    }
    
    colnames(cand_reg) <- colnames
    rownames(cand_reg) <- seq_len(nrow(cand_reg))
    
  } else{
    # window_size = 1 => report subset of scan with individual extremal markers
    cand_reg <-
      scan[score >= threshold,
           c("CHR", "POSITION", "POSITION", colnames(scan)[score_column_nr])]
    
    colnames(cand_reg) <- c("CHR",
                            "START",
                            "END",
                            "EXTR_MRK")
  }
  
  ## output
  
  return(cand_reg)
}
