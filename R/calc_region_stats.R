#'Calculate score statistics for given regions
#'@description Calculate score statistics (extremal values) for given regions. This function
#'is intended for the comparison of different scores for the same chromosomal regions.
#'@param scan a data frame containing scores (output of \code{\link{ihh2ihs}}, \code{\link{ines2rsb}} or \code{\link{ies2xpehh}}).
#'@param regions a data frame with column names \code{CHR}, \code{START} and \code{END},
#'specifying chromosomal regions (e.g. as obtained by function \code{\link{calc_candidate_regions}}).
#'@param threshold boundary score above which markers are defined as "extreme".
#'@param pval logical. If \code{TRUE} use the (negative log-) p-value instead of the score.
#'@param ignore_sign logical. If \code{TRUE} (default), take absolute values of score.
#'@param right logical, indicating if the regions should be closed on the right (and open on the left) or vice versa.
#'@return A data frame with chromosomal regions.
#'For each region the overall number of markers, their mean and maximum,
#'the number of markers with extremal values, their percentage of all markers
#'and their average are reported.
#'@seealso \code{\link{calc_candidate_regions}}
#'@export
calc_region_stats <- function(scan,
                              regions,
                              threshold = NA,
                              pval = FALSE,
                              ignore_sign = FALSE,
                              right = TRUE) {
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
  if (!("CHR" %in% colnames(regions))) {
    stop("Regions do not contain a column named 'CHR'.", call. = FALSE)
  }
  if (!("START" %in% colnames(regions))) {
    stop("Regions do not contain a column named 'START'.", call. = FALSE)
  }
  if (!("END" %in% colnames(regions))) {
    stop("Regions do not contain a column named 'END'.", call. = FALSE)
  }
  if (is.na(threshold) | !is.numeric(threshold)) {
    stop("Specify a numeric value for 'threshold'.",
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
    warning("Could not identify column with score/p-value. Trying third column.")
    score_column_nr <- 3
  }
  
  if (ignore_sign) {
    score <- abs(scan[[score_column_nr]])
  } else{
    score <- scan[[score_column_nr]]
  }
  
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
  
  #create empty data frames for integer and numeric values
  region_stats_num <-
    data.frame(matrix(ncol = 4, nrow = length(regions$CHR)))
  region_stats_int <-
    data.frame(matrix(ncol = 2, nrow = length(regions$CHR)))
  
  # cannot be vectorized by "apply" because vectors don't allow different data types
  for (i in seq_along(regions$CHR)) {
    # subset score to chromosomal region
    if (right) {
      # left open, right closed intervals
      region_score <- score[scan$CHR == regions[["CHR"]][i] &
                              regions[["START"]][i] < scan$POSITION &
                              scan$POSITION <= regions[["END"]][i]]
    } else{
      # left closed, right open intervals
      region_score <- score[scan$CHR == regions[["CHR"]][i] &
                              regions[["START"]][i] <= scan$POSITION &
                              scan$POSITION < regions[["END"]][i]]
    }
    # subset of extremal markers in region
    region_score_extr <- region_score[region_score >= threshold]
    
    # calculate integer values
    region_stats_int[i, ] <- c(length(region_score),
                               length(region_score_extr))
    # calculate numeric values
    if (length(region_score) > 0) {
      perc <- length(region_score_extr) / length(region_score) * 100
      mean_mrk <- mean(region_score)
      max_mrk <- max(region_score)
    } else{
      perc <- NA
      mean_mrk <- NA
      max_mrk <- NA
    }
    
    if (length(region_score_extr) > 0) {
      mean_mrk_extr <- mean(region_score_extr)
    } else{
      mean_mrk_extr <- NA
    }
    
    region_stats_num[i, ] <-
      c(mean_mrk, max_mrk, round(perc, 2), mean_mrk_extr)
  }
  region_stats <-
    cbind(
      regions[, c("CHR", "START", "END")],
      region_stats_int[, 1],
      region_stats_num[, c(1, 2)],
      region_stats_int[, 2],
      region_stats_num[, c(3, 4)]
    )
  
  colnames(region_stats) <- colnames
  
  ## output
  return(region_stats)
}
