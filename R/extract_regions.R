#'Extract regions from a scan
#'@description Extract regions from a scan data frame.
#'@param scan A data frame with chromosomal positions like obtained
#'by \code{\link{scan_hh}}, \code{\link{ihh2ihs}}, \code{\link{ines2rsb}} or \code{\link{ies2xpehh}}.
#'@param regions A data frame with genomic regions like the output of \code{\link{calc_candidate_regions}}.
#'@param right logical, indicating if the intervals should be closed on the right (and open on the left) or vice versa.
#'@return A subset of data frame \code{scan}, retaining only positions belonging to
#'the regions specified in data frame \code{regions}.
#'@examples library(rehh.data)
#'data(wgscan.cgu)
#'regions <- data.frame(CHR = 12, START = 2.88e+7, END = 2.92e+7)
#'extract_regions(wgscan.cgu, regions)
#'@export
extract_regions <- function(scan, regions, right = TRUE) {
  # check parameters
  if (!is.data.frame(scan)) {
    stop("'scan' is not a data frame.", call. = FALSE)
  }
  if (is.null(scan$CHR)) {
    stop("'scan' has no column with name 'CHR'.", call. = FALSE)
  }
  if (is.null(scan$POSITION)) {
    stop("'scan' has no column with name 'POSITION'.", call. = FALSE)
  }
  if (!is.data.frame(regions)) {
    stop("'regions' is not a data frame.", call. = FALSE)
  }
  if (is.null(regions$CHR)) {
    stop("'regions' has no column with name 'CHR'.", call. = FALSE)
  }
  if (is.null(regions$START)) {
    stop("'regions' has no column with name 'START'.",
         call. = FALSE)
  }
  if (is.null(regions$END)) {
    stop("'regions' has no column with name 'END'.", call. = FALSE)
  }
  
  # perform subsetting
  
  # compare characters, not factors, but don't change 'scan'
  scan_CHR <- as.character(scan$CHR)
  
  select <- rep(FALSE, nrow(scan))
  
  for (i in seq_len(nrow(regions))) {
    select <- select | (scan_CHR == as.character(regions$CHR[i]) &
                          ((
                            right &
                              scan$POSITION > regions$START[i] &
                              scan$POSITION <= regions$END[i]
                          ) |
                            (
                              !right &
                                scan$POSITION >= regions$START[i] &
                                scan$POSITION < regions$END[i]
                            )
                          ))
  }
  return(scan[select, ])
}
