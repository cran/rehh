#'Plot EHH around a focal marker
#'@description Plot curve of EHH values around a focal marker.
#'@param x data (output of \code{\link{calc_ehh}}).
#'@param ylim the y limits of the plot
#'@param type plot type (see \code{\link[graphics]{matplot}}).
#'@param main title for the plot (default \code{NA}, i.e. none).
#'@param xlab title for the x-axis.
#'@param ylab title for the y-axis.
#'@param bty box type around plot (see \code{\link[graphics]{par}}).
#'@param col color for the ancestral and derived alleles (respectively) curves.
#'@param mrk.col color of the vertical line at the focal marker position.
#'@param lty line type for the ancestral and derived allele EHH (respectively) curves.
#'@param legend legend text.
#'@param legend.xy.coords if \code{"automatic"} (default) places legend either top left or top right;
#'if \code{"none"}, no legend is drawn; otherwise the argument is passed to \code{\link[graphics]{legend}}.
#'@param ... further arguments to be passed to function \code{\link[graphics]{matplot}}.
#'@seealso \code{\link{data2haplohh}}, \code{\link{calc_ehh}}, \code{\link{plot.ehhs}}, \code{\link{scan_hh}}.
#'@examples
#'#example haplohh object (280 haplotypes, 1424 SNPs)
#'#see ?haplohh_cgu_bta12 for details
#'data(haplohh_cgu_bta12)
#'#computing EHH statistics for the marker "F1205400"
#'#which displays a strong signal of selection
#'ehh <- calc_ehh(haplohh_cgu_bta12, mrk = "F1205400")
#'plot(ehh)
#'@export
#'@importFrom graphics matplot
plot.ehh <-
  function(x,
           ylim = c(0, 1),
           type = "l",
           main = paste0("EHH around '", x$mrk.name, "'"),
           xlab = "Position",
           ylab = "Extended Haplotype Homozygosity",
           col = c("blue", "red", "violet", "orange"),
           mrk.col = "gray",
           bty = "n",
           lty = 1,
           legend = NA,
           legend.xy.coords = "automatic",
           ...) {
    # check parameters
    if (is.null(x$ehh)) {
      stop("Data does not contain a data frame with EHH values.", call. = FALSE)
    }
    if (is.null(x$ehh$POSITION)) {
      stop("Data frame seems not to contain a column with marker positions.",
           call. = FALSE)
    }
    if (nrow(x$ehh) == 0) {
      stop("Empty data frame.", call. = FALSE)
    }
    if (is.null(x$mrk.name)) {
      stop("No marker id found.", call. = FALSE)
    }
    foc_pos <- x$ehh[x$mrk.name, "POSITION"]
    if (is.na(foc_pos)) {
      stop(paste0("No position information found for marker '", x$mrk.name, "'."),
           call. = FALSE)
    }
    ehh <- x$ehh[grepl("EHH", colnames(x$ehh))]
    if (length(ehh) == 0) {
      stop("Data frame seems not to contain a column with EHH values.",
           call. = FALSE)
    }
    
    description <- get_description(colnames(ehh))
    description_colors <- get_description_colors(description, col)
    
    # perform plot
    p <- floor(log(max(x$ehh$POSITION), 1000))
    ## only shrink big scales, but never magnify small ones (p<0)
    scale <- 1000 ** max(0, p)
    ## no unit if p < 0
    unit <- c("", "(bp)", "(kb)", "(Mb)", "(Gb)")[max(-1, p)  + 2]
    
    dot_args <- list(...)
    if (!is.null(dot_args$xlim)) {
      dot_args$xlim <- dot_args$xlim / scale
    }
    
    do.call("matplot", c(
      list(
        x = x$ehh$POSITION / scale,
        y = ehh,
        ylim = ylim,
        type = type,
        main = main,
        xlab = paste(xlab, unit),
        ylab = ylab,
        bty = bty,
        col = description_colors,
        lty = lty
      ),
      dot_args
    ))
    
    abline(v = foc_pos / scale,
           lty = 2,
           col = mrk.col)
    
    # no argument legend, but argument legend.xy.coords
    if (is.na(legend) & !anyNA(legend.xy.coords)) {
      if (is.numeric(legend.xy.coords)) {
        legend.xy.coords[1] <- legend.xy.coords[1] / scale
      }
      if (length(legend.xy.coords) == 2) {
        legend(
          legend.xy.coords[1],
          legend.xy.coords[2],
          legend = description,
          bty = bty,
          col = description_colors,
          lty = lty,
          xpd = TRUE
        )
      } else{
        if (legend.xy.coords != "none") {
          if (legend.xy.coords == "automatic") {
            if (is.null(dot_args$xlim)) {
              picturemiddle <- mean(range(x$ehh$POSITION) / scale)
            } else{
              picturemiddle <- mean(dot_args$xlim)
            }
            legend.xy.coords <-
              ifelse(foc_pos / scale > picturemiddle, "topleft", "topright")
          }
          legend(
            legend.xy.coords,
            legend = description,
            bty = bty,
            col = description_colors,
            lty = lty,
            xpd = TRUE
          )
        }
      }
    }
  }

# translate column names to descriptions
get_description <- function(colnames) {
  if (length(colnames) > 2) {
    index_derived <- seq_len(length(colnames) - 1)
  } else{
    index_derived <- ""
  }
  
  if ("EHH_A" %in% colnames) {
    return(c("Ancestral", paste0("Derived", index_derived)))
  }
  else{
    return(c("Major", paste0("Minor", index_derived)))
  }
}

# get the associated color for each meaning, e.g. blue for ancestral, etc
get_description_colors <- function(description, colors) {
  # repeat color vector if necessary
  if (length(colors) < length(description)) {
    colors <- rep_len(colors, length(description))
  }
  
  vapply(description, function(x) {
    x <- toupper(x)
    
    if (x == "ANCESTRAL" | x == "MAJOR") {
      return(colors[1])
    } else if (x == "DERIVED" | x == "MINOR") {
      return(colors[2])
    } else if (substring(x, 1, 7) == "DERIVED" |
               substring(x, 1, 5) == "MINOR") {
      other_nr <-
        suppressWarnings(as.integer(regmatches(x, gregexpr("[0-9]+$", x))[[1]]))
      
      if (other_nr > 0 & other_nr < length(colors)) {
        return(colors[other_nr + 1L])
      }
    }
    # default, if something went wrong
    return("gray")
  }, USE.NAMES = FALSE, FUN.VALUE = "gray")
}
