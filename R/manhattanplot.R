#'Manhattan plot of iHS, XP-EHH or Rsb over a genome.
#'@description Manhattanplot of iHS, XP-EHH or Rsb over a genome.
#'@param data output of either \code{\link{ihh2ihs}}, \code{\link{ies2xpehh}} or \code{\link{ines2rsb}}.
#'@param pval logical. If \code{TRUE}, the p-value is plotted, otherwise the score itself.
#'@param threshold a horizontal line is added at the corresponding value(s),
#'for instance to represent a significance threshold.
#'A single value (upper or lower threshold) or two values (upper and lower) can be specified.
#'@param chr.name if \code{NA} (default), all chromosomes are plotted, otherwise only those specified.
#'@param cr highlight "candidate regions" specified by a data.frame with three
#'columns: the first containing the chromosome, the other begin and end of the region as
#'obtained by the function \code{\link{calc_candidate_regions}}.
#'@param cr.col the color for highlighting
#'@param cr.opacity a value between 0 (invisible) and 1 (opaque).
#'@param cr.lab.cex text size of candidate region labels.
#'@param cr.lab.offset offset of candidate region labels.
#'@param cr.lab.pos if \code{"top"} (default) or \code{"bottom"}, candidate regions are labeled by numbers; to turn off, use \code{"none"}
#'@param main main title of the plot.
#'@param xlim set x coordinate range of the plot. Ignored, if more than one chromosome is depicted.
#'@param cex size of the points representing markers in the plot(s) (see \code{\link[graphics]{par}}).
#'@param las orientation of axis labels (see \code{\link[graphics]{par}}).
#'@param pch type of the points representing markers in the plot(s) (see \code{\link[graphics]{points}}).
#'@param ... further arguments to be passed to \code{\link[graphics]{plot.default}}.
#'@details The color of chromosomes is taken from the "Graphics Palette", see \code{\link[grDevices]{palette}}.
#'@return The function returns a plot.
#'@seealso \code{\link{ihh2ihs}}, \code{\link{ies2xpehh}}, \code{\link{ines2rsb}}, \code{\link{calc_candidate_regions}}.
#'@examples library(rehh.data)
#'data(wgscan.cgu)
#'## results from a genome scan (44,057 SNPs)
#'## see ?wgscan.eut and ?wgscan.cgu for details
#'wgscan.ihs <- ihh2ihs(wgscan.cgu)
#'manhattanplot(wgscan.ihs)
#'@export
#'@importFrom graphics axis rect
#'@importFrom grDevices adjustcolor
manhattanplot <-
  function(data,
           pval = FALSE,
           threshold = c(-2, 2),
           chr.name = NA,
           cr = NULL,
           cr.col = "gray",
           cr.opacity = 0.5,
           cr.lab.cex = 0.6,
           cr.lab.offset = 0,
           cr.lab.pos = "top",
           main = NA,
           xlim = NULL,
           cex = 0.5,
           las = 1,
           pch = 20,
           ...) {
    # check parameters
    
    ## output of ihh2ihs is a list; extract data frame
    if (is.list(data) & !is.data.frame(data)) {
      if (!is.null(data$iHS)) {
        data <- data$iHS
      } else if (!is.null(data$ihs)) {
        data <- data$ihs
      }
    }
    
    if (!("CHR" %in% colnames(data))) {
      stop("Data does not contain a column named 'CHR'.", call. = FALSE)
    }
    if (!("POSITION" %in% colnames(data))) {
      stop("Data does not contain a column named 'POSITION'.", call. = FALSE)
    }
    if (anyDuplicated(data[c("CHR" , "POSITION")])) {
      stop("Data contains duplicated chromosomal positions.", call. = FALSE)
    }
    
    if (!anyNA(threshold)) {
      if (!is.numeric(threshold) | length(threshold) > 2) {
        stop("Threshold has to be given by one or two numbers.", call. = FALSE)
      }
    }
    
    uppercolnames <- toupper(colnames(data))
    
    if ("IHS" %in% uppercolnames) {
      score_colname <- colnames(data)[which(uppercolnames == "IHS")]
      statistic <- "iHS"
    }  else if ("UNIHS" %in% uppercolnames) {
      score_colname <- colnames(data)[which(uppercolnames == "UNIHS")]
      statistic <- "unst. iHS"
    } else if ("XPEHH" %in% uppercolnames) {
      score_colname <- colnames(data)[which(uppercolnames == "XPEHH")]
      statistic <- "XP-EHH"
    } else if ("UNXPEHH" %in% uppercolnames) {
      score_colname <- colnames(data)[which(uppercolnames == "UNXPEHH")]
      statistic <- "unst. XP-EHH"
    } else if ("RSB" %in% uppercolnames) {
      score_colname <- colnames(data)[which(uppercolnames == "RSB")]
      statistic <- "Rsb"
    } else if ("UNRSB" %in% uppercolnames) {
      score_colname <- colnames(data)[which(uppercolnames == "UNRSB")]
      statistic <- "unst. Rsb"
    } else {
      ## best guess that score is in column 3
      score_colname <- colnames(data)[3]
      statistic <- score_colname
    }
    
    if (pval) {
      #check only on "VALUE" since write and read table may change column name
      col_nr <- which(grepl("VALUE", uppercolnames))
      if (length(col_nr) != 1) {
        stop("Could not determine column with p-values.", call. = FALSE)
      }
      
      if (grepl("LEFT", uppercolnames[col_nr])) {
        ylab <-
          bquote("-" ~ log[10] ~ "[" * Phi[scriptstyle(italic(.(statistic)))] *
                   "]")
      } else if (grepl("RIGHT", uppercolnames[col_nr])) {
        ylab <-
          bquote("-" ~ log[10] ~ "[1" ~ "-" ~ Phi[scriptstyle(italic(.(statistic)))] *
                   "]")
      } else {
        ylab <-
          bquote("-" ~ log[10] ~ "[2" * Phi[ ~ "-" ~ "|" ~ scriptstyle(italic(.(statistic))) ~
                                               "|"] * "]")
      }
      
      score_colname <- colnames(data)[col_nr]
    } else{
      ylab <- bquote(italic(.(statistic)))
    }
    
    #refactor with factor level order as occurring in file
    data$CHR <- factor(data$CHR, levels = unique(data$CHR))
    
    chromosomes <- as.character(unique(data$CHR))
    
    if (!is.character(chr.name)) {
      chr.name <- as.character(chr.name)
    }
    
    if (!anyNA(chr.name)) {
      if (!all(chr.name %in% chromosomes)) {
        stop("Specified chromosomes not contained in data.", call. = FALSE)
      }
      chromosomes <- chr.name
      data <- data[data$CHR %in% chromosomes,]
    }
    
    chr_max <-
      vapply(split(data, data$CHR, drop = TRUE), function(x) {
        max(x$POSITION)
      }, FUN.VALUE =  0)
    cum <- cumsum(c(0, chr_max[chromosomes]))
    label_pos <- (cum[-length(cum)] + cum[-1]) / 2
    
    cum <- cum[-length(cum)]
    names(cum) <- chromosomes
    names(label_pos) <- chromosomes
    
    if (length(chromosomes) > 1) {
      scale <- 1
      xlab <- "Chromosome"
      xaxt <- "n"
      xlim <- NULL
    } else{
      p <- floor(log(chr_max[chromosomes], 1000))
      ## only shrink big scales, but never magnify small ones (p<0)
      scale <- 1000 ** max(0, p)
      ## no unit if p < 0
      unit <- c("", "(bp)", "(kb)", "(Mb)", "(Gb)")[max(-1, p)  + 2]
      
      xlab <- paste("Position", unit)
      xaxt <- "s"
      if (!is.null(xlim)) {
        xlim <- xlim / scale
      }
    }
    
    plot(
      (data$POSITION + cum[as.character(data$CHR)]) / scale,
      data[[score_colname]],
      main = main,
      xlim = xlim,
      xlab = xlab,
      ylab = ylab,
      cex = cex,
      col = data$CHR,
      las = las,
      pch = pch,
      xaxt = xaxt,
      ...
    )
    
    if (length(chromosomes) > 1) {
      axis(1,
           at = label_pos[chromosomes],
           labels = chromosomes,
           las = las)
    }
    
    if (!anyNA(threshold)) {
      abline(h = threshold, lty = 2)
    }
    
    if (!is.null(cr)) {
      cr <- cr[cr[[1]] %in% chromosomes,]
      
      if (nrow(cr) > 0) {
        col <-  adjustcolor(cr.col, alpha.f = cr.opacity)
        
        if (is.null(list(...)$ylim)) {
          ymin <- min(data[[score_colname]], na.rm = TRUE)
          ymax <- max(data[[score_colname]], na.rm = TRUE)
        } else{
          ymin <- list(...)$ylim[1]
          ymax <- list(...)$ylim[2]
        }
        
        xmin <- (cr[[2]] + cum[as.character(cr[[1]])]) /
          scale
        
        xmax <- (cr[[3]] + cum[as.character(cr[[1]])]) /
          scale
        
        rect(xmin,
             ymin,
             xmax,
             ymax,
             col = col,
             border = col)
        
        if (cr.lab.pos == "top" |
            cr.lab.pos == "bottom") {
          text((xmin + xmax) / 2,
               ifelse(cr.lab.pos == "top", ymax, ymin),
               rownames(cr),
               pos = ifelse(cr.lab.pos == "top", 3, 1),
               xpd = TRUE,
               cex = cr.lab.cex,
               offset = cr.lab.offset
          )
        }
      }
    }
  }
