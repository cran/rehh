#'Manhattan plot of iHS, XP-EHH or Rsb over a genome.
#'@description Manhattanplot of iHS, XP-EHH or Rsb over a genome.
#'@param data output of either \code{\link{ihh2ihs}}, \code{\link{ies2xpehh}} or \code{\link{ines2rsb}}.
#'@param pval logical. If \code{TRUE}, the p-value is plotted, otherwise the score itself.
#'@param threshold a horizontal line is added at the corresponding value(s),
#'for instance to represent a significance threshold.
#'A single value (upper or lower threshold) or two values (upper and lower) can be specified.
#'@param chr.name if \code{NA} (default), all chromosomes are plotted, otherwise only those specified.
#'@param cr highlight "candidate regions" specified by a data.frame with columns \code{CHR}, \code{START}
#'and \code{END} as obtained by the function \code{\link{calc_candidate_regions}}.
#'@param cr.col the color for highlighting
#'@param cr.opacity a value between 0 (invisible) and 1 (opaque).
#'@param cr.lab.cex text size of candidate region labels.
#'@param cr.lab.offset offset of candidate region labels.
#'@param cr.lab.pos if \code{"top"} (default) or \code{"bottom"}, candidate regions are labeled by numbers; to turn off, use \code{"none"}
#'@param mrk highlight marker specified by a data.frame containing the
#'colums \code{CHR} and \code{POSITION}. The row names of that data frame are taken as labels.
#'Alternatively a vector with marker IDs can be specified. In the latter case the ID is used as label.
#'@param mrk.cex size of marker label.
#'@param mrk.col color of the highlighted points.
#'@param mrk.pch type of the highlighted points.
#'@param mrk.lab.cex text size of marker label. If zero, no labels are printed.
#'@param mrk.lab.pos a position specifier for the text.
#'Values of 1, 2, 3 and 4, respectively indicate positions below, to the left of, above and to the right of the highlighted marker.
#'@param ignore_sign logical. If \code{TRUE}, absolute values are plotted.
#'@param cex size of the points representing markers in the plot(s) (see \code{\link[graphics]{par}}).
#'@param las orientation of axis labels (see \code{\link[graphics]{par}}).
#'@param pch type of the points representing markers in the plot(s) (see \code{\link[graphics]{points}}).
#'@param resolution Rasterize data points to the specified resolution and remove
#'duplicate points. Defaults to NULL, i.e. no rasterization. A typical value might be \code{c(1E5, 0.01)},
#'meaning that resolution on the x-axis (chromosomal position) is 100000 and on the y-axis (score or p-value) is 0.01.
#'@param ... further arguments to be passed to \code{\link[graphics]{plot.default}}.
#'@details The color of chromosomes is taken from the "Graphics Palette", see \code{\link[grDevices]{palette}}.
#'@details If a single chromosome is plotted, a genomic region can be specified by
#'argument \code{xlim}.
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
           mrk = NULL,
           mrk.cex = 1,
           mrk.col = "gray",
           mrk.pch = 1,
           mrk.lab.cex = 0.4,
           mrk.lab.pos = 4,
           ignore_sign = FALSE,
           cex = 0.5,
           las = 1,
           pch = 20,
           resolution = NULL,
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
    if (anyDuplicated(data[c("CHR", "POSITION")])) {
      stop("Data contains duplicated chromosomal positions.", call. = FALSE)
    }
    
    if (!anyNA(threshold)) {
      if (!is.numeric(threshold) | length(threshold) > 2) {
        stop("Threshold has to be given by one or two numbers.", call. = FALSE)
      }
    }
    
    if (!is.null(cr)) {
      if (!("CHR" %in% colnames(cr))) {
        stop("Region table does not contain a column named 'CHR'.", call. = FALSE)
      }
      if (!("START" %in% colnames(cr))) {
        stop("Region table does not contain a column named 'START'.",
             call. = FALSE)
      }
      if (!("END" %in% colnames(cr))) {
        stop("Region table does not contain a column named 'END'.", call. = FALSE)
      }
      
    }
    
    if (!is.null(mrk)) {
      if (is.data.frame(mrk)) {
        if (!("CHR" %in% colnames(mrk))) {
          stop("Marker table does not contain a column named 'CHR'.",
               call. = FALSE)
        }
        if (!("POSITION" %in% colnames(mrk))) {
          stop("Marker table does not contain a column named 'POSITION'.",
               call. = FALSE)
        }
      }
      else if (!(is.vector(mrk) & is.character(mrk))) {
        stop("Markers have to be provided as a table with positions or as a vector with IDs.")
      }
    }
    
    if (!is.null(resolution) &
        (length(resolution) != 2 |
         !is.numeric(resolution) | any(resolution <= 0))) {
      stop("Resolution has to be specified by a vector of two positive numbers.",
           call. = FALSE)
    }
    
    # perform plot
    
    ## try to identify statistic by column name
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
      ## check only on "VALUE" since write and read table may change column name
      col_nr <- which(grepl("VALUE", uppercolnames))
      if (length(col_nr) != 1) {
        stop("Could not determine column with p-values.", call. = FALSE)
      }
      score_colname <- colnames(data)[col_nr]
    }
    
    dot.args <- list(...)
    
    if (is.null(dot.args$ylab)) {
      if (pval) {
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
      } else{
        if (ignore_sign) {
          ylab <- bquote(italic("|" ~ .(statistic) ~ "|"))
        } else
          ylab <- bquote(italic(.(statistic)))
      }
      dot.args$ylab <- quote(ylab)
    }
    
    ## remove unused columns
    data <- data[c("CHR", "POSITION", score_colname)]
    
    ## refactor with factor level order as occurring in file
    data$CHR <- factor(data$CHR, levels = unique(data$CHR))
    
    ## take absolute values, if specified
    if (ignore_sign) {
      data[[score_colname]] <- abs(data[[score_colname]])
    }
    
    chromosomes <- as.character(unique(data$CHR))
    
    if (!is.character(chr.name)) {
      chr.name <- as.character(chr.name)
    }
    
    ## check if all highlighted markers can be found in complete data set
    if (!is.null(mrk)) {
      if (is.vector(mrk)) {
        nmrk <- length(mrk)
        data_highlighted <- data[mrk, ]
      } else{
        nmrk <- nrow(mrk)
        ## merge erases row.names; duplicate them as column
        mrk$Row.names <- row.names(mrk)
        
        data_highlighted <-
          suppressWarnings(merge(data, mrk, by = c("CHR", "POSITION")))
        
        ## set column back to row.names
        row.names(data_highlighted) <- data_highlighted$Row.names
        data_highlighted$Row.names <- NULL
      }
      
      ## remove rows with NAs (arising by empty subset)
      data_highlighted <-
        data_highlighted[!is.na(data_highlighted$CHR), ]
      
      if (nrow(data_highlighted) < nmrk) {
        warning(paste(
          "Could not find",
          nmrk - nrow(data_highlighted),
          "markers to be highlighted in data."
        ),
        call. = FALSE)
      }
      
      if (!is.null(resolution)) {
        # rasterize highlighted markers but do not remove duplicated
        data_highlighted$POSITION <-
          round(data_highlighted$POSITION / resolution[1]) * resolution[1]
        data_highlighted[[score_colname]] <-
          round(data_highlighted[[score_colname]] / resolution[2]) * resolution[2]
      }
    }
    
    if (!anyNA(chr.name)) {
      if (!all(chr.name %in% chromosomes)) {
        stop("Specified chromosomes not contained in data.", call. = FALSE)
      }
      chromosomes <- chr.name
      data <- data[data$CHR %in% chromosomes, ]
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
      if (is.null(dot.args$xlab)) {
        dot.args$xlab <- "Chromosome"
      }
      xaxt <- "n"
      dot.args$xlim <- NULL
    } else{
      p <- floor(log(chr_max[chromosomes], 1000))
      ## only shrink big scales, but never magnify small ones (p<0)
      scale <- 1000 ** max(0, p)
      ## no unit if p < 0
      unit <- c("", "(bp)", "(kb)", "(Mb)", "(Gb)")[max(-1, p)  + 2]
      
      if (is.null(dot.args$main)) {
        dot.args$main <- chromosomes[1]
      }
      if (is.null(dot.args$xlab)) {
        dot.args$xlab <- paste("Position", unit)
      }
      if (!is.null(dot.args$xlim)) {
        # subset to specified positions
        data <- data[data$POSITION >= dot.args$xlim[1] &
                       data$POSITION <= dot.args$xlim[2], ]
        
        dot.args$xlim <- dot.args$xlim / scale
      }
      xaxt <- "s"
    }
    
    if (!is.null(resolution)) {
      original_nrow <- nrow(data)
      
      data$POSITION <- round(data$POSITION / resolution[1])
      data[[score_colname]] <-
        round(data[[score_colname]] / resolution[2])
      
      data <- unique(data)
      
      data$POSITION <- data$POSITION * resolution[1]
      data[[score_colname]] <- data[[score_colname]] * resolution[2]
      
      if (nrow(data) < original_nrow) {
        cat("Rasterization reduced",
            original_nrow,
            "data points to",
            nrow(data),
            ".\n")
      }
    }
    
    do.call(plot,
            c(
              list(
                (data$POSITION + cum[as.character(data$CHR)]) / scale,
                data[[score_colname]],
                cex = cex,
                col = data$CHR,
                las = las,
                pch = pch,
                xaxt = xaxt
              ),
              dot.args
            ))
    
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
      cr <- cr[cr$CHR %in% chromosomes, ]
      
      if (nrow(cr) > 0) {
        col <-  adjustcolor(cr.col, alpha.f = cr.opacity)
        
        if (is.null(list(...)$ylim)) {
          ymin <- min(data[[score_colname]], na.rm = TRUE)
          ymax <- max(data[[score_colname]], na.rm = TRUE)
        } else{
          ymin <- list(...)$ylim[1]
          ymax <- list(...)$ylim[2]
        }
        
        xmin <- (cr$START + cum[as.character(cr$CHR)]) /
          scale
        
        xmax <- (cr$END + cum[as.character(cr$CHR)]) /
          scale
        
        rect(xmin,
             ymin,
             xmax,
             ymax,
             col = col,
             border = col)
        
        if ((cr.lab.pos == "top" |
             cr.lab.pos == "bottom") & cr.lab.cex != 0) {
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
    
    if (!is.null(mrk)) {
      if (nrow(data_highlighted) > 0) {
        points((data_highlighted$POSITION + cum[as.character(data_highlighted$CHR)]) / scale,
               data_highlighted[[score_colname]],
               cex = mrk.cex,
               col = mrk.col,
               pch = mrk.pch
        )
        
        if (mrk.lab.cex != 0) {
          text(
            (data_highlighted$POSITION + cum[as.character(data_highlighted$CHR)]) / scale,
            data_highlighted[[score_colname]],
            row.names(data_highlighted),
            pos = mrk.lab.pos,
            cex = mrk.lab.cex,
            xpd = TRUE
          )
        }
      }
    }
  }
