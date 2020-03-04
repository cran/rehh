#'Plots furcation trees around a focal marker
#'
#'@param x an object of class furcation (see \code{\link{calc_furcation}}).
#'@param allele If \code{NA} (default), furcation trees for all alleles of the focal marker are plotted,
#'otherwise for the specified alleles. Alleles must be specified by their
#'internal coding, i.e. '0' for ancestral resp. major allele, etc.
#'@param col color for each allele (as coded internally).
#'@param mrk.col color of the vertical line at the focal marker position.
#'@param lwd controls the relative width of the diagram lines on the plot (default 0.1).
#'@param hap.names a vector containing names of chromosomes.
#'@param cex.lab relative size of labels. See \code{\link[graphics]{par}}.
#'@param family.lab font family for labels. See \code{\link[graphics]{par}}.
#'@param offset.lab offset of labels. See \code{\link[graphics]{par}}.
#'@param legend legend text.
#'@param legend.xy.coords if \code{"automatic"} (default) places legend either top left or top right;
#'if \code{"none"}, no legend is drawn; otherwise argument is passed to \code{\link[graphics]{legend}}.
#'@param ... other arguments to be passed to \code{\link[graphics]{plot.default}}.
#'@seealso \code{\link{plot.haplen}}.
#'@examples #example haplohh object (280 haplotypes, 1424 SNPs)
#'#see ?haplohh_cgu_bta12 for details
#'data(haplohh_cgu_bta12)
#'#plotting furcation diagram for both ancestral and derived allele
#'#from the marker "F1205400"
#'#which display a strong signal of selection
#'f <- calc_furcation(haplohh_cgu_bta12, mrk = "F1205400")
#'plot(f)
#'plot(f, xlim = c(2e+07,3.5e+07))
#'plot(f, xlim = c(2.7e+07,3.1e+07))
#'plot(f, xlim = c(2.7e+07,3.1e+07), hap.names = hap.names(haplohh_cgu_bta12), cex.lab=0.3)
#'@export
#'@importFrom stats na.omit
plot.furcation <-
  function(x,
           allele = NA,
           col = c("blue", "red", "violet", "orange"),
           mrk.col = "gray",
           lwd = 0.1,
           hap.names = NULL,
           cex.lab = 1.0,
           family.lab = "",
           offset.lab = 0.5,
           legend = NA,
           legend.xy.coords = "automatic",
           ...) {
    ##check parameters
    if (!is.furcation(x)) {
      stop("The data is not a valid object of a furcation.", call. = FALSE)
    }
    
    if (!is.null(hap.names)) {
      if (length(hap.names) != x@nhap) {
        stop(
          "Number of specified haplotype names has to be equal to number of haplotypes.",
          call. = FALSE
        )
      }
    }
    
    if (!anyNA(allele) & !is.numeric(allele)) {
      stop("Alleles have to be specified by integers 0,1,2,...", call. = FALSE)
    }
    
    allele_in_furcation <- as.integer(names(x))
    if (!anyNA(allele)) {
      if (anyDuplicated(allele)) {
        stop("Repeated specification of the same allele.", call. = FALSE)
      }
      for (i in allele) {
        if (!(i %in% allele_in_furcation)) {
          stop(paste0("No allele ", i, " found in furcation."), call. = FALSE)
        }
      }
    } else{
      allele <- allele_in_furcation
    }
    
    dot_args <- list(...)
    
    if (!is.null(dot_args$xlim)) {
      if (!is.numeric(dot_args$xlim) |
          length(dot_args$xlim) != 2 |
          dot_args$xlim[2] < dot_args$xlim[1]) {
        stop("Incorrect specification of xlim.", call. = FALSE)
      }
      if (!((dot_args$xlim)[1] <= x@position &
            x@position <= (dot_args$xlim)[2])) {
        stop("Focal marker must lie in specified xlim.", call. = FALSE)
      }
    }
    
    if (!anyNA(legend.xy.coords)) {
      if (is.numeric(legend.xy.coords) & length(legend.xy.coords) == 2) {
        legend.x <- legend.xy.coords[1]
        legend.y <- legend.xy.coords[2]
      } else if (length(legend.xy.coords) == 1) {
        legend.x <- legend.xy.coords
        legend.y <- NULL
      } else{
        stop("Invalid legend coordinate specification.", call. = FALSE)
      }
    } else{
      legend.x <- "none"
    }
    
    # perform plot
    
    # NOTE that the allele in object furcation has to be addressed by its
    # name and not its index! The name is always a character!
    
    #number of chromosomes without missing values at focal marker
    
    nhap_without_missing <- 0
    for (i in as.character(allele)) {
      nhap_without_missing <- nhap_without_missing + x[[i]]@count
    }
    
    xmin <- x@position
    for (i in as.character(allele)) {
      #if there are unresolved nodes (nodes with more than one chromosome),
      #than set region to be showen to first resp. last marker
      if (anyDuplicated(na.omit(x[[i]]@left@label_parent))) {
        xmin <- x@xlim[1]
        break
      }
      else{
        if (min(x[[i]]@left@node_pos) < xmin) {
          xmin <- min(x[[i]]@left@node_pos)
        }
      }
    }
    
    xmax <- x@position
    for (i in as.character(allele)) {
      if (anyDuplicated(na.omit(x[[i]]@right@label_parent))) {
        xmax <- x@xlim[2]
        break
      }
      else{
        if (max(x[[i]]@right@node_pos) > xmax) {
          xmax <- max(x[[i]]@right@node_pos)
        }
      }
    }
    
    if (is.null(dot_args$xlim)) {
      dot_args$xlim <- c(xmin, xmax)
    } else{
      if (dot_args$xlim[1] > xmin) {
        for (i in as.character(allele)) {
          f <- x[[i]]
          f@left <- prune_tree(f@left, dot_args$xlim[1])
          x[[i]] <- f
        }
      }
      if (dot_args$xlim[2] < xmax) {
        for (i in as.character(allele)) {
          f <- x[[i]]
          f@right <- prune_tree(f@right, dot_args$xlim[2])
          x[[i]] <- f
        }
      }
    }
    
    y_l <- list()
    y_r <- list()
    
    #y coordinates of nodes for each furcation
    for (i in seq_along(allele)) {
      y_l[[i]] <- calc_y(x[[as.character(allele[i])]]@left@node_parent,
                         x[[as.character(allele[i])]]@left@label_parent)
      y_r[[i]] <-
        calc_y(x[[as.character(allele[i])]]@right@node_parent,
               x[[as.character(allele[i])]]@right@label_parent)
    }
    
    #transform y values of each plot to interval [0,1]
    for (i in seq_along(allele)) {
      y_min <- min(y_l[[i]], y_r[[i]])
      y_l[[i]] <- y_l[[i]] - y_min
      y_r[[i]] <- y_r[[i]] - y_min
      y_max <- max(y_l[[i]], y_r[[i]])
      if (y_max != 0) {
        y_l[[i]] <- y_l[[i]] / y_max
        y_r[[i]] <- y_r[[i]] / y_max
      }
    }
    
    if (length(allele) > 1) {
      #weight plots by their number of chromosomes
      for (i in seq_along(allele)) {
        y_l[[i]] <-
          y_l[[i]] * x[[as.character(allele[i])]]@count / nhap_without_missing
        y_r[[i]] <-
          y_r[[i]] * x[[as.character(allele[i])]]@count / nhap_without_missing
      }
      
      #add some space between sub-plots
      y_max <- max(y_l[[length(allele)]], y_r[[length(allele)]])
      offset <- 1 / nhap_without_missing
      for (i in (length(allele) - 1):1) {
        y_l[[i]] <- y_l[[i]] + y_max + offset
        y_r[[i]] <- y_r[[i]] + y_max + offset
        y_max <- max(y_l[[i]], y_r[[i]])
      }
      
      #rescale whole plot to interval [0,1]
      for (i in seq_along(allele)) {
        y_l[[i]] <- y_l[[i]] / y_max
        y_r[[i]] <- y_r[[i]] / y_max
      }
    } else{
      if (max(abs(y_l[[i]]), abs(y_r[[i]])) == 0) {
        # if only one plot and no furcation
        y_l[[i]] <- y_l[[i]] + 0.5   # center vertically
        y_r[[i]] <- y_r[[i]] + 0.5   # center vertically
      }
    }
    
    p <- floor(log(dot_args$xlim[2], 1000))
    ## only shrink big scales, but never magnify small ones (p<0)
    scale <- 1000 ** max(0, p)
    ## no unit if p < 0
    unit <- c("", "(bp)", "(kb)", "(Mb)", "(Gb)")[max(-1, p)  + 2]
    
    dot_args$xlim <- dot_args$xlim / scale
    
    dot_args$ylim <- c(0, 1)
    
    if (is.null(dot_args$xlab)) {
      dot_args$xlab <- paste("Position", unit)
    }
    
    if (is.null(dot_args$ylab)) {
      dot_args$ylab <- ""
    }
    
    if (is.null(dot_args$bty)) {
      dot_args$bty <- "n"
    }
    
    if (is.null(dot_args$yaxt)) {
      dot_args$yaxt <- "n"
    }
    
    if (is.null(dot_args$main)) {
      dot_args$main <-
        paste0("Haplotype furcations around '", x@mrk.name, "'")
    }
    
    # get descriptions of alleles
    descriptions <- vapply(allele, function(a) {
      x[[as.character(a)]]@description
    }, USE.NAMES = FALSE, FUN.VALUE = "")
    description_colors <- get_description_colors(descriptions, col)
    
    #invisible plot to get coordinate system
    do.call("plot", c(list(NULL),
                      dot_args))
    
    #dashed vertical line at focal mrk
    abline(v = x@position / scale,
           lty = 2,
           col = mrk.col)
    
    for (i in seq_along(allele)) {
      #left and right tree
      allelefurcation <- x[[as.character(allele[i])]]
      if (allelefurcation@count > 0) {
        draw_tree(
          allelefurcation@left,
          y_l[[i]],
          scale,
          max(dot_args$xlim[1], x@xlim[1] / scale),
          hap.names,
          lwd,
          description_colors[i],
          cex.lab,
          family.lab,
          offset.lab
        )
        draw_tree(
          allelefurcation@right,
          y_r[[i]],
          scale,
          min(dot_args$xlim[2], x@xlim[2] / scale),
          hap.names,
          lwd,
          description_colors[i],
          cex.lab,
          family.lab,
          offset.lab
        )
      }
    }
    if (legend.x != "none") {
      if (legend.x == "automatic") {
        picturemiddle <- mean(dot_args$xlim)
        legend.x <-
          ifelse(x@position / scale > picturemiddle, "topleft", "topright")
      } else if (is.numeric(legend.x)) {
        legend.x <- legend.x / scale
      }
      
      if (is.na(legend)) {
        legend <- descriptions
      }
      
      legend(
        legend.x,
        legend.y,
        legend = legend,
        bty = dot_args$bty,
        col = description_colors,
        lwd = 1,
        xpd = TRUE
      )
    }
    
  }

#draws a furcation tree
draw_tree <- function(ftree,
                      y,
                      scale,
                      xcutoff,
                      hap.names = NULL,
                      lwd,
                      col,
                      cex.lab,
                      family.lab,
                      offset.lab) {
  node_size <- calc_node_size(ftree@node_parent, ftree@label_parent)
  #are we left or right of the root node?
  direction <- sign(ftree@node_pos[1] / scale - xcutoff)
  node_pos <- ftree@node_pos / scale
  
  #lines from border til leaves of tree
  for (node in ftree@label_parent) {
    if (!is.na(node)) {
      lines(
        c(xcutoff, node_pos[node]),
        c(y[node], y[node]),
        col = col,
        lwd = lwd * node_size[node],
        lty = 1 + (ftree@node_with_missing_data[node]) * 1
      )
    }
  }
  
  if (length(ftree@node_parent) > 1) {
    already_visited_parent <- vector("integer")
    for (node in (seq_along(ftree@node_parent)[-1])) {
      parent <- ftree@node_parent[node]
      # vertical line at furcation
      lines(
        c(node_pos[node], node_pos[node]),
        c(y[node], y[parent]),
        lwd = lwd * node_size[node],
        col = col,
        lty = 1 + (ftree@node_with_missing_data[node]) * 1
      )
      if (!(parent %in% already_visited_parent)) {
        # horizontal line from furcation til parent
        lines(
          c(node_pos[node], node_pos[parent]),
          c(y[parent], y[parent]),
          lwd = lwd * node_size[parent],
          col = col
        )
        # draw horizontal line only once
        already_visited_parent <-
          append(already_visited_parent, parent)
      }
    }
  }
  
  if (!is.null(hap.names) & direction != 0) {
    labelled <- vector("numeric")
    for (i in which(!is.na(ftree@label_parent))) {
      hap.is.singleton <- node_size[ftree@label_parent[i]] == 1
      if (!hap.is.singleton) {
        if (ftree@label_parent[i] %in% labelled) {
          next
        } else{
          labelled <- append(labelled, ftree@label_parent[i])
        }
      }
      pos <- 3 - direction
      text(
        x = xcutoff,
        y = y[ftree@label_parent[i]],
        labels = hap.names[i],
        pos = pos,
        cex = cex.lab,
        family = family.lab,
        offset = offset.lab,
        font = 2 - hap.is.singleton * 1,
        xpd = TRUE
      )
    }
  }
  
}

prune_tree <- function(ftree, xcutoff) {
  if (length(ftree@node_parent) > 1) {
    direction <- sign(ftree@node_pos[1] - mean(ftree@node_pos[-1]))
    if (direction != 0) {
      for (node in length(ftree@node_parent):2) {
        if (ftree@node_pos[node] * direction < xcutoff * direction) {
          for (chr in which(!is.na(ftree@label_parent))) {
            if (ftree@label_parent[chr] == node) {
              ftree@label_parent[chr] <- ftree@node_parent[node]
            }
          }
          ftree@node_parent <- ftree@node_parent[-node]
          ftree@node_pos <- ftree@node_pos[-node]
          ftree@node_with_missing_data <-
            ftree@node_with_missing_data[-node]
        }
      }
    }
  }
  return(ftree)
}

#calculates how many haplotypes are present at a node
calc_node_size <- function(node_parent, label_parent) {
  node_size <- rep(0, length(node_parent))
  #start with leaves
  for (node in seq_along(label_parent)) {
    node_size[label_parent[node]] <- node_size[label_parent[node]] + 1
  }
  #continue with inner nodes
  for (node1 in (length(node_parent) - 1):0) {
    #assumes that index of parent node is less than index of node itself
    for (node2 in (node1 + 1):length(node_parent)) {
      if (!is.na(node_parent[node2]) & node_parent[node2] == node1) {
        node_size[node1] <- node_size[node1] + node_size[node2]
      }
    }
  }
  return(node_size)
}

#calculate vertical position of (inner) nodes
calc_y <- function(node_parent, label_parent) {
  y <- rep(0, length(node_parent))
  if (length(node_parent) > 1) {
    nodes_with_chr_label <- unique(label_parent[!is.na(label_parent)])
    children <- vector("list", length(node_parent))
    for (node in seq_along(node_parent)) {
      children[[node]] <- which(node_parent == node)
    }
    
    node <- 1
    child_counter <- 0
    children_index <- rep(1, length(node_parent))
    while (children_index[1] <= length(children[[1]])) {
      while (children_index[node] <= length(children[[node]])) {
        node <- children[[node]][children_index[node]]
      }
      if (node %in% nodes_with_chr_label) {
        #equal space between siblings; first node gets highest y value
        y[node] <-
          (length(nodes_with_chr_label) - child_counter) / length(nodes_with_chr_label)
        child_counter <- child_counter + 1
      } else{
        #parent gets middle height of children
        for (child in children[[node]]) {
          y[node] <- y[node] + y[child] / length(children[[node]])
        }
      }
      node <- node_parent[node]
      children_index[node] <- children_index[node] + 1
    }
    for (child in children[[1]]) {
      y[1] <- y[1] + y[child] / length(children[[1]])
    }
    #transform to set root at zero (necessary for aligning left and right tree)
    y <- y - y[1]
  }
  return(y)
}
