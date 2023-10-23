#!Rscript
""" Accesory functions.

    Author: Valentin Maurer <valentin.maurer@dkfz-heidelberg>
"""

library(RColorBrewer)

translate_pids = function(x){
  trans = factor(
    gsub(x, pattern = ".*([D|I]\\d{3}).*", replacement = "\\1"),
    levels = c("D117", "D129", "D217", "I217", "D213", "D124", "D123", "D360"),
    labels = c("P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8")
    )
  return(as.vector(trans))
}
translate_type = function(x){
  trans = factor(
    gsub(x, pattern = ".*_[cC]{1}(\\d{1,2}).*", replacement = "\\1"),
    levels = c("0", "1", "25"),
    labels = c("Negative", "Case", "Positive")
  )
  return(as.vector(trans))
}

translate_ss2names = function(x){
  x = gsub(basename(x), pattern = "(AS-\\d+-LR-\\d+).*", replacement = "\\1")
  x = data.table(NAME = x)
  metadata = fread("~/icgc/Valentin/JMMLT/processing/metadata/smartseq.csv")
  x = merge(x, metadata, by = "NAME", all.x = T)
  x[, PID   := translate_pids(PID)]
  x[, TYPE  := translate_type(TYPE)]
  x[, PLATE := as.numeric(PLATE)]
  x[, WELL  := toupper(WELL)]
  return(list(PID = x$PID, TYPE = x$TYPE, PLATE = x$PLATE, WELL = x$WELL))
}
translate_bsnames = function(x){
  pid = translate_pids(x)
  type = translate_type(x)
  plate = as.numeric(gsub(x, pattern = ".*_p(\\d{1,2}).*", replacement = "\\1"))
  well = toupper(gsub(x, pattern = ".*p\\d{1,2}-([a-z]\\d{1,2}).*", replacement = "\\1"))
  
  return(list(PID = pid, TYPE = type, PLATE = plate, WELL = well))
}
translate_names = function(x, modality){
  x = as.vector(x)
  if (modality == "bisulfiteseq"){
    return(translate_bsnames(x))
  }
  x = gsub(basename(x), pattern = "(AS-\\d+-LR-\\d+).*", replacement = "\\1")
  x = data.table(NAME = x)
  metadata = fread(
    file.path("~/icgc/Valentin/JMMLT/processing/metadata/", paste0(modality, ".csv"))
  )
  x[, N := 1:nrow(x)]
  x = merge(x, metadata, by = "NAME", all.x = T)
  x = x[order(N, decreasing = F)]
  x$N = NULL
  x[, PID   := translate_pids(PID)]
  x[, TYPE  := translate_type(TYPE)]
  x[, PLATE := as.numeric(PLATE)]
  x[, WELL  := toupper(WELL)]
  return(list(PID = x$PID, TYPE = x$TYPE, PLATE = x$PLATE, WELL = x$WELL))
}


cluster_na = function(m){
  m = t(apply(X = m, MARGIN = 1, FUN = function(d){
    d_median = median(d, na.rm = T)
    if(is.na(d_median)){d_median = 0}
    d[is.na(d)] = d_median
    d
  }))
  dist(m)
}

factor_chrom = function(vec){
  factor(
    vec,
    levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14",
               "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y", "MT"),
    labels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14",
               "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y", "MT")
  )
}


JMMLT_COLORSCHEME_FULL = list(
  Patient = c(  
    P1 = "#0058b4", P2 = "#2188c9", 
    P3 = "#fbbb25", P4 = "#fca349", 
    P5 = "#ff6b36", P6 = "#e34e2e",
    P7 = "#c33126", P8 = "#a41220"
  ),
  Control = c(
    Negative = "#0000FF", Case = "#FFFFFF", Positive = "#FF0000"
  ),
  Plate = c(
    "1" = "#FFFFFF", "2" = "#666666"
  ),
  Status = c(
    "Unknown" = "#BABABA",
    "WT" = "#2188c9",
    "Mut" = "#c33126"
  )
)



JMMLT_COLORSCHEME = list(
  Patient = c(  
    P1 = "#0058b4", P2 = "#2188c9", 
    P3 = "#fbbb25", P4 = "#fca349", 
    P5 = "#ff6b36", P6 = "#BABABA",
    P7 = "#c33126", P8 = "#a41220"
  ),
  Control = c(
    Negative = "#0000FF", Case = "#FFFFFF", Positive = "#FF0000"
  ),
  Plate = c(
    "1" = "#FFFFFF", "2" = "#666666"
  ),
  Status = c(
    "Unknown" = "#BABABA",
    "WT" = "#2188c9",
    "Mut" = "#c33126"
  )
)

SEQUENCED_PATIENTS = c("D123", "D124", "D129", "D217")

mymean = function(x, na.rm = T){
  return(mean(as.vector(unlist(x)), na.rm = na.rm))
}

split_vector = function(vec, chunks = NA, percent = NA, size = NA){
    if (sum(is.na(c(chunks, percent, size))) != 2) 
      stop("Invalid input. Must contain 1 of either chunks, percent, or size")
    if (!is.na(percent)) 
      chunks = 100/percent
    if (!is.na(size)) 
      chunks = length(vec)/ceiling(size)
    chunks = max(1, chunks)
    return(unname(split(vec, sort(rep_len(1:ceiling(chunks), 
                                          length(vec))))))
}

scaleNULLONE = function(x){
  (x - min(x))/(max(x) - min(x))
}

replaceRowNA = function(mat, FUN = median){
  t(
    apply(X = mat, MARGIN = 1, function(d){
      stat = FUN(na.omit(d))
      if(is.na(stat)){stat = 0}
      d[is.na(d)] = stat
      d
    }))
}

cor_na = function(mat, method = "pearson"){
  mat = t(apply(X = mat, MARGIN = 1, FUN = function(d){
    d_median = median(d, na.rm = T)
    if(is.na(d_median)){d_median = 0}
    d[is.na(d)] = d_median
    d
  }))
  cor(mat, method = method)
}

is_mut = function(x){
  # Return Mut if any mut; WT if WT and ./. in vector and ./. for ./.
  if (length(x) > 1){
    temp = unique(x)
    if("Mut" %in% temp){
      return("Mut")
    }
    if("Other mutation" %in% temp){
      return("Other mutation")
    }
    if("WT" %in% temp){
      return("WT")
    }
    if("./." %in% temp){
      return("./.")
    }
  }else if(length(x) == 0){
    return("./.")
  }
  return(x)
}

plot_tree = function(tree, col_df = NULL, col_id = 1, pch_id = NULL, alpha = 1, point_size = 1, type = "cladogram", ...){
  #tree = ape::read.tree(file = tree)
  if(class(tree) == "multiPhylo"){
    tree = tree[[1]]
  }
  
  if(!is.null(col_df)){
    cluster_cols = col_df[tree$tip.label, col_id]
    cluster_cols = grDevices::adjustcolor(col = cluster_cols, alpha.f = alpha)
  }
  
  xypos = plot.phylo2(x = tree, type = type, tip.color = cluster_cols, ...)
  
  if(is.null(pch_id)){
    ape::tiplabels(pch = 21, adj = c(0.6, 0.5), bg = cluster_cols, col = "black", cex = point_size)
  }else{
    ape::tiplabels(pch = as.numeric(as.character(col_df[tree$tip.label, pch_id])), adj = c(0.6, 0.5), bg = cluster_cols, col = 'black', cex = point_size)
  }
  
  invisible(xypos)
  
  #ape::tiplabels(pch = 21, adj = c(0.6, 0.5), bg = cluster_cols, col = "black", cex = point_size)
}



#This is the source code for ape::plot.phylo but modified to return invisible x and y coordinates of tips
plot.phylo2 = function (x, type = "phylogram", use.edge.length = TRUE, node.pos = NULL,
                        show.tip.label = TRUE, show.node.label = FALSE, edge.color = "black",
                        edge.width = 1, edge.lty = 1, font = 3, cex = par("cex"),
                        adj = NULL, srt = 0, no.margin = FALSE, root.edge = FALSE,
                        label.offset = 0, underscore = FALSE, x.lim = NULL, y.lim = NULL,
                        direction = "rightwards", lab4ut = NULL, tip.color = "black",
                        plot = TRUE, rotate.tree = 0, open.angle = 0, node.depth = 1,
                        align.tip.label = FALSE, ...)
{
  Ntip <- length(x$tip.label)
  if (Ntip < 2) {
    warning("found less than 2 tips in the tree")
    return(NULL)
  }
  .nodeHeight <- function(edge, Nedge, yy) .C(node_height,
                                              as.integer(edge[, 1]), as.integer(edge[, 2]), as.integer(Nedge),
                                              as.double(yy))[[4]]
  .nodeDepth <- function(Ntip, Nnode, edge, Nedge, node.depth) .C(node_depth,
                                                                  as.integer(Ntip), as.integer(edge[, 1]), as.integer(edge[,
                                                                                                                           2]), as.integer(Nedge), double(Ntip + Nnode), as.integer(node.depth))[[5]]
  .nodeDepthEdgelength <- function(Ntip, Nnode, edge, Nedge,
                                   edge.length) .C(node_depth_edgelength, as.integer(edge[,
                                                                                          1]), as.integer(edge[, 2]), as.integer(Nedge), as.double(edge.length),
                                                   double(Ntip + Nnode))[[5]]
  Nedge <- dim(x$edge)[1]
  Nnode <- x$Nnode
  if (any(x$edge < 1) || any(x$edge > Ntip + Nnode))
    stop("tree badly conformed; cannot plot. Check the edge matrix.")
  ROOT <- Ntip + 1
  type <- match.arg(type, c("phylogram", "cladogram", "fan",
                            "unrooted", "radial"))
  direction <- match.arg(direction, c("rightwards", "leftwards",
                                      "upwards", "downwards"))
  if (is.null(x$edge.length)) {
    use.edge.length <- FALSE
  }
  else {
    if (use.edge.length && type != "radial") {
      tmp <- sum(is.na(x$edge.length))
      if (tmp) {
        warning(paste(tmp, "branch length(s) NA(s): branch lengths ignored in the plot"))
        use.edge.length <- FALSE
      }
    }
  }
  if (is.numeric(align.tip.label)) {
    align.tip.label.lty <- align.tip.label
    align.tip.label <- TRUE
  }
  else {
    if (align.tip.label)
      align.tip.label.lty <- 3
  }
  if (align.tip.label) {
    if (type %in% c("unrooted", "radial") || !use.edge.length ||
        is.ultrametric(x))
      align.tip.label <- FALSE
  }
  if (type %in% c("unrooted", "radial") || !use.edge.length ||
      is.null(x$root.edge) || !x$root.edge)
    root.edge <- FALSE
  phyloORclado <- type %in% c("phylogram", "cladogram")
  horizontal <- direction %in% c("rightwards", "leftwards")
  xe <- x$edge
  if (phyloORclado) {
    phyOrder <- attr(x, "order")
    if (is.null(phyOrder) || phyOrder != "cladewise") {
      x <- reorder(x)
      if (!identical(x$edge, xe)) {
        ereorder <- match(x$edge[, 2], xe[, 2])
        if (length(edge.color) > 1) {
          edge.color <- rep(edge.color, length.out = Nedge)
          edge.color <- edge.color[ereorder]
        }
        if (length(edge.width) > 1) {
          edge.width <- rep(edge.width, length.out = Nedge)
          edge.width <- edge.width[ereorder]
        }
        if (length(edge.lty) > 1) {
          edge.lty <- rep(edge.lty, length.out = Nedge)
          edge.lty <- edge.lty[ereorder]
        }
      }
    }
    yy <- numeric(Ntip + Nnode)
    TIPS <- x$edge[x$edge[, 2] <= Ntip, 2]
    yy[TIPS] <- 1:Ntip
  }
  z <- reorder(x, order = "postorder")
  if (phyloORclado) {
    if (is.null(node.pos))
      node.pos <- if (type == "cladogram" && !use.edge.length)
        2
    else 1
    if (node.pos == 1)
      yy <- .nodeHeight(z$edge, Nedge, yy)
    else {
      ans <- .C(node_height_clado, as.integer(Ntip), as.integer(z$edge[,
                                                                       1]), as.integer(z$edge[, 2]), as.integer(Nedge),
                double(Ntip + Nnode), as.double(yy))
      xx <- ans[[5]] - 1
      yy <- ans[[6]]
    }
    if (!use.edge.length) {
      if (node.pos != 2)
        xx <- .nodeDepth(Ntip, Nnode, z$edge, Nedge,
                         node.depth) - 1
      xx <- max(xx) - xx
    }
    else {
      xx <- .nodeDepthEdgelength(Ntip, Nnode, z$edge, Nedge,
                                 z$edge.length)
    }
  }
  else {
    twopi <- 2 * pi
    rotate.tree <- twopi * rotate.tree/360
    if (type != "unrooted") {
      TIPS <- x$edge[which(x$edge[, 2] <= Ntip), 2]
      xx <- seq(0, twopi * (1 - 1/Ntip) - twopi * open.angle/360,
                length.out = Ntip)
      theta <- double(Ntip)
      theta[TIPS] <- xx
      theta <- c(theta, numeric(Nnode))
    }
    switch(type, fan = {
      theta <- .nodeHeight(z$edge, Nedge, theta)
      if (use.edge.length) {
        r <- .nodeDepthEdgelength(Ntip, Nnode, z$edge,
                                  Nedge, z$edge.length)
      } else {
        r <- .nodeDepth(Ntip, Nnode, z$edge, Nedge, node.depth)
        r <- 1/r
      }
      theta <- theta + rotate.tree
      if (root.edge) r <- r + x$root.edge
      xx <- r * cos(theta)
      yy <- r * sin(theta)
    }, unrooted = {
      nb.sp <- .nodeDepth(Ntip, Nnode, z$edge, Nedge, node.depth)
      XY <- if (use.edge.length) unrooted.xy(Ntip, Nnode,
                                             z$edge, z$edge.length, nb.sp, rotate.tree) else unrooted.xy(Ntip,
                                                                                                         Nnode, z$edge, rep(1, Nedge), nb.sp, rotate.tree)
      xx <- XY$M[, 1] - min(XY$M[, 1])
      yy <- XY$M[, 2] - min(XY$M[, 2])
    }, radial = {
      r <- .nodeDepth(Ntip, Nnode, z$edge, Nedge, node.depth)
      r[r == 1] <- 0
      r <- 1 - r/Ntip
      theta <- .nodeHeight(z$edge, Nedge, theta) + rotate.tree
      xx <- r * cos(theta)
      yy <- r * sin(theta)
    })
  }
  if (phyloORclado) {
    if (!horizontal) {
      tmp <- yy
      yy <- xx
      xx <- tmp - min(tmp) + 1
    }
    if (root.edge) {
      if (direction == "rightwards")
        xx <- xx + x$root.edge
      if (direction == "upwards")
        yy <- yy + x$root.edge
    }
  }
  if (no.margin)
    par(mai = rep(0, 4))
  if (show.tip.label)
    nchar.tip.label <- nchar(x$tip.label)
  max.yy <- max(yy)
  getLimit <- function(x, lab, sin, cex) {
    s <- strwidth(lab, "inches", cex = cex)
    if (any(s > sin))
      return(1.5 * max(x))
    Limit <- 0
    while (any(x > Limit)) {
      i <- which.max(x)
      alp <- x[i]/(sin - s[i])
      Limit <- x[i] + alp * s[i]
      x <- x + alp * s
    }
    Limit
  }
  if (is.null(x.lim)) {
    if (phyloORclado) {
      if (horizontal) {
        xx.tips <- xx[1:Ntip]
        if (show.tip.label) {
          pin1 <- par("pin")[1]
          tmp <- getLimit(xx.tips, x$tip.label, pin1,
                          cex)
          tmp <- tmp + label.offset
        }
        else tmp <- max(xx.tips)
        x.lim <- c(0, tmp)
      }
      else x.lim <- c(1, Ntip)
    }
    else switch(type, fan = {
      if (show.tip.label) {
        offset <- max(nchar.tip.label * 0.018 * max.yy *
                        cex)
        x.lim <- range(xx) + c(-offset, offset)
      } else x.lim <- range(xx)
    }, unrooted = {
      if (show.tip.label) {
        offset <- max(nchar.tip.label * 0.018 * max.yy *
                        cex)
        x.lim <- c(0 - offset, max(xx) + offset)
      } else x.lim <- c(0, max(xx))
    }, radial = {
      if (show.tip.label) {
        offset <- max(nchar.tip.label * 0.03 * cex)
        x.lim <- c(-1 - offset, 1 + offset)
      } else x.lim <- c(-1, 1)
    })
  }
  else if (length(x.lim) == 1) {
    x.lim <- c(0, x.lim)
    if (phyloORclado && !horizontal)
      x.lim[1] <- 1
    if (type %in% c("fan", "unrooted") && show.tip.label)
      x.lim[1] <- -max(nchar.tip.label * 0.018 * max.yy *
                         cex)
    if (type == "radial")
      x.lim[1] <- if (show.tip.label)
        -1 - max(nchar.tip.label * 0.03 * cex)
    else -1
  }
  if (phyloORclado && direction == "leftwards")
    xx <- x.lim[2] - xx
  if (is.null(y.lim)) {
    if (phyloORclado) {
      if (horizontal)
        y.lim <- c(1, Ntip)
      else {
        pin2 <- par("pin")[2]
        yy.tips <- yy[1:Ntip]
        if (show.tip.label) {
          tmp <- getLimit(yy.tips, x$tip.label, pin2,
                          cex)
          tmp <- tmp + label.offset
        }
        else tmp <- max(yy.tips)
        y.lim <- c(0, tmp)
      }
    }
    else switch(type, fan = {
      if (show.tip.label) {
        offset <- max(nchar.tip.label * 0.018 * max.yy *
                        cex)
        y.lim <- c(min(yy) - offset, max.yy + offset)
      } else y.lim <- c(min(yy), max.yy)
    }, unrooted = {
      if (show.tip.label) {
        offset <- max(nchar.tip.label * 0.018 * max.yy *
                        cex)
        y.lim <- c(0 - offset, max.yy + offset)
      } else y.lim <- c(0, max.yy)
    }, radial = {
      if (show.tip.label) {
        offset <- max(nchar.tip.label * 0.03 * cex)
        y.lim <- c(-1 - offset, 1 + offset)
      } else y.lim <- c(-1, 1)
    })
  }
  else if (length(y.lim) == 1) {
    y.lim <- c(0, y.lim)
    if (phyloORclado && horizontal)
      y.lim[1] <- 1
    if (type %in% c("fan", "unrooted") && show.tip.label)
      y.lim[1] <- -max(nchar.tip.label * 0.018 * max.yy *
                         cex)
    if (type == "radial")
      y.lim[1] <- if (show.tip.label)
        -1 - max(nchar.tip.label * 0.018 * max.yy * cex)
    else -1
  }
  if (phyloORclado && direction == "downwards")
    yy <- y.lim[2] - yy
  if (phyloORclado && root.edge) {
    if (direction == "leftwards")
      x.lim[2] <- x.lim[2] + x$root.edge
    if (direction == "downwards")
      y.lim[2] <- y.lim[2] + x$root.edge
  }
  asp <- if (type %in% c("fan", "radial", "unrooted"))
    1
  else NA
  plot.default(0, type = "n", xlim = x.lim, ylim = y.lim, xlab = "",
               ylab = "", axes = FALSE, asp = asp, ...)
  if (plot) {
    if (is.null(adj))
      adj <- if (phyloORclado && direction == "leftwards")
        1
    else 0
    if (phyloORclado && show.tip.label) {
      MAXSTRING <- max(strwidth(x$tip.label, cex = cex))
      loy <- 0
      if (direction == "rightwards") {
        lox <- label.offset + MAXSTRING * 1.05 * adj
      }
      if (direction == "leftwards") {
        lox <- -label.offset - MAXSTRING * 1.05 * (1 -
                                                     adj)
      }
      if (!horizontal) {
        psr <- par("usr")
        MAXSTRING <- MAXSTRING * 1.09 * (psr[4] - psr[3])/(psr[2] -
                                                             psr[1])
        loy <- label.offset + MAXSTRING * 1.05 * adj
        lox <- 0
        srt <- 90 + srt
        if (direction == "downwards") {
          loy <- -loy
          srt <- 180 + srt
        }
      }
    }
    if (type == "phylogram") {
      phylogram.plot(x$edge, Ntip, Nnode, xx, yy, horizontal,
                     edge.color, edge.width, edge.lty)
    }
    else {
      if (type == "fan") {
        ereorder <- match(z$edge[, 2], x$edge[, 2])
        if (length(edge.color) > 1) {
          edge.color <- rep(edge.color, length.out = Nedge)
          edge.color <- edge.color[ereorder]
        }
        if (length(edge.width) > 1) {
          edge.width <- rep(edge.width, length.out = Nedge)
          edge.width <- edge.width[ereorder]
        }
        if (length(edge.lty) > 1) {
          edge.lty <- rep(edge.lty, length.out = Nedge)
          edge.lty <- edge.lty[ereorder]
        }
        circular.plot(z$edge, Ntip, Nnode, xx, yy, theta,
                      r, edge.color, edge.width, edge.lty)
      }
      else cladogram.plot(x$edge, xx, yy, edge.color, edge.width,
                          edge.lty)
    }
    if (root.edge) {
      rootcol <- if (length(edge.color) == 1)
        edge.color
      else "black"
      rootw <- if (length(edge.width) == 1)
        edge.width
      else 1
      rootlty <- if (length(edge.lty) == 1)
        edge.lty
      else 1
      if (type == "fan") {
        tmp <- polar2rect(x$root.edge, theta[ROOT])
        segments(0, 0, tmp$x, tmp$y, col = rootcol, lwd = rootw,
                 lty = rootlty)
      }
      else {
        switch(direction, rightwards = segments(0, yy[ROOT],
                                                x$root.edge, yy[ROOT], col = rootcol, lwd = rootw,
                                                lty = rootlty), leftwards = segments(xx[ROOT],
                                                                                     yy[ROOT], xx[ROOT] + x$root.edge, yy[ROOT],
                                                                                     col = rootcol, lwd = rootw, lty = rootlty),
               upwards = segments(xx[ROOT], 0, xx[ROOT], x$root.edge,
                                  col = rootcol, lwd = rootw, lty = rootlty),
               downwards = segments(xx[ROOT], yy[ROOT], xx[ROOT],
                                    yy[ROOT] + x$root.edge, col = rootcol, lwd = rootw,
                                    lty = rootlty))
      }
    }
    if (show.tip.label) {
      if (is.expression(x$tip.label))
        underscore <- TRUE
      if (!underscore)
        x$tip.label <- gsub("_", " ", x$tip.label)
      if (phyloORclado) {
        if (align.tip.label) {
          xx.tmp <- switch(direction, rightwards = max(xx[1:Ntip]),
                           leftwards = min(xx[1:Ntip]), upwards = xx[1:Ntip],
                           downwards = xx[1:Ntip])
          yy.tmp <- switch(direction, rightwards = yy[1:Ntip],
                           leftwards = yy[1:Ntip], upwards = max(yy[1:Ntip]),
                           downwards = min(yy[1:Ntip]))
          segments(xx[1:Ntip], yy[1:Ntip], xx.tmp, yy.tmp,
                   lty = align.tip.label.lty)
        }
        else {
          xx.tmp <- xx[1:Ntip]
          yy.tmp <- yy[1:Ntip]
        }
        text(xx.tmp + lox, yy.tmp + loy, x$tip.label,
             adj = adj, font = font, srt = srt, cex = cex,
             col = tip.color)
      }
      else {
        angle <- if (type == "unrooted")
          XY$axe
        else atan2(yy[1:Ntip], xx[1:Ntip])
        lab4ut <- if (is.null(lab4ut)) {
          if (type == "unrooted")
            "horizontal"
          else "axial"
        }
        else match.arg(lab4ut, c("horizontal", "axial"))
        xx.tips <- xx[1:Ntip]
        yy.tips <- yy[1:Ntip]
        if (label.offset) {
          xx.tips <- xx.tips + label.offset * cos(angle)
          yy.tips <- yy.tips + label.offset * sin(angle)
        }
        if (lab4ut == "horizontal") {
          y.adj <- x.adj <- numeric(Ntip)
          sel <- abs(angle) > 0.75 * pi
          x.adj[sel] <- -strwidth(x$tip.label)[sel] *
            1.05
          sel <- abs(angle) > pi/4 & abs(angle) < 0.75 *
            pi
          x.adj[sel] <- -strwidth(x$tip.label)[sel] *
            (2 * abs(angle)[sel]/pi - 0.5)
          sel <- angle > pi/4 & angle < 0.75 * pi
          y.adj[sel] <- strheight(x$tip.label)[sel]/2
          sel <- angle < -pi/4 & angle > -0.75 * pi
          y.adj[sel] <- -strheight(x$tip.label)[sel] *
            0.75
          text(xx.tips + x.adj * cex, yy.tips + y.adj *
                 cex, x$tip.label, adj = c(adj, 0), font = font,
               srt = srt, cex = cex, col = tip.color)
        }
        else {
          if (align.tip.label) {
            POL <- rect2polar(xx.tips, yy.tips)
            POL$r[] <- max(POL$r)
            REC <- polar2rect(POL$r, POL$angle)
            xx.tips <- REC$x
            yy.tips <- REC$y
            segments(xx[1:Ntip], yy[1:Ntip], xx.tips,
                     yy.tips, lty = align.tip.label.lty)
          }
          if (type == "unrooted") {
            adj <- abs(angle) > pi/2
            angle <- angle * 180/pi
            angle[adj] <- angle[adj] - 180
            adj <- as.numeric(adj)
          }
          else {
            s <- xx.tips < 0
            angle <- angle * 180/pi
            angle[s] <- angle[s] + 180
            adj <- as.numeric(s)
          }
          font <- rep(font, length.out = Ntip)
          tip.color <- rep(tip.color, length.out = Ntip)
          cex <- rep(cex, length.out = Ntip)
          for (i in 1:Ntip) text(xx.tips[i], yy.tips[i],
                                 x$tip.label[i], font = font[i], cex = cex[i],
                                 srt = angle[i], adj = adj[i], col = tip.color[i])
        }
      }
    }
    if (show.node.label)
      text(xx[ROOT:length(xx)] + label.offset, yy[ROOT:length(yy)],
           x$node.label, adj = adj, font = font, srt = srt,
           cex = cex)
  }
  L <- list(type = type, use.edge.length = use.edge.length,
            node.pos = node.pos, node.depth = node.depth, show.tip.label = show.tip.label,
            show.node.label = show.node.label, font = font, cex = cex,
            adj = adj, srt = srt, no.margin = no.margin, label.offset = label.offset,
            x.lim = x.lim, y.lim = y.lim, direction = direction,
            tip.color = tip.color, Ntip = Ntip, Nnode = Nnode, root.time = x$root.time,
            align.tip.label = align.tip.label)
  assign("last_plot.phylo", c(L, list(edge = xe, xx = xx, yy = yy)),
         envir = .PlotPhyloEnv)
  #invisible(L)
  #Return x and y position of tips
  list(xpos = xx.tips, ypos = yy.tips)
}


write_bed = function (m, output_dir = NULL, rm_NA = TRUE, force = FALSE, 
                      n_thr = 4, compress = TRUE, SeqStyle = "UCSC", multiBed = NULL, 
                      metilene = FALSE, phenoCol = NULL, add_coverage = FALSE) 
{
  if (!dir.exists(output_dir)) {
    dir.create(path = output_dir, showWarnings = FALSE, recursive = TRUE)
  }
  mat_gr <- methrix::get_matrix(m = m, type = "M", add_loci = TRUE, 
                                in_granges = TRUE)
  GenomeInfoDb::seqlevelsStyle(mat_gr) <- SeqStyle
  mat <- as.data.table(mat_gr)
  mat <- mat[, c("seqnames", "start", "end", "strand", rownames(colData(x = m))), 
             with = FALSE]
  if (add_coverage){
    cov_gr <- methrix::get_matrix(m = m, type = "C", add_loci = TRUE, 
                                  in_granges = TRUE)
    GenomeInfoDb::seqlevelsStyle(cov_gr) <- SeqStyle
    cov_mat <- as.data.table(cov_gr)
    cov_mat <- cov_mat[, c("seqnames", "start", "end", "strand", rownames(colData(x = m))),
                       with = FALSE] 
  }
  
  if (is.null(output_dir)) {
    output_dir <- getwd()
  }
  parameters <- c(color = "255,0,0", visibility = "full", altColor = "128,128,128", 
                  autoScale = "on", viewLimits = "0:1", windowingFunction = "mean")
  parameters <- paste0(" ", paste(names(parameters), parameters, 
                                  sep = "=", collapse = " "))
  if (!is.null(multiBed)) {
    if (compress) {
      op_bdg <- paste0(output_dir, "/", multiBed, ".bedGraph.gz")
    }
    else {
      op_bdg <- paste0(output_dir, "/", multiBed, ".bedGraph")
    }
    message("----------------------")
    message("*Writing ", op_bdg, " ")
    if (metilene) {
      if (is.null(phenoCol)) {
        stop("Please provide a value to phenoCol.")
      }
      if (!phenoCol %in% colnames(colData(m))) {
        stop(phenoCol, " is not a valid column name in colData().\n Available column name are: ", 
             paste(colnames(colData(m)), collapse = ", "))
      }
      colnames(mat)[1:2] = c("chr", "pos")
      colnames(mat)[5:ncol(mat)] = paste(as.character(colData(m)[, 
                                                                 phenoCol]), rownames(colData(m)), sep = "_")
      if (rm_NA) {
        mat = mat[complete.cases(mat), ]
      }
      data.table::fwrite(x = mat[, c(1, 2, 5:ncol(mat)), 
                                 with = FALSE], file = op_bdg, sep = "\t", na = ".", 
                         scipen = 7, nThread = n_thr, compress = "auto", 
                         showProgress = TRUE, quote = FALSE)
    }
    else {
      if (rm_NA) {
        mat = mat[complete.cases(mat), ]
      }
      colnames(mat)[1] = paste0("#", colnames(mat)[1])
      data.table::fwrite(x = mat, file = op_bdg, sep = "\t", 
                         na = ".", scipen = 7, nThread = n_thr, compress = "auto", 
                         showProgress = TRUE, quote = FALSE)
    }
  }
  else {
    message("----------------------")
    message("*Writing bedGraphs:")
    op_bdgs <- lapply(seq_len(nrow(colData(m))), function(i) {
      mat_i <- mat[, c(seq_len(3), i + 4), with = FALSE]
      if (add_coverage) {
        mat_i <- cbind(mat_i, cov_mat[, c(i + 4), with = FALSE])
      }
      if (rm_NA) {
        mat_i <- mat_i[complete.cases(mat_i), , ]
      }
      if (compress) {
        op_bdg <- paste0(output_dir, "/", rownames(colData(m))[i], 
                         ".bedGraph.gz")
      }
      else {
        op_bdg <- paste0(output_dir, "/", rownames(colData(m))[i], 
                         ".bedGraph")
      }
      
      header <- data.table(paste0("track type=bedGraph name=\"", 
                                  rownames(colData(m))[i], "\"", parameters))
      if (add_coverage) {
        header <- data.table(paste0("track type=bedGraph name=\"", 
                                    rownames(colData(m))[i], "\"", parameters, " coverage"))
      }
      
      if (file.exists(op_bdg)) {
        if (force) {
          message(paste0("**Writing ", rownames(colData(m))[i]))
          colnames(mat_i) <- paste0("V", seq_len(ncol(mat_i)))
          data.table::fwrite(x = header, file = op_bdg, 
                             sep = "\t", append = FALSE, quote = FALSE, 
                             col.names = FALSE, nThread = n_thr, scipen = 7, 
                             compress = "auto")
          data.table::fwrite(x = mat_i, file = op_bdg, 
                             sep = "\t", append = TRUE, col.names = FALSE, 
                             nThread = n_thr, scipen = 7, compress = "auto")
        }
        else {
          message(paste0("**File ", basename(op_bdg), 
                         " already exists. Skipped re-writing"))
        }
      }
      else {
        message(paste0("**Writing ", rownames(colData(m))[i]))
        colnames(mat_i) <- paste0("V", seq_len(ncol(mat_i)))
        data.table::fwrite(x = header, file = op_bdg, 
                           sep = "\t", col.names = FALSE, append = FALSE, 
                           quote = FALSE, nThread = n_thr)
        data.table::fwrite(x = mat_i, file = op_bdg, 
                           sep = "\t", col.names = FALSE, append = TRUE, 
                           nThread = n_thr)
      }
      op_bdg
    })
  }
  message("----------------------")
}

#' Writes bedGraphs from methrix object
#' @param m \code{\link{methrix}} object
#' @param output_dir Output directory name where the files should be saved.
#' If \code{NULL} creats a \code{tempdir}
#' @param n_thr Default 4.
#' @param rm_NA remove NAs
#' @param force forces to create files if they are existing
#' @param compress Whether to compress the output. Default TRUE
#' @param SeqStyle Default `UCSC` with `chr` prefix.
#' @param multiBed Default NULL. If provided a filename, a single bedGraph file with all samples included is generated.
#' @param metilene Default FALSE. If TRUE outputs bedgraphs in `metilene` format that can be directly used for DMR calling with `metilene`. This option works only when \code{multiBed = TRUE}.
#' @param phenoCol Default NULL. `condition` column from colData. Only applicable if \code{metilene = TRUE}
#' @param add_coverage Default FALSE. Whether to add a coverage column to the output. Only applicable if \code{multiBed = NULL}
#' @examples
#' data('methrix_data')
#' write_bedgraphs(m = methrix_data, output_dir = './temp')
#' #Export to metline format for DMR calling with metline
#' write_bedgraphs(m = methrix_data, output_dir = "./temp", rm_NA = FALSE, metilene = TRUE,multiBed = "metline_ip", phenoCol = "Condition")
#' @return writes bedgraph files to output
#' @export


write_bedgraphs <- function(m, output_dir = NULL, rm_NA = TRUE, force = FALSE, 
                            n_thr = 4, compress = TRUE, SeqStyle="UCSC", multiBed = NULL, metilene = FALSE, 
                            phenoCol = NULL, add_coverage = FALSE) {
  
  if (!dir.exists(output_dir)) {
    dir.create(path = output_dir, showWarnings = FALSE, recursive = TRUE)
  } else if (is.null(output_dir)) {
    output_dir <- getwd()
  }
  
  mat_gr <- methrix::get_matrix(m = m, type = "M", add_loci = TRUE, in_granges = TRUE)
  GenomeInfoDb::seqlevelsStyle(mat_gr)<- SeqStyle
  mat <- as.data.table(mat_gr) 
  mat <- mat[, c("seqnames", "start", "end", "strand", rownames(colData(x = m))), 
             with = FALSE]
  
  
  if (!is.null(multiBed)) {    
    op_bdg = file.path(output_dir, paste0(
      multiBed, ".bedGraph", ifelse(compress, yes = ".gz", no = ""))
    )
    if (rm_NA) {
      mat = mat[complete.cases(mat),]
    }
    
    if (metilene) {
      if (is.null(phenoCol)) {
        stop("Please provide a value to phenoCol.")
      }
      else if (!phenoCol %in% colnames(colData(m))) {
        stop(phenoCol, 
             " is not a valid column name in colData().\n Available column name are: ",
             paste(colnames(colData(m)), collapse = ", "))
      }
      
      colnames(mat)[1:2] = c("chr", "pos")
      colnames(mat)[5:ncol(mat)] = paste(as.character(colData(m)[,phenoCol]),
                                         rownames(colData(m)), sep = "_")
      mat = mat[,c(1,2,5:ncol(mat)), with = FALSE]
    } else {
      colnames(mat)[1] = paste0("#", colnames(mat)[1])
    }
    
    message("*Writing ", op_bdg, " ")
    data.table::fwrite(x = mat, file = op_bdg, sep = "\t", na = ".", scipen = 7,
                       nThread = n_thr, compress = "auto", showProgress = TRUE, quote = FALSE)
    
  } else {
    
    if (add_coverage){
      cov_gr <- methrix::get_matrix(m = m, type = "C", add_loci = TRUE, 
                                    in_granges = TRUE)
      GenomeInfoDb::seqlevelsStyle(cov_gr) <- SeqStyle
      cov_mat <- as.data.table(cov_gr)
      cov_mat <- cov_mat[, c("seqnames", "start", "end", "strand", 
                             rownames(colData(x = m))), with = FALSE] 
    }
    op_bdgs <- lapply(seq_len(nrow(colData(m))), function(i) {
      mat_i <- mat[, c(seq_len(3), i + 4), with = FALSE]
      if (add_coverage) {
        mat_i <- cbind(mat_i, cov_mat[, c(i + 4), with = FALSE])
      }
      if (rm_NA) {
        mat_i <- mat_i[complete.cases(mat_i), , ]
      }
      
      op_bdg = file.path(output_dir, paste0(
        rownames(colData(m))[i], ".bedGraph", ifelse(compress, yes = ".gz", no = ""))
      )
      
      # bedGraph header line
      parameters <-c("color" = "255,0,0", 
                     "visibility"="full", 
                     "altColor" = "128,128,128",
                     "autoScale"="on", 
                     "viewLimits"="0:1", 
                     "windowingFunction"="mean")  
      parameters <- paste0(" ", paste(names(parameters), parameters, 
                                      sep = "=", collapse = " "))
      header <- data.table(paste0('track type=bedGraph name="', 
                                  rownames(colData(m))[i], '"', parameters))
      
      if (file.exists(op_bdg) & !force) {
        message(paste0("*File ", basename(op_bdg), " already exists. Skipped re-writing"))
        return(op_bdg)
      }
      
      message(paste0("*Writing ", rownames(colData(m))[i]))
      colnames(mat_i) <- paste0("V", seq_len(ncol(mat_i)))
      data.table::fwrite(x = header, file = op_bdg, sep = "\t", append=FALSE, quote=FALSE,
                         col.names = FALSE, nThread = n_thr, scipen = 7, compress = "auto")
      data.table::fwrite(x = mat_i, file = op_bdg, sep = "\t", append=TRUE,
                         col.names = FALSE, nThread = n_thr, scipen = 7, compress = "auto")
      return(op_bdg)
    })
    
  }
  
  message("----------------------")
}

scale_na = function(x){ (x - mean(x, na.rm = T)/(sd(x, na.rm = T)))}

COLORSCHEME <- c(
  ADU = "#E5E5E5", JUV = "#A7D9CB", NEO = "#978474", FEL = "#6D7A9F", FES = "#595959",
  P1 = "#0058b4",  P2 = "#2188c9",  P3 = "#fbbb25",  P4 = "#fca349", 
  P5 = "#ff6b36",  P6 = "#e34e2e",  P7 = "#c33126",  P8 = "#a41220",
  LML = "#23395d", HML = "#6f0000", FET = "#BABABA"
)
scDMRMeth = function(liftover_path){
  dmr_anno = fread(liftover_path, col.names = c("chr", "chromStart", "chromEnd", "index"))
  dmr_anno$index = NULL
  dmr_anno[, chr := gsub(pattern = "chr", replacement = "", chr)]
  setkey(dmr_anno, chr, chromStart, chromEnd)
  bed_files = list.files("~/icgc/Valentin/JMMLT/processing/bisulfiteseq/filteredMethCallsPE/",
                         full.names = T)
  bed_files = bed_files[!grepl(pattern = "i7-only", bed_files)]
  
  rbindlist(mclapply(bed_files, function(i){
    data = fread(i)
    setkey(data, chr, start, end)
    data = foverlaps(data, dmr_anno,
                     type = "within",
                     by.x = c("chr", "start", "end"),
                     by.y = c("chr", "chromStart", "chromEnd"),
                     nomatch = 0L
    )
    data[, feature := paste0(chr, ":", chromStart, "-", chromEnd)]
    data = data[, .(n_total = sum(n_total, na.rm = T),
                    n_meth = sum(n_meth, na.rm = T),
                    avg_beta = mean(beta_value, na.rm = T)),
                by = .(feature)]
    data[, sample := paste(translate_bsnames(i)[c("PID", "PLATE", "WELL")], collapse = "_")]
  }, mc.cores = 10))
}
hg19ToHg38 = function(regions){
  seqlevelsStyle(regions) = "UCSC" 
  path = system.file(package="liftOver", "extdata", "hg19ToHg38.over.chain")
  ch = import.chain("~/icgc/Valentin/genomes/hg38_ensembl/hg19ToHg38.over.chain")
  unlist(liftOver(regions, ch))
}
tree_cols = function(data){
  pheno = readRDS(
    "~/icgc/Valentin/JMMLT/processing/bisulfiteseq/stats/bsseq_HSC_comb_snpRemoved_repMerged_pheno.rds"
  )
  cluster_cols = data.table(sample = colnames(data), donor = NA)
  cluster_cols[, donor := factor(sample, levels = rownames(pheno), labels = pheno$Donor)]
  cluster_cols[is.na(donor), donor := sample]
  cluster_cols[, Patient := translate_pids(donor)]
  cluster_cols[!is.na(Patient), donor := Patient][, Patient := NULL]
  cluster_cols[, donor := as.character(donor)]
  cluster_cols[!grepl(pattern = "P\\d", donor) & ! donor %in% c("LML", "HML"), donor := factor(
    donor, levels = c("cordblood", "adult_bonemarrow"), labels = c("NEO", "ADU"))]
  cluster_cols[grepl(pattern = "JU", sample), donor := "JUV"]
  cluster_cols[grepl(pattern = "FS", sample), donor := "FES"]
  cluster_cols[grepl(pattern = "FL", sample), donor := "FEL"]
  cluster_cols[, color := COLORSCHEME[donor]]
  cluster_cols[, donorCtype := donor]
  cluster_cols[donor == "ADU", donorCtype := paste(
    donor, gsub(sample, pattern = "(.*)_(.*)_(.*)\\d_NORMAL", replacement = "\\3"),
    sep = "_")]
  cluster_cols[sample == donor, name := donorCtype]
  cluster_cols[sample != donor, name := paste0(donorCtype, ".", 1:.N), by = donorCtype]
  cluster_cols[sample %in% c("P2", "P3"), donorCtype := paste(donorCtype, "scIGMT", sep = "_")]
  return(cluster_cols)
}
tree_cols2 = function(data_colnames){
  pheno = readRDS(
    "~/icgc/Valentin/JMMLT/processing/bisulfiteseq/stats/bsseq_HSC_comb_snpRemoved_repMerged_pheno.rds"
  )
  cluster_cols = data.table(sample = data_colnames, donor = NA)
  cluster_cols[, donor := factor(sample, levels = rownames(pheno), labels = pheno$Donor)]
  cluster_cols[is.na(donor), donor := sample]
  cluster_cols[, Patient := translate_pids(donor)]
  cluster_cols[!is.na(Patient), donor := Patient][, Patient := NULL]
  cluster_cols[, donor := as.character(donor)]
  cluster_cols[!grepl(pattern = "P\\d", donor) & ! donor %in% c("LML", "HML"), donor := factor(
    donor, levels = c("cordblood", "adult_bonemarrow"), labels = c("NEO", "ADU"))]
  cluster_cols[grepl(pattern = "JU", sample), donor := "JUV"]
  cluster_cols[grepl(pattern = "FS", sample), donor := "FES"]
  cluster_cols[grepl(pattern = "FL", sample), donor := "FEL"]
  cluster_cols[, color := COLORSCHEME[donor]]
  cluster_cols[, donorCtype := donor]
  cluster_cols[donor == "ADU", donorCtype := paste(
    donor, gsub(sample, pattern = "(.*)_(.*)_(.*)\\d_NORMAL", replacement = "\\3"),
    sep = "_")]
  cluster_cols[sample == donor, name := donorCtype]
  cluster_cols[sample != donor, name := paste0(donorCtype, ".", 1:.N), by = donorCtype]
  cluster_cols[sample %in% c("P2", "P3"), donorCtype := paste(donorCtype, "scIGMT", sep = "_")]
  return(cluster_cols)
}

EXPRESSION_COLOR = viridis::viridis(100)
ANNOTATION_COLORSCALE = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100)
ANNOTATION_COLORSCALE_W = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100)
ANNOTATION_COLORSCALE_FUN = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))
FACS_COLOR = colorRampPalette(c("#FFFFFF", "#AAAAAA","#000000"), bias = 1)(100)
FACS_COLOR = viridis::cividis(100)