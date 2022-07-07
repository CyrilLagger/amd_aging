####################################################
##
## Project: amd_aging
##
## D.Dhirachaikulpanich@liverpool.ac.uk
## cyril.lagger@liverpool.ac.uk
##
####################################################
##

## Libraries ####

library(Seurat)
library(scDiffCom)
library(CellChat)
library(data.table)
library(ggplot2)
library(future)
library(ComplexUpset)
library(escape)
library(ggVennDiagram)
library(ggExtra)
library(limma)
library(ggpubr)
library(cowplot)
library(ComplexHeatmap)
library(GEOquery)
library(GSVA)

## Options ####

plan(multicore, workers = 20)

## Directory paths ####

path_data <- "data/"
path_results <- "results/"

## Custom CellChat gnetAnalysis_signalingRole_heatmap functions #####

netAnalysis_signalingRole_heatmap_custom <- function(
  object,
  signaling = NULL,
  pattern = c("outgoing", "incoming", "all"),
  slot.name = "netP",
  color.use = NULL,
  color.heatmap = "BuGn",
  title = NULL,
  width = 10,
  height = 8,
  font.size = 8,
  font.size.title = 10,
  cluster.rows = FALSE,
  cluster.cols = FALSE
  ) {
  pattern <- match.arg(pattern)
  if (length(slot(object, slot.name)$centr) == 0) {
    stop("Please run `netAnalysis_computeCentrality`")
  }
  centr <- slot(object, slot.name)$centr
  outgoing <- matrix(0, nrow = nlevels(object@idents), ncol = length(centr))
  incoming <- matrix(0, nrow = nlevels(object@idents), ncol = length(centr))
  dimnames(outgoing) <- list(levels(object@idents), names(centr))
  dimnames(incoming) <- dimnames(outgoing)
  for (i in 1:length(centr)) {
    outgoing[, i] <- centr[[i]]$outdeg
    incoming[, i] <- centr[[i]]$indeg
  }
  if (pattern == "outgoing") {
    mat <- t(outgoing)
    legend.name <- "Outgoing"
  } else if (pattern == "incoming") {
    mat <- t(incoming)
    legend.name <- "Incoming"
  } else if (pattern == "all") {
    mat <- t(outgoing + incoming)
    legend.name <- "Overall"
  }
  if (is.null(title)) {
    title <- paste0(legend.name, " signaling patterns")
  } else {
    title <- paste0(paste0(legend.name, " signaling patterns"), " - ", title)
  }
  if (!is.null(signaling)) {
    mat1 <- mat[rownames(mat) %in% signaling, , drop = FALSE]
    mat <- matrix(0, nrow = length(signaling), ncol = ncol(mat))
    idx <- match(rownames(mat1), signaling)
    mat[idx[!is.na(idx)], ] <- mat1
    dimnames(mat) <- list(signaling, colnames(mat1))
  }
  mat.ori <- mat
  mat <- sweep(mat, 1L, apply(mat, 1, max), '/', check.margin = FALSE)
  mat[mat == 0] <- NA
  if (is.null(color.use)) {
    color.use <- scPalette(length(colnames(mat)))
  }
  color.heatmap.use <- grDevices::colorRampPalette((
    RColorBrewer::brewer.pal(n = 9, name = color.heatmap)
  ))(100)
  df <- data.frame(group = colnames(mat))
  rownames(df) <- colnames(mat)
  names(color.use) <- colnames(mat)
  col_annotation <- HeatmapAnnotation(
    df = df,
    col = list(group = color.use),
    which = "column",
    show_legend = FALSE,
    show_annotation_name = FALSE,
    simple_anno_size = grid::unit(0.2, "cm")
  )
  ha2 <- HeatmapAnnotation(
    Strength = anno_barplot(
      colSums(mat.ori),
      border = FALSE,
      gp = gpar(fill = color.use, col = color.use),
      axis_param = list(
        at = c(
          0,
          signif(max(colSums(mat.ori)), 1) / 2,
          signif(max(colSums(mat.ori)), 1)
        ),
        gp = gpar(fontsize = 10)
      )
    ),
    show_annotation_name = FALSE
  )
  pSum <- rowSums(mat.ori)
  pSum.original <- pSum
  pSum <- -1 / log(pSum)
  pSum[is.na(pSum)] <- 0
  idx1 <- which(is.infinite(pSum) | pSum < 0)
  if (length(idx1) > 0) {
    values.assign <- seq(
      max(pSum) * 1.1, max(pSum) * 1.5, length.out = length(idx1)
    )
    position <- sort(pSum.original[idx1], index.return = TRUE)$ix
    pSum[idx1] <- values.assign[match(1:length(idx1), position)]
  }
  ha1 <- rowAnnotation(
    Strength = anno_barplot(
      pSum,
      border = FALSE,
      axis_param = list(
        at = c(
          0,
          signif(max(pSum), 1) / 2,
          signif(max(pSum), 1)
        ),
        gp = gpar(fontsize = 10)
      )
    ),
    show_annotation_name = FALSE
  )
  if (min(mat, na.rm = T) == max(mat, na.rm = T)) {
    legend.break <- max(mat, na.rm = T)
  } else {
    legend.break <- c(
      round(min(mat, na.rm = T), digits = 1),
      round(max(mat, na.rm = T), digits = 1)
    )
  }
  ht1 <- Heatmap(
    mat,
    col = color.heatmap.use,
    na_col = "white",
    name = "Relative strength",
    bottom_annotation = col_annotation,
    top_annotation = ha2,
    right_annotation = ha1,
    cluster_rows = cluster.rows,
    cluster_columns = cluster.rows,
    row_names_side = "left",
    row_names_rot = 0,
    row_names_gp = gpar(fontsize = font.size),
    column_names_gp = gpar(fontsize = font.size),
    width = unit(width, "cm"),
    height = unit(height, "cm"),
    column_title = title,
    column_title_gp = gpar(fontsize = font.size.title),
    column_names_rot = 90,
    heatmap_legend_param = list(
      title_gp = gpar(fontsize = 14, fontface = "plain"),
      title_position = "leftcenter-rot",
      border = NA,
      at = legend.break,
      legend_height = unit(40, "mm"),
      labels_gp = gpar(fontsize = 18),
      grid_width = unit(2, "mm")
    )
  )
  #  draw(ht1)
  return(ht1)
}

png(
  paste0(
    path_results,
    "images/test.png"
  ),
  width = 800
)
netAnalysis_signalingRole_heatmap_custom(
    cc_list[[1]],
    signaling = c("VEGF", "BMP", "PTN", "GDF", "GRN", "HGF"),
    pattern = "all",
    width = 10,
    height = 8,
    title = names(cc_list)[1],
    color.heatmap = "OrRd",
    font.size = 16,
    font.size.title = 18
  )
dev.off()

## Custom CellChat netVisual_heatmap functions #####

netVisual_heatmap_custom <- function(
  object,
  comparison = c(1, 2),
  measure = c("count", "weight"),
  signaling = NULL,
  slot.name = c("netP", "net"),
  color.use = NULL,
  color.heatmap = c("#2166ac","#b2182b"),
  title.name = NULL, width = NULL,
  height = NULL,
  font.size = 8,
  font.size.title = 10,
  cluster.rows = FALSE,
  cluster.cols = FALSE,
  sources.use = NULL,
  targets.use = NULL,
  remove.isolate = FALSE,
  row.show = NULL,
  col.show = NULL
) {
  if (!is.null(measure)) {
    measure <- match.arg(measure)
  }
  slot.name <- match.arg(slot.name)
  if (is.list(object@net[[1]])) {
    message("Do heatmap based on a merged object \n")
    obj1 <- object@net[[comparison[1]]][[measure]]
    obj2 <- object@net[[comparison[2]]][[measure]]
    net.diff <- obj2 - obj1
    if (measure == "count") {
      if (is.null(title.name)) {
        title.name <= "Differential number of interactions"
      }
    } else if (measure == "weight") {
      if (is.null(title.name)) {
        title.name <- "Differential interaction strength"
      }
    }
    legend.name <- "Relative values"
  } else {
    message("Do heatmap based on a single object \n")
    if (!is.null(signaling)) {
      net.diff <- slot(object, slot.name)$prob[,,signaling]
      if (is.null(title.name)) {
        title.name <- paste0(signaling, " signaling network")
      }
      legend.name <- "Communication Prob."
    } else if (!is.null(measure)) {
      net.diff <- object@net[[measure]]
      if (measure == "count") {
        if (is.null(title.name)) {
          title.name <- "Number of interactions"
        }
      } else if (measure == "weight") {
        if (is.null(title.name)) {
          title.name <- "Interaction strength"
        }
      }
      legend.name <- title.name
    }
  }
  net <- net.diff
  if ((!is.null(sources.use)) | (!is.null(targets.use))) {
    df.net <- reshape2::melt(net, value.name = "value")
    colnames(df.net)[1:2] <- c("source", "target")
    # keep the interactions associated with sources and targets of interest
    if (!is.null(sources.use)){
      if (is.numeric(sources.use)) {
        sources.use <- rownames(net.diff)[sources.use]
      }
      df.net <- subset(df.net, source %in% sources.use)
    }
    if (!is.null(targets.use)){
      if (is.numeric(targets.use)) {
        targets.use <- rownames(net.diff)[targets.use]
      }
      df.net <- subset(df.net, target %in% targets.use)
    }
    cells.level <- rownames(net.diff)
    df.net$source <- factor(df.net$source, levels = cells.level)
    df.net$target <- factor(df.net$target, levels = cells.level)
    df.net$value[is.na(df.net$value)] <- 0
    net <- tapply(
      df.net[["value"]],
      list(df.net[["source"]], df.net[["target"]]),
      sum
    )
  }
  net[is.na(net)] <- 0
  if (remove.isolate) {
    idx1 <- which(Matrix::rowSums(net) == 0)
    idx2 <- which(Matrix::colSums(net) == 0)
    idx <- intersect(idx1, idx2)
    net <- net[-idx, ]
    net <- net[, -idx]
    if (length(idx) > 0) {
      net <- net[-idx, ]
      net <- net[, -idx]
    }
  }
  mat <- net
  if (is.null(color.use)) {
    color.use <- scPalette(ncol(mat))
  }
  names(color.use) <- colnames(mat)
  if (!is.null(row.show)) {
    mat <- mat[row.show, ]
  }
  if (!is.null(col.show)) {
    mat <- mat[ ,col.show]
    color.use <- color.use[col.show]
  }
  if (min(mat) < 0) {
    color.heatmap.use <- colorRamp3(
      c(min(mat), 0, max(mat)),
      c(color.heatmap[1], "#f7f7f7", color.heatmap[2])
    )
    colorbar.break <- c(
      round(
        min(mat, na.rm = T),
        digits = nchar(sub(".*\\.(0*).*","\\1",
        min(mat, na.rm = T))) + 1
      ),
      0,
      round(
        max(mat, na.rm = T),
        digits = nchar(sub(".*\\.(0*).*","\\1",max(mat, na.rm = T))) + 1
        )
      )
  } else {
    if (length(color.heatmap) == 3) {
      color.heatmap.use = colorRamp3(c(0, min(mat), max(mat)), color.heatmap)
    } else if (length(color.heatmap) == 2) {
      color.heatmap.use = colorRamp3(c(min(mat), max(mat)), color.heatmap)
    } else if (length(color.heatmap) == 1) {
      color.heatmap.use = grDevices::colorRampPalette((
        RColorBrewer::brewer.pal(n = 9, name = color.heatmap
      )))(100)
    }
    colorbar.break <- c(
      round(
        min(mat, na.rm = T),
        digits = nchar(sub(".*\\.(0*).*","\\1",min(mat, na.rm = T))) + 1
      ),
      round(
        max(mat, na.rm = T),
        digits = nchar(sub(".*\\.(0*).*","\\1",
        max(mat, na.rm = T))) + 1
      )
    )
  }
  df <- data.frame(group = colnames(mat))
  rownames(df) <- colnames(mat)
  col_annotation <- HeatmapAnnotation(
    df = df,
    col = list(group = color.use),
    which = "column",
    show_legend = FALSE,
    show_annotation_name = FALSE,
    simple_anno_size = grid::unit(0.2, "cm")
  )
  row_annotation <- HeatmapAnnotation(
    df = df,
    col = list(group = color.use),
    which = "row",
    show_legend = FALSE,
    show_annotation_name = FALSE,
    simple_anno_size = grid::unit(0.2, "cm")
  )
  ha1 <- rowAnnotation(
    Strength = anno_barplot(
      rowSums(abs(mat)),
      border = FALSE,
      gp = gpar(fill = color.use, col = color.use),
      axis_param = list(
        at = c(
          0,
          signif(max(rowSums(abs(mat))), 1) / 2,
          signif(max(rowSums(abs(mat))), 1)
        ),
        gp = gpar(fontsize = 10)
      )
    ),
    show_annotation_name = FALSE
  )
  ha2 <- HeatmapAnnotation(
    Strength = anno_barplot(
      colSums(abs(mat)),
      border = FALSE,
      gp = gpar(fill = color.use, col = color.use),
      axis_param = list(
        at = c(
          0,
          signif(max(colSums(abs(mat))), 1) / 2,
          signif(max(colSums(abs(mat))), 1)
        ),
        gp = gpar(fontsize = 10)
      )
    ),
    show_annotation_name = FALSE
  )
  if (sum(abs(mat) > 0) == 1) {
    color.heatmap.use = c("white", color.heatmap.use)
  } else {
    mat[mat == 0] <- NA
  }
  ht1 = Heatmap(
    mat,
    col = color.heatmap.use,
    na_col = "white",
    name = legend.name,
    bottom_annotation = col_annotation,
    left_annotation = row_annotation,
    top_annotation = ha2,
    right_annotation = ha1,
    cluster_rows = cluster.rows,
    cluster_columns = cluster.rows,
    row_names_side = "left",
    row_names_rot = 0,
    row_names_gp = gpar(fontsize = font.size),
    column_names_gp = gpar(fontsize = font.size),
    # width = unit(width, "cm"), height = unit(height, "cm"),
    column_title = title.name,
    column_title_gp = gpar(fontsize = font.size.title),
    column_names_rot = 90,
    row_title = "Sources (Sender)",
    row_title_gp = gpar(fontsize = font.size.title),
    row_title_rot = 90,
    heatmap_legend_param = list(
      title_gp = gpar(fontsize = 14, fontface = "plain"),
      title_position = "leftcenter-rot",
      border = NA, #at = colorbar.break,
      legend_height = unit(40, "mm"),
      labels_gp = gpar(fontsize = 18),
      grid_width = unit(2, "mm")
    )
  )
  #  draw(ht1)
  return(ht1)
}

png(
  paste0(
    path_results,
    "images/test2.png"
  ),
  width = 800
)
netVisual_heatmap_custom(
  cellchat_healthy,
  signaling = "VEGF",
  color.heatmap = "Reds",
  title.name = paste("VEGF signaling ", names(cc_list)[1]),
  width = 10,
  height = 8,
  font.size = 16,
  font.size.title = 18
)
dev.off()

## Custom CellChat netVisual_aggregate functions #####

searchPair_custom <- function(
  signaling = c(),
  pairLR.use,
  key = c("pathway_name","ligand"),
  matching.exact = FALSE, pair.only = TRUE
) {
  key <- match.arg(key)
  pairLR = future.apply::future_sapply(
    X = 1:length(signaling),
    FUN = function(x) {
      if (!matching.exact) {
        index <- grep(signaling[x], pairLR.use[[key]])
      } else {
        index <- which(pairLR.use[[key]] %in% signaling[x])
      }
      if (length(index) > 0) {
        if (pair.only) {
          pairLR <- dplyr::select(pairLR.use[index, ], interaction_name, pathway_name, ligand, receptor)
        } else {
          pairLR <- pairLR.use[index, ]
        }
        return(pairLR)
      } else {
        stop(cat(paste("Cannot find ", signaling[x], ".", "Please input a correct name!"),'\n'))
      }
    }
  )
  if (pair.only) {
    pairLR0 <- vector("list", length(signaling))
    for (i in 1:length(signaling)) {
      pairLR0[[i]] <- matrix(unlist(pairLR[c(4*i-3, 4*i-2, 4*i-1, 4*i)]), ncol=4, byrow=F)
    }
    pairLR <- do.call(rbind, pairLR0)
    dimnames(pairLR)[[2]] <- dimnames(pairLR.use)[[2]][1:4]
    rownames(pairLR) <- pairLR[,1]
  } else {
    pairLR0 <- vector("list", length(signaling))
    for (i in 1:length(signaling)) {
      pairLR0[[i]] <- matrix(unlist(pairLR[(i*ncol(pairLR.use)-(ncol(pairLR.use)-1)):(i*ncol(pairLR.use))]), ncol=ncol(pairLR.use), byrow=F)
    }
    pairLR <- do.call(rbind, pairLR0)
    dimnames(pairLR)[[2]] <- dimnames(pairLR.use)[[2]]
    rownames(pairLR) <- pairLR[,1]
  }
  return(as.data.frame(pairLR, stringsAsFactors = FALSE))
}

netVisual_circle_custom <- function(
  net,
  color.use = NULL,
  title.name = NULL,
  sources.use = NULL,
  targets.use = NULL,
  idents.use = NULL,
  remove.isolate = FALSE,
  top = 1,
  weight.scale = FALSE,
  vertex.weight = 20,
  vertex.weight.max = NULL,
  vertex.size.max = NULL,
  vertex.label.cex = 1,
  vertex.label.color = "black",
  edge.weight.max = NULL,
  edge.width.max=8,
  alpha.edge = 0.6,
  label.edge = FALSE,
  edge.label.color = 'black',
  edge.label.cex = 0.8,
  edge.curved = 0.2,
  shape = 'circle',
  layout = in_circle(),
  margin = 0.4,
  vertex.size = NULL,
  arrow.width=1,
  arrow.size = 0.2
) {
  if (!is.null(vertex.size)) {
    warning("'vertex.size' is deprecated. Use `vertex.weight`")
  }
  if (is.null(vertex.size.max)) {
    if (length(unique(vertex.weight)) == 1) {
      vertex.size.max <- 5
    } else {
      vertex.size.max <- 15
    }
  }
  options(warn = -1)
  thresh <- stats::quantile(net, probs = 1-top)
  net[net < thresh] <- 0
  if ((!is.null(sources.use)) | (!is.null(targets.use)) | (!is.null(idents.use)) ) {
    if (is.null(rownames(net))) {
      stop("The input weighted matrix should have rownames!")
    }
    cells.level <- rownames(net)
    df.net <- reshape2::melt(net, value.name = "value")
    colnames(df.net)[1:2] <- c("source","target")
    # keep the interactions associated with sources and targets of interest
    if (!is.null(sources.use)){
      if (is.numeric(sources.use)) {
        sources.use <- cells.level[sources.use]
      }
      df.net <- subset(df.net, source %in% sources.use)
    }
    if (!is.null(targets.use)){
      if (is.numeric(targets.use)) {
        targets.use <- cells.level[targets.use]
      }
      df.net <- subset(df.net, target %in% targets.use)
    }
    if (!is.null(idents.use)) {
      if (is.numeric(idents.use)) {
        idents.use <- cells.level[idents.use]
      }
      df.net <- filter(df.net, (source %in% idents.use) | (target %in% idents.use))
    }
    df.net$source <- factor(df.net$source, levels = cells.level)
    df.net$target <- factor(df.net$target, levels = cells.level)
    df.net$value[is.na(df.net$value)] <- 0
    net <- tapply(df.net[["value"]], list(df.net[["source"]], df.net[["target"]]), sum)
  }
  net[is.na(net)] <- 0
  if (remove.isolate) {
    idx1 <- which(Matrix::rowSums(net) == 0)
    idx2 <- which(Matrix::colSums(net) == 0)
    idx <- intersect(idx1, idx2)
    net <- net[-idx, ]
    net <- net[, -idx]
  }
  g <- graph_from_adjacency_matrix(net, mode = "directed", weighted = T)
  edge.start <- igraph::ends(g, es=igraph::E(g), names=FALSE)
  coords<-layout_(g,layout)
  if(nrow(coords)!=1){
    coords_scale=scale(coords)
  }else{
    coords_scale<-coords
  }
  if (is.null(color.use)) {
    color.use = scPalette(length(igraph::V(g)))
  }
  if (is.null(vertex.weight.max)) {
    vertex.weight.max <- max(vertex.weight)
  }
  vertex.weight <- vertex.weight/vertex.weight.max*vertex.size.max+5
  loop.angle<-ifelse(coords_scale[igraph::V(g),1]>0,-atan(coords_scale[igraph::V(g),2]/coords_scale[igraph::V(g),1]),pi-atan(coords_scale[igraph::V(g),2]/coords_scale[igraph::V(g),1]))
  igraph::V(g)$size<-vertex.weight
  igraph::V(g)$color<-color.use[igraph::V(g)]
  igraph::V(g)$frame.color <- color.use[igraph::V(g)]
  igraph::V(g)$label.color <- vertex.label.color
  igraph::V(g)$label.cex<-vertex.label.cex
  if(label.edge){
    igraph::E(g)$label<-igraph::E(g)$weight
    igraph::E(g)$label <- round(igraph::E(g)$label, digits = 1)
  }
  if (is.null(edge.weight.max)) {
    edge.weight.max <- max(igraph::E(g)$weight)
  }
  if (weight.scale == TRUE) {
    #E(g)$width<-0.3+edge.width.max/(max(E(g)$weight)-min(E(g)$weight))*(E(g)$weight-min(E(g)$weight))
    igraph::E(g)$width<- 0.3+igraph::E(g)$weight/edge.weight.max*edge.width.max
  }else{
    igraph::E(g)$width<-0.3+edge.width.max*igraph::E(g)$weight
  }
  igraph::E(g)$arrow.width<-arrow.width
  igraph::E(g)$arrow.size<-arrow.size
  igraph::E(g)$label.color<-edge.label.color
  igraph::E(g)$label.cex<-edge.label.cex
  igraph::E(g)$color<- grDevices::adjustcolor(igraph::V(g)$color[edge.start[,1]],alpha.edge)
  if(sum(edge.start[,2]==edge.start[,1])!=0){
    igraph::E(g)$loop.angle[which(edge.start[,2]==edge.start[,1])]<-loop.angle[edge.start[which(edge.start[,2]==edge.start[,1]),1]]
  }
  radian.rescale <- function(x, start=0, direction=1) {
    c.rotate <- function(x) (x + start) %% (2 * pi) * direction
    c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
  }
  label.locs <- radian.rescale(x=1:length(igraph::V(g)), direction=-1, start=0)
  label.dist <- vertex.weight/max(vertex.weight)+2
  plot(
    g,
    edge.curved = edge.curved,
    vertex.shape = shape,
    layout = coords_scale,
    margin = margin,
    vertex.label.dist = label.dist,
    vertex.label.degree = label.locs,
    vertex.label.family = "Helvetica",
    edge.label.family = "Helvetica"
  ) # "sans"
  if (!is.null(title.name)) {
    text(0, 1.5, title.name, cex = 1.5)
  }
  # https://www.andrewheiss.com/blog/2016/12/08/save-base-graphics-as-pseudo-objects-in-r/
  # grid.echo()
  # gg <-  grid.grab()
  gg <- recordPlot()
  return(gg)
}

netVisual_aggregate_custom <- function(
  object,
  signaling,
  signaling.name = NULL,
  color.use = NULL,
  vertex.receiver = NULL,
  sources.use = NULL,
  targets.use = NULL,
  idents.use = NULL,
  top = 1,
  remove.isolate = FALSE,
  vertex.weight = 1,
  vertex.weight.max = NULL,
  vertex.size.max = NULL,
  weight.scale = TRUE,
  edge.weight.max = NULL,
  edge.width.max=8,
  thresh = 0.05,
  from = NULL,
  to = NULL,
  bidirection = NULL,
  pt.title = 12,
  title.space = 6,
  vertex.label.cex = 0.8,
  group = NULL,
  cell.order = NULL,
  small.gap = 1,
  big.gap = 10,
  scale = FALSE,
  reduce = -1,
  show.legend = FALSE,
  legend.pos.x = 20,
  legend.pos.y = 20,
  ...
) {
  if (is.null(vertex.weight)) {
    vertex.weight <- as.numeric(table(object@idents))
  }
  if (is.null(vertex.size.max)) {
    if (length(unique(vertex.weight)) == 1) {
      vertex.size.max <- 5
    } else {
      vertex.size.max <- 15
    }
  }
  pairLR <- searchPair_custom(
    signaling = signaling,
    pairLR.use = object@LR$LRsig,
    key = "pathway_name",
    matching.exact = T, pair.only = T
  )
  if (is.null(signaling.name)) {
    signaling.name <- signaling
  }
  net <- object@net
  pairLR.use.name <- dimnames(net$prob)[[3]]
  pairLR.name <- intersect(rownames(pairLR), pairLR.use.name)
  pairLR <- pairLR[pairLR.name, ]
  prob <- net$prob
  pval <- net$pval
  prob[pval > thresh] <- 0
  if (length(pairLR.name) > 1) {
    pairLR.name.use <- pairLR.name[apply(prob[,,pairLR.name], 3, sum) != 0]
  } else {
    pairLR.name.use <- pairLR.name[sum(prob[,,pairLR.name]) != 0]
  }
  if (length(pairLR.name.use) == 0) {
    stop(paste0('There is no significant communication of ', signaling.name))
  } else {
    pairLR <- pairLR[pairLR.name.use,]
  }
  nRow <- length(pairLR.name.use)
  prob <- prob[, ,pairLR.name.use]
  pval <- pval[, ,pairLR.name.use]
  if (length(dim(prob)) == 2) {
    prob <- replicate(1, prob, simplify = "array")
    pval <- replicate(1, pval, simplify = "array")
  }
  prob.sum <- apply(prob, c(1,2), sum)
  gg <- netVisual_circle_custom(
    prob.sum,
    sources.use = sources.use,
    targets.use = targets.use,
    idents.use = idents.use,
    remove.isolate = remove.isolate,
    top = top,
    color.use = color.use,
    vertex.weight = vertex.weight,
    vertex.weight.max = vertex.weight.max,
    vertex.size.max = vertex.size.max,
    weight.scale = weight.scale,
    edge.weight.max = edge.weight.max,
    edge.width.max = edge.width.max,
    title.name = paste0(
      signaling.name,
      " signaling pathway network"
    ),
    vertex.label.cex = vertex.label.cex,
    ...
  )
  return(gg)
}

png(
  paste0(
    path_results,
    "images/test3.png"
  ),
  width = 500
)
netVisual_aggregate_custom(
  cellchat_sen,
  signaling = "TENASCIN",
  vertex.label.cex = 1.2
)
dev.off()


## Custom CellChat netAnalysis_signalingRole_network functions #####

netAnalysis_signalingRole_network_custom <- function(
  object,
  signaling,
  slot.name = "netP",
  measure = c("outdeg","indeg","flowbet","info"),
  measure.name = c("Sender","Receiver","Mediator","Influencer"),
  color.use = NULL,
  color.heatmap = "BuGn",
  width = 6.5,
  height = 1.4,
  font.size = 8,
  font.size.title = 10,
  cluster.rows = FALSE,
  cluster.cols = FALSE
) {
  if (length(slot(object, slot.name)$centr) == 0) {
    stop("Please run `netAnalysis_computeCentrality`")
  }
  centr <- slot(object, slot.name)$centr[signaling]
  for(i in 1:length(centr)) {
    centr0 <- centr[[i]]
    mat <- matrix(unlist(centr0), ncol = length(centr0), byrow = FALSE)
    mat <- t(mat)
    rownames(mat) <- names(centr0); colnames(mat) <- names(centr0$outdeg)
    if (!is.null(measure)) {
      mat <- mat[measure,]
      if (!is.null(measure.name)) {
        rownames(mat) <- measure.name
      }
    }
    mat <- sweep(mat, 1L, apply(mat, 1, max), '/', check.margin = FALSE)
    if (is.null(color.use)) {
      color.use <- scPalette(length(colnames(mat)))
    }
    color.heatmap.use = grDevices::colorRampPalette((
      RColorBrewer::brewer.pal(n = 9, name = color.heatmap)
      ))(100)
    df <- data.frame(group = colnames(mat))
    rownames(df) <- colnames(mat)
    cell.cols.assigned <- setNames(color.use, unique(as.character(df$group)))
    col_annotation <- HeatmapAnnotation(
      df = df,
      col = list(group = cell.cols.assigned),
      which = "column",
      show_legend = FALSE,
      show_annotation_name = FALSE,
      simple_anno_size = grid::unit(0.2, "cm")
    )
    ht1 = Heatmap(
      mat,
      col = color.heatmap.use,
      na_col = "white",
      name = "Importance",
      bottom_annotation = col_annotation,
      cluster_rows = cluster.rows,
      cluster_columns = cluster.rows,
      row_names_side = "left",
      row_names_rot = 0,
      row_names_gp = gpar(fontsize = font.size),
      column_names_gp = gpar(fontsize = font.size),
      width = unit(width, "cm"), height = unit(height, "cm"),
      column_title = paste0(
        names(centr[i]),
        " signaling pathway network"
      ),
      column_title_gp = gpar(fontsize = font.size.title),
      column_names_rot = 45,
      heatmap_legend_param = list(
        title = "Importance",
        title_gp = gpar(fontsize = 14, fontface = "plain"),
        title_position = "leftcenter-rot",
        border = NA,
        at = c(
          round(min(mat, na.rm = T),
          digits = 1),
          round(max(mat, na.rm = T), digits = 1)
          ),
        legend_height = unit(35, "mm"),
        labels_gp = gpar(fontsize = 16),
        grid_width = unit(2, "mm")
      )
    )
    draw(ht1)
  }
}
