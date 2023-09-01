library(CellTrails)
library(tidyverse)
library(reshape)
library(pheatmap)
library(RColorBrewer)

#' Fitting dynamics 
#' Fits expression as a function of pseudotime using generalized additive models
#' @param x Pseudotime values
#' @param y Expression values
#' @param z Sample states
#' @param k Dimensions of the bases used to represent the smooth term
#' @param n.out Number of predicted expression values within pseudotime range
#' @param x.out Predicted expression values for given set of pseudotime values
#' @return A list containing the following components:
#' \item{\code{x}}{Pseudotime}
#' \item{\code{y}}{Fitted values}
#' \item{\code{mod}}{GAM}
#' @importFrom mgcv gam
fitDynamic_def <- function(x, y, z, k=5, n.out=NULL, x.out=NULL) {
  x.ptime <- x
  y.expr <- y
  z.states <- z
  
  if(var(y.expr) == 0) {
    #return(list(x = NA, y = NA, mod = NA))
    if(is.null(x.out)) {
      x.out <- x.ptime
    }
    return(list(x=x.out, y=rep(mean(y.expr), length(x.out)), mod=NA))
  }
  
  # dynamic nd weight
  w <- aggregate(y.expr, list(z.states), function(x) {sum(x == 0) / length(x)})
  w <- w[match(z.states, w[,1]), 2]
  w[y.expr > 0] <- 1 #set measured value weight always to 1
  
  y <- y.expr
  x <- x.ptime
  
  eq <- formula(paste0("y ~ te(x, k = ", k, ", bs = \"tp\", fx=TRUE)"))
  mod <- mgcv::gam(eq, family="gaussian", weights=w, select=TRUE, method="REML",
                   gamma=1)
  
  y.pred <- NA
  if(!is.null(n.out)) {
    s <- seq(from=0, to=1, length.out=n.out)
    y.pred <- stats::predict(mod, type = "response", newdata = data.frame(x = s))
    x.ptime <- s
  } else {
    y.pred <- predict(mod, type = "response")
  }
  
  y.pred[y.pred < 0] <- 0
  o <- order(x.ptime)
  list(x = x.ptime[o], y = y.pred[o], mod = mod)
}

## Function to get smoothed expression of in vivo trails
#' Fitting dynamics
#' Fits expression as a function of pseudotime using generalized
#' additive models
#' @param sce A SingleCellExperiment object containing expression matrix
#' @param time_df Data frame of all trail pseudotimes from sce object along with a column of cell state
#' @param gene_list Vector of all dynamically expressed genes along the trail
#' @param trail_name Trail to be test
#' @return A matrix of smoothed expression
get_smoothed_expr <- function(sce, time_df, gene_list, trail_name) {
  smoothed_expr = NULL
  trail_df <- arrange(time_df[!is.na(time_df[[trail_name]]),],trail_name)
  x.ptime <- trail_df[[trail_name]]
  x <- x.ptime/max(x.ptime)
  name = rownames(trail_df)
  num = as.integer(vapply(strsplit(name, "_", fixed = TRUE), "[", "", 2))
  z = trail_df$state
  
  for (i in gene_list){
    y <- sce@assays@data@listData[["logcounts"]][i,num]
    fit <- fitDynamic_def(x,y,z,n.out = length(x))
    smoothed_expr = rbind(smoothed_expr,fit$y)
  }
  rownames(smoothed_expr) = gene_list
  colnames(smoothed_expr) = fit$x
  return (smoothed_expr)
}

## Function to get smoothed expression of in vitro trails
#' Fitting dynamics
#' Fits expression as a function of pseudotime using generalized
#' additive models
#' @param sce A SingleCellExperiment object containing expression matrix
#' @param time_df Data frame of all trail pseudotimes from sce object along with a column of cell state
#' @param gene_list Vector of all dynamically expressed genes along the trail
#' @param trail_name Trail to be test
#' @return A matrix of smoothed expression
get_invitro_smoothed_expr <- function(sce, time_df, gene_list, trail_name) {
  smoothed_expr = NULL
  trail_df <- arrange(time_df[!is.na(time_df[[trail_name]]),],trail_name)
  x.ptime <- trail_df[[trail_name]]
  x <- x.ptime/max(x.ptime)
  name = rownames(trail_df)
  z = trail_df$state
  
  for (i in gene_list){
    y <- sce@assays@data@listData[["logcounts"]][i,name]
    fit <- fitDynamic_def(x,y,z,n.out = length(x))
    smoothed_expr = rbind(smoothed_expr,fit$y)
  }
  rownames(smoothed_expr) = gene_list
  colnames(smoothed_expr) = fit$x
  return (smoothed_expr)
}

# Import in vitro reclustering object
HC_SC_sub <- readRDS("~/Invitro_HC_SC_rmS8_11states_sce.rds")

# Define subtrails for in vitro cells
HC_SC_sub <- addTrail(HC_SC_sub, from="B5", to="H8", name="trailA")
HC_SC_sub <- addTrail(HC_SC_sub, from="B5", to="H10", name="trailB")
invitro_time_df <- as.data.frame(trails(HC_SC_sub))

Jan_sce <- readRDS("~/external_data/Jan_2021_sce_wt.mouse.utricle_805_cells.rds")

Jan_sce <- connectStates(Jan_obj, l=7)
Jan_sce <- fitTrajectory(Jan_obj)

tl <- read.ygraphml(file="~/celltrails/Jan_12states_layout.graphml")
trajLayout(Jan_sce, adjust=TRUE) <- tl

Jan_sce <- addTrail(Jan_sce, from="U1", to="H3", name="trailA")
Jan_sce <- addTrail(Jan_sce, from="U1", to="H5", name="trailB")
invivo_time_df <- as.data.frame(trails(Jan_sce))

trailA_gene <- read.table("~/Jan_trailA_dynamic_gene.txt",header=TRUE)
trailB_gene <- read.table("~/Jan_trailB_dynamic_gene.txt",header=TRUE)

# In vivo trailA pseudotime heatmap
smoothed_expr = get_smoothed_expr(Jan_sce,invivo_time_df,trailA_gene$x, 'trailA')
heatdata <- as.matrix(smoothed_expr)
palette.breaks <- seq(-3, 3, 0.01)
pheatmap(t(scale(t(heatdata))), cluster_rows=TRUE, cluster_cols=FALSE, scale = 'none',show_rownames=FALSE,show_colnames = FALSE,legend = TRUE, clustering_distance_rows = "correlation", clustering_method = "ward.D", breaks = palette.breaks, color = colorRampPalette(rev(brewer.pal(n = 9, name ="RdYlBu")))(length(palette.breaks)-1),treeheight_row=0)

# In vitro trailA pseudotime heatmap
invivo_trailA_gene = rownames(t(scale(t(heatdata)))[out$tree_row[["order"]],])
invitro_gene = data.frame(x = invivo_trailA_gene) %>% inner_join(data.frame(x = as.vector(featureNames(HC_SC_sub)))) %>% pull(x)
invitro_smoothed_expr = get_invitro_smoothed_expr(HC_SC_sub,invitro_time_df,invitro_gene, 'trailA')
heatdata <- as.matrix(invitro_smoothed_expr)
pheatmap(t(scale(t(heatdata))), cluster_rows=FALSE, cluster_cols=FALSE, scale = 'none',show_rownames=FALSE,show_colnames = FALSE,legend = TRUE, breaks = palette.breaks, color = colorRampPalette(rev(brewer.pal(n = 9, name ="RdYlBu")))(length(palette.breaks)-1),treeheight_row=0)

# In vivo trailB pseudotime heatmap
smoothed_expr = get_smoothed_expr(Jan_sce,invivo_time_df,trailB_gene$x, 'trailB')

heatdata <- as.matrix(smoothed_expr)
pheatmap(t(scale(t(heatdata))), cluster_rows=TRUE, cluster_cols=FALSE, scale = 'none',show_rownames=FALSE,show_colnames = FALSE,legend = TRUE, clustering_distance_rows = "correlation", clustering_method = "ward.D", breaks = palette.breaks, color = colorRampPalette(rev(brewer.pal(n = 9, name ="RdYlBu")))(length(palette.breaks)-1),treeheight_row=0, cutree_rows = 6)

# In vitro trailB pseudotime heatmap
invivo_trailB_gene = rownames(t(scale(t(heatdata)))[out$tree_row[["order"]],])
invitro_gene = data.frame(x = invivo_trailB_gene) %>% inner_join(data.frame(x = as.vector(featureNames(HC_SC_sub)))) %>% pull(x)
invitro_smoothed_expr = get_invitro_smoothed_expr(HC_SC_sub,invitro_time_df,invitro_gene, 'trailB')

heatdata <- as.matrix(invitro_smoothed_expr)
pheatmap(t(scale(t(heatdata))), cluster_rows=FALSE, cluster_cols=FALSE, scale = 'none',show_rownames=FALSE,show_colnames = FALSE,legend = TRUE, breaks = palette.breaks, color = colorRampPalette(rev(brewer.pal(n = 9, name ="RdYlBu")))(length(palette.breaks)-1),treeheight_row=0)


