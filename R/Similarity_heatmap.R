library(dplyr)
library(Seurat)
library(ggplot2)

# Import in vitro reclustering object
HC_SC_sub <- readRDS("~/Invitro_HC_SC_rmS8_11states_sce.rds")

# Define subtrails for in vitro cells
plotMap(HC_SC_sub, color_by="phenoName", name="landmark")
HC_SC_sub <- addTrail(HC_SC_sub, from="B5", to="H8", name="trailA")
HC_SC_sub <- addTrail(HC_SC_sub, from="B5", to="H10", name="trailB")

# Import Jan data 805 cells
Jan_sce <- readRDS("~/external_data/Jan_2021_sce_wt.mouse.utricle_805_cells.rds")

# Rerun Celltrails for Jan data
Jan_sce <- connectStates(Jan_sce, l=7)
Jan_sce <- fitTrajectory(Jan_sce)
write.ygraphml(Jan_sce, file="~/celltrails/Jan_12states.graphml", color_by="phenoName",name="state", node_label="state")
tl <- read.ygraphml(file="~/celltrails/Jan_12states_layout.graphml")
trajLayout(Jan_sce, adjust=TRUE) <- tl
plotMap(Jan_sce, color_by="phenoName", name="state")

# Define subtrails
Jan_sce <- addTrail(Jan_sce, from="U1", to="H3", name="trailA")
Jan_sce <- addTrail(Jan_sce, from="U1", to="H5", name="trailB")
Jan_sce <- addTrail(Jan_sce, from="U1", to="H7", name="trailC")

plotTrail(Jan_sce, name="trailA")
plotTrail(Jan_sce, name="trailB")
plotTrail(Jan_sce, name="trailC")

plotDynamic(Jan_sce, feature_name=c("Sox2", "Atoh1","Tmc1"), trail_name="trailA")
plotDynamic(Jan_sce, feature_name=c("Sox2", "Atoh1","Tmc1"), trail_name="trailB")

# 
invitro_time_df <- as.data.frame(trails(HC_SC_sub))
invitro_time_df$state = HC_SC_sub@colData@listData[["CellTrails.state"]]
invivo_time_df <- as.data.frame(trails(Jan_sce))
invivo_time_df$state = Jan_sce@colData@listData[["CellTrails.state"]]

expr_df <- as.data.frame(Jan_sce@assays@data@listData[["logcounts"]])
Jan_object <- as.Seurat(Jan_sce, counts = "counts", data = "logcounts")

#' find_top_variable() is to fine top 2000 variable features in vitro and in vivo datasets
#' @param invitro_df A dataframe of in vitro trails
#' @param invivo_df A dataframe of in vivo trails
#' @param invitro_obj A Seurat object of in vitro HCs and SCs
#' @param invivo_obj A Seurat object of in vivo HCs and SCs
#' @param trail Name of trails
find_top_variable = function(invitro_df,invivo_df,invitro_obj,invivo_obj,trail){
  invivo_cell <- invivo_time_df[!is.na(invivo_time_df[[trail]]),] %>% arrange(!!trail) %>% rownames()
  invivo_cell_name <- as.integer(vapply(strsplit(invivo_cell, "_", fixed = TRUE), "[", "", 2))
  invitro_cell <- invitro_df[!is.na(invitro_df[[trail]]),] %>% arrange(!!trail) %>% rownames()
  invitro_trail <- subset(invitro_obj,cells = cell)
  invivo_trail <- subset(invivo_obj,cells = invivo_cell_name)
  invivo_trail <- FindVariableFeatures(invivo_trail, selection.method = "vst", nfeatures = 2000)
  
  list <- list(invitro_trail,invivo_trail)
  
  all_datasets.anchors <- FindIntegrationAnchors(object.list = list, dims = 1:15, anchor.features = 2000, normalization.method = "LogNormalize")
  all_datasets.combined <- IntegrateData(anchorset = all_datasets.anchors, dims = 1:15, normalization.method = "LogNormalize")
  DefaultAssay(all_datasets.combined) <- "integrated"
  all_genes1 = rownames(all_datasets.combined)
  all_datasets.combined <- ScaleData(all_datasets.combined, features = all_genes1)
  all_datasets.combined <- RunPCA(all_datasets.combined, npcs = 30, verbose = FALSE)
  top_variables <- SelectIntegrationFeatures(object.list = list, nfeatures = 2000)
  
  return(list(top_variables,invitro_cell,invivo_cell,invivo_cell_name))
}

# TrailA in vitro and in vivo Cosine distance matrix
result <- find_top_variable(invitro_time_df,invivo_time_df,invitro_HC_SC,Jan_object,'trailA')
top_variables <- result[[1]]
invitro_cell <- result[[2]]
invivo_cell <- result[[3]]
invivo_cell_name <- result[[4]]
invitro_trailA_df <- as.data.frame(t(HC_SC_sub@assays@data@listData[["logcounts"]][rownames(HC_SC_sub@assays@data@listData[["logcounts"]]) %in% top_variables,invitro_cell]))
order <- c(1:dim(invitro_trailA_df)[1])
aggr_invitro_df <- generate_meta_cell(invitro_trailA_df,order,ncells=6)

invivo_trailA_df <- as.data.frame(t(expr_df[rownames(expr_df) %in% top_variables,invivo_cell_name]))
order <- c(1:dim(invivo_trailA_df)[1])
aggr_invivo_df <- generate_meta_cell(invivo_trailA_df,order,ncells=6)
aggr_invivo_df <- aggr_invivo_df[,colnames(aggr_invitro_df)]

similarity_mtx <- calculate_cell_similarity(aggr_invitro_df,aggr_invivo_df, method = "cosine")
distance_mtx <- distance_to_similarity(similarity_mtx,method = "cosine")
categories_col = data.frame(invitro_state = factor(add_label(invitro_time_df[invitro_cell,]$state,6)))
categories_row = data.frame(invivo_state = factor(add_label(invivo_time_df[invivo_cell,]$state,6)))

annoCol <- list(invitro_state = c(S2 = "#DF8747", S3="#B44439", S4="#866DB2", S5="#C683BB"),invivo_state = c(S5="#A75E30",S7="#82BB67",S8="#4D82B9"))
heatmap_visual(distance_mtx,normalize ="both",annotation_col = categories_col, annotation_row = categories_row, annotation_colors = annoCol, annotation_legend = F, legend = FALSE, color = col.pal)

# TrailA in vitro and in vivo Euclidean distance matrix
similarity_mtx <- calculate_cell_similarity(aggr_invitro_df,aggr_invivo_df, method = "Euclidean")
distance_mtx <- distance_to_similarity(similarity_mtx,method = "inverse")

heatmap_visual(distance_mtx,normalize ="both",annotation_col = categories_col, annotation_row = categories_row, annotation_colors = annoCol, annotation_legend = F, legend = FALSE, color = col.pal)

# TrailB in vitro and in vivo Cosine distance matrix
result <- find_top_variable(invitro_time_df,invivo_time_df,invitro_HC_SC,Jan_object,'trail2')
top_variables <- result[[1]]
invitro_cell <- result[[2]]
invivo_cell <- result[[3]]
invivo_cell_name <- result[[4]]
invitro_trailB_df <- as.data.frame(t(HC_SC_sub@assays@data@listData[["logcounts"]][rownames(HC_SC_sub@assays@data@listData[["logcounts"]]) %in% top_variables,invitro_cell]))
order <- c(1:dim(invitro_trailB_df)[1])
aggr_invitro_df <- generate_meta_cell(invitro_trailB_df,order,ncells=6)

invivo_trailB_df <- as.data.frame(t(expr_df[rownames(expr_df) %in% top_variables,invivo_cell_name]))
order <- c(1:dim(invivo_trailB_df)[1])
aggr_invivo_df <- generate_meta_cell(invivo_trailB_df,order,ncells=6)
aggr_invivo_df <- aggr_invivo_df[,colnames(aggr_invitro_df)]

similarity_mtx <- calculate_cell_similarity(aggr_invitro_df,aggr_invivo_df, method = "cosine")
distance_mtx <- distance_to_similarity(similarity_mtx,method = "cosine")
categories_col = data.frame(invitro_state = factor(add_label(invitro_time_df[invitro_cell,]$state,6)))
categories_row = data.frame(invivo_state = factor(add_label(invivo_time_df[invivo_cell,]$state,6)))

annoCol <- list(invitro_state = c(S2="#DF8747",S3="#B44439",S4="#866DB2",S5="#C683BB",S6="#7E7E7E",S7="#BBBA58",S8="#B6C4E1",S9="#EABD89",S11="#E4A09A"),invivo_state = c(S2="#AECDE0",S3="#F3E495",S4="#D3A086",S5="#A75E30",S6="#EEB678",S7="#82BB67",S8="#4D82B9"))
heatmap_visual(distance_mtx,normalize ="both",annotation_col = categories_col, annotation_row = categories_row, annotation_colors = annoCol, annotation_legend = F, legend = FALSE, color = col.pal)

# TrailB in vitro and in vivo Euclidean distance matrix
similarity_mtx <- calculate_cell_similarity(aggr_invitro_df,aggr_invivo_df, method = "Euclidean")
distance_mtx <- distance_to_similarity(similarity_mtx,method = "inverse")

heatmap_visual(distance_mtx,normalize ="both",annotation_col = categories_col, annotation_row = categories_row, annotation_colors = annoCol, annotation_legend = F, legend = TRUE, color = col.pal)


