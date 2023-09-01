library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

# Import in vitro data and subset HCs & SCs
invitro_haircell <- readRDS("~/aggr/Invitro_22clusters.rds")
invitro_HC_SC = subset(invitro_haircell,idents = c(3,4,5,14,15))

# Import Jan data
Jan_sce <- readRDS("~/external_data/Jan_2021_sce_wt.mouse.utricle.rds")
Jan_object <- as.Seurat(Jan_sce, counts = "counts", data = "logcounts")
Jan_object <- FindVariableFeatures(Jan_object, selection.method = "vst", nfeatures = 2000)

# Import P2 data and subset relevant cell types
P2_object <- readRDS("~/P2_scRNAseq_seurat_object.RDS")
P2.updated = UpdateSeuratObject(object = P2_object)
P2.relevant = subset(P2.updated,idents = c("HC","medSC","PC/DC/IBC","latSC","RF"))
P2.relevant <- FindVariableFeatures(P2.relevant, selection.method = "vst", nfeatures = 2000)

all_list <- list(P2.relevant,Jan_object,invitro_HC_SC)

# Perform integration
all_datasets.anchors <- FindIntegrationAnchors(object.list = all_list, dims = 1:25, anchor.features = 3500, normalization.method = "LogNormalize")
all_datasets.combined <- IntegrateData(anchorset = all_datasets.anchors, dims = 1:25, normalization.method = "LogNormalize")

# Run the standard workflow for visualization
DefaultAssay(all_datasets.combined) <- "integrated"

all_genes = rownames(all_datasets.combined)
all_datasets.combined <- ScaleData(all_datasets.combined, features = all_genes)
all_datasets.combined <- RunPCA(all_datasets.combined, npcs = 30, verbose = FALSE)
all_datasets.combined <- RunUMAP(all_datasets.combined, dims = 1:12, seed = 20)

