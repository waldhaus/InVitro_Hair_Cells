## Seurat scRNA-seq cluster analysis pipeline
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(cowplot)
library(ggpubr)

## load scRNA-seq data
invitro_haircell.data <- Read10X(data.dir = "~/aggr/filtered_feature_bc_matrix/")
invitro_haircell <- CreateSeuratObject(counts = invitro_haircell.data, project = "invitro_scRNAseq_seurat_object")
invitro_haircell[["percent.mt"]] <- PercentageFeatureSet(invitro_haircell, pattern = "^[Mm][Tt]-")

# filtering cells
invitro_haircell <- subset(invitro_haircell, subset = nFeature_RNA > 600 & nFeature_RNA < 8000 & percent.mt < 10)

# Normalization
invitro_haircell <- NormalizeData(invitro_haircell, normalization.method = "LogNormalize", scale.factor = 10000)

# Identification of highly variable features (feature selection)
invitro_haircell <- FindVariableFeatures(invitro_haircell, selection.method = "vst", nfeatures = 3500)

# Scaling the data
all_genes <- rownames(invitro_haircell)
invitro_haircell <- ScaleData(invitro_haircell, features = all_genes)

# Perform linear dimensional reduction
invitro_haircell <- RunPCA(invitro_haircell, features = VariableFeatures(object = invitro_haircell),verbose = F)

# Clustering
invitro_haircell <- FindNeighbors(object = invitro_haircell, dims = 1:15, k.param = 25)
invitro_haircell <- FindClusters(invitro_haircell, resolution = 0.6, algorithm = 4, random.seed = 10)

invitro_haircell <- RunUMAP(invitro_haircell, dims = 1:15, seed.use = 20)

color = c("#B14231","#E6242E","#F28429","#D9B03B","#F1E911","#F8E0A7","#97C83E","#5B683F","#978837","#65C0A5","#43BBEC","#3952A1","#863F93","#C17EA5","#E61B60","#FF969E","#965125","#751428","#A48F7E","#A4A3A3","#8B98AA","black")
DimPlot(invitro_haircell, reduction = "umap", label = FALSE, pt.size = 0.68, cols = color) + guides(color = guide_legend(override.aes = list(size=3),label.hjust = 0.7))

saveRDS(invitro_haircell,"~/aggr/Invitro_22clusters.rds")

# DE analysis across cell types
invitro_haircell@active.ident = factor(invitro_haircell@active.ident,levels = c("3","4","14","5","15","7","8","10","13","16","21","2","19","11","6","1","9","17","22","18","12","20"))
invitro_haircell.markers <- FindAllMarkers(invitro_haircell, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25,test.use = "wilcox")
top100_degs <- invitro_haircell.markers[invitro_haircell.markers$p_val_adj<0.05,] %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC)
