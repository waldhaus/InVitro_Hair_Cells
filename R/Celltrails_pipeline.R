library(CellTrails)
library(ggplotify)
library(tidyverse)

# Import in vitro data and subset HCs & SCs
invitro_haircell <- readRDS("~/aggr/Invitro_22clusters.rds")
invitro_HC_SC = subset(invitro_haircell,idents = c(3,4,5,14,15))
data <- as.matrix(invitro_HC_SC@assays$RNA@data)
data_sub = data[apply(data,1,FUN = function(x) sum(x) != 0),]
HC_SC <- SingleCellExperiment(assays=list(logcounts=data_sub))

# Select trajectory features
trajFeatureNames(HC_SC) <- featureNames(HC_SC)
trajFeatureNames(HC_SC) <- filterTrajFeaturesByFF(HC_SC, threshold=0.7)
tfeat <- featureNames(HC_SC[trajFeatureNames(HC_SC), ])
trajFeatureNames(HC_SC) <- tfeat
showTrajInfo(HC_SC)

# Perform spectral embedding
se <- embedSamples(HC_SC)
d <- findSpectrum(se$eigenvalues, frac=100)
latentSpace(HC_SC) <- se$components[, d]

# Perform clustering
cl <- findStates(HC_SC, min_size=0.02, min_feat=5, max_pval=1e-4, min_fc=1.8)
states(HC_SC) <- cl

# Sample ordering
HC_SC <- connectStates(HC_SC, l=10)
HC_SC <- fitTrajectory(HC_SC)

# Export trajectory graph structure to graphml and then import layout computed from yEd
write.ygraphml(HC_SC, file="~/celltrails/HC_SC_8states.graphml", color_by="phenoName",name="state", node_label="state")
tl <- read.ygraphml(file="~/celltrails/HC_SC_8states_layout.graphml")

# Adjust layout and store to object
trajLayout(HC_SC, adjust=TRUE) <- tl

plotMap(HC_SC, color_by="phenoName", name="state")

saveRDS(HC_SC, "/home/rstudio/keith_data/Invitro_HC_SC_sce.rds")

# Remove cells from S8 and reclustering
Idents(invitro_HC_SC) = as.vector(HC_SC@colData@listData$CellTrails.state)
invitro_HC_SC_rmS8 = subset(invitro_HC_SC,idents = c('S8'), invert = TRUE)
data <- as.matrix(invitro_HC_SC_rmS8@assays$RNA@data)
HC_SC_sub <- SingleCellExperiment(assays=list(logcounts=data[apply(data,1,FUN = function(x) sum(x) != 0),]))

trajFeatureNames(HC_SC_sub) <- featureNames(HC_SC_sub)
trajFeatureNames(HC_SC_sub) <- filterTrajFeaturesByFF(HC_SC_sub, threshold=0.7)
tfeat <- featureNames(HC_SC_sub[trajFeatureNames(HC_SC_sub), ])
trajFeatureNames(HC_SC_sub) <- tfeat
se <- embedSamples(HC_SC_sub)
d <- findSpectrum(se$eigenvalues, frac=100)
latentSpace(HC_SC_sub) <- se$components[, d]
cl <- findStates(HC_SC_sub, min_size=0.01, min_feat=5, max_pval=1e-4, min_fc=1.7)
states(HC_SC_sub) <- cl
HC_SC_sub <- connectStates(HC_SC_sub, l=10)
HC_SC_sub <- fitTrajectory(HC_SC_sub)
write.ygraphml(HC_SC_sub, file="~/HC_SC_rmS8_11states.graphml", color_by="phenoName",name="state", node_label="state")
tl <- read.ygraphml(file="~/HC_SC_rmS8_11states_layout.graphml")
trajLayout(HC_SC_sub, adjust=TRUE) <- tl

saveRDS(HC_SC_sub, "~/Invitro_HC_SC_rmS8_11states_sce.rds")



