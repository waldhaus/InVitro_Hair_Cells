# InVitro_Hair_Cells

Abstract
============================
Inner ear organoids recapitulate development and are intended to generate cell types of the otic lineage for applications such as basic science research and cell replacement strategies. Here, we use single-cell sequencing to study the cellular heterogeneity of late-stage mouse inner ear organoid sensory epithelia, which we validated by comparison with data sets of the mouse cochlea and vestibular epithelia. We resolved supporting cell sub-types, cochlear like hair cells, and vestibular Type I and Type II like hair cells.  While cochlear like hair cells aligned best with an outer hair cell trajectory, vestibular like hair cells followed developmental trajectories similar to in vivo programs branching into Type II and then Type I extrastriolar hair cells. These results highlight the transcriptional accuracy of the organoid developmental program but will also inform future strategies to improve synaptic connectivity and regional specification.

Analysis code
============================
Here we include analysis code for our manuscript.
* Seurat pipeline: Seurat pipeline is to analyze three scRNA-seq datasets derived from in vitro generated mouse inner ear organoids. Datasets of 16, 20 and 21 days-in-vitro were combined and analyzed together. We identified distinct clusters, annotated cell types for each cluster based on marker genes, and determined differentially expressed genes across various cell types.
* scRNA-seq integration pipeline: We used Seurat to perform joint analysis of in vitro data and two in vivo scRNA-seq datasets. We co-embedded these three datasets onto the same UMAP. We followed the documentation from https://satijalab.org/seurat/articles/integration_introduction.html.
* Celltrails pipeline: CellTrails pipeline is to determine sub-populations of in vitro generated cells. Here we focused on HCLCs and SCLCs. Following the initial clustering into 8 states, we determined the lineage identity of each state by calculating cell type specific enrichment score with SingleCellSignatureExplorer. Then, we removed cells from S8 which correspond to cochlear HCs and reran clustering. This resulted in a total of 11 states. We followed the documentation from https://hellerlab.stanford.edu/celltrails/.
Single-Cell Signature Explorer is a package of four packages to compute and visualize gene set enrichment score. We used Single-Cell Signature Scorer to compute gene set enrichment score in each single cell and then Single-Cell Signature Merger to collate the signature scores table with UMAP coordinates.Single-Cell Signature Viewer enabled us to display enrichment score on a UMAP.
* Elkon_microarray.R: We performed differential expression analysis on auditory and vestibular HCs and SCs from a microarray dataset
using limma. We followed limma pipeline from https://bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf
* find_dynamic_gene.R: `find_dynamic_gene.R` script includes the function to identify genes dynamically expressed along a trajectory trail from scRNA-seq data.
* Similarity_heatmap.R: `Similarity_heatmap.R` is to generate similarity heatmap in Figure 6I, 6J that calculates pairwise cell similarity score between in vivo and in vitro trails. The first step is to find top 2000 variable features within in vivo and in vitro cells from the corresponding trail. Next, we aggregated 6 neighbor cells into metacells in order to denoise the dataset for both in vivo and in vitro data. By leveraging highly variable genes, we calcualted the distance for each pair of metacells between in vivo and in vitro and convert distance to similarity, followed by visualization as heatmap.
* Pseudotime_heatmap.R: `Pseudotime_heatmap.R` includes function to fit the average gene expression as a function of pseudotime using generalized additive model within a trail. The smoothed expression profiles are visualized in heatmaps and genes are ordered based on the similarity of their expression patterns.
* Utilities.R: `Utilities.R` script contains utility functions.
* scVelo.ipynb: It contains the code how we use scVelo to perform RNA velocity analysis in in vitro scRNA-seq data. It enables us to infer transcriptional dynamics that are not directly observable in a scRNA-seq experiment, accomplished through a mathematical model of transcriptional kinetics. We computed RNA velovity and visualized them in a 2D space.
