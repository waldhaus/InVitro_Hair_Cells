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
* Utilities.R: Utilities.R script contains utility functions.
* scVelo.ipynb: It contains the code how we use scVelo to perform RNA velocity analysis in in vitro scRNA-seq data. It enables us to infer transcriptional dynamics that are not directly observable in a scRNA-seq experiment, accomplished through a mathematical model of transcriptional kinetics. We computed RNA velovity and visualized them in a 2D space.
