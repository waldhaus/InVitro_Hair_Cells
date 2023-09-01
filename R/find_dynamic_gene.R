#' find_dynamic_gene() is to find genes dynamically expressed along a trajectory trail
#' @param obj A SingleCellExperiment object
#' @param gene_ls A vector contains all genes to be test
#' @param trail A string of trail name
#' @param threshold corrected p-value threshold for significance. By default, `threshold` is 0.01.
#' @return The function returns 1 vector of dynamically expressed genes
find_dynamic_gene <- function(obj, gene_ls, trail, threshold = 0.01){
  pvalue = NULL
  for (gene in gene_ls){
    fit <- fitDynamic(obj,gene,trail)
    if (is.na(fit$gam)){
      pvalue = c(pvalue,NA)
    } else{
    pvalue = c(pvalue,summary(fit$gam)$s.p)
    }
  }
  result = data.frame(x=gene_ls,y=pvalue)
  result = result[!is.na(result$y),]
  result$fdr <- p.adjust(result$y, method = "BH", n = length(result$y))
  dynamic_gene = result[result$fdr<threshold,]$x %>% sort()
  
  return(dynamic_gene)
}

library(CellTrails)
library(tidyverse)

Jan_sce <- readRDS("~/external_data/Jan_2021_sce_wt.mouse.utricle_805_cells.rds")

Jan_sce <- connectStates(Jan_sce, l=7)
Jan_sce <- fitTrajectory(Jan_sce)

tl <- read.ygraphml(file="~/celltrails/Jan_12states_layout.graphml")
trajLayout(Jan_sce, adjust=TRUE) <- tl

Jan_sce <- addTrail(Jan_sce, from="U1", to="H3", name="trailA")
Jan_sce <- addTrail(Jan_sce, from="U1", to="H5", name="trailB")

data <- as.matrix(Jan_sce@assays@data@listData[["logcounts"]])
genes = rownames(data[apply(data,1,FUN = function(x) sum(x) != 0),])

trailA_gene <- find_dynamic_gene(Jan_sce,genes,'trailA')
write.table(data.frame(x=trailA_gene),"~/Jan_trailA_dynamic_genes.txt",row.names = FALSE)

trailB_gene <- find_dynamic_gene(Jan_sce,genes,'trailB')
write.table(data.frame(x=trailB_gene),"~/Jan_trailB_dynamic_genes.txt",row.names = FALSE)

