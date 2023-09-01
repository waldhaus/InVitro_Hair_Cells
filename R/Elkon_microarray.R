library(limma)
library(illuminaio)
library(tidyverse)

x <- read.ilmn(files = "/home/rstudio/keith_data/microarray/GSE64543_non-normalized_3.txt",probeid="ID_REF",expr="signal",other.columns = "Detection")
expressed <- rowSums(x$other$Detection < 0.01) > 2
x <- x[expressed,]
y <- normalizeBetweenArrays(x,method = "quantile")

plotMDS(y,labels=c("CHC","CHC","CHC","CSC","CSC","CSC","CMES","CMES","CMES","VHC","VHC","VHC","VSC","VSC","VSC","VMES","VMES","VMES"))

probe_id <- readBGX("/home/rstudio/keith_data/microarray/GPL6887_MouseWG-6_V2_0_R0_11278593_A.bgx")
gene_id <- data.frame(x=probe_id[["probes"]][["Probe_Id"]],y=probe_id[["probes"]][["Symbol"]])
colnames(gene_id) <- c("id","gene")

# Compare between cochlear HC & vestibular
CHC_VHC <- y[,c("C.hc.rep1","C.hc.rep2","C.hc.rep3","V.hc.rep1","V.hc.rep2","V.hc.rep3")]
design <- cbind(CHC=c(1,1,1,0,0,0),VHC=c(0,0,0,1,1,1))
fit <- lmFit(CHC_VHC$E,design)
cont.matrix <- makeContrasts(CHCvsVHC=CHC-VHC, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
CHC_VHC_DEG = topTable(fit2, adjust="BH",n=Inf,sort="logFC",resort.by="logFC",p.value=0.05)

CHC_DEG = CHC_VHC_DEG %>% mutate(id=rownames(CHC_VHC_DEG)) %>% left_join(gene_id) %>% filter(logFC>0) %>% select(gene)
VHC_DEG = CHC_VHC_DEG %>% mutate(id=rownames(CHC_VHC_DEG)) %>% left_join(gene_id) %>% filter(logFC<0) %>% select(gene)

# Compare between cochlear SC & vestibular SC
CSC_VSC <- y[,c("C.sc.rep1","C.sc.rep2","C.sc.rep3","V.sc.rep1","V.sc.rep2","V.sc.rep3")]
design <- cbind(CSC=c(1,1,1,0,0,0),VSC=c(0,0,0,1,1,1))
fit <- lmFit(CSC_VSC$E,design)
cont.matrix <- makeContrasts(CSCvsVSC=CSC-VSC, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
CSC_VSC_DEG = topTable(fit2, adjust="BH",n=Inf,sort="logFC",resort.by="logFC",p.value=0.05)

CSC_DEG = CSC_VSC_DEG %>% mutate(id=rownames(CSC_VSC_DEG)) %>% left_join(gene_id) %>% filter(logFC>0) %>% select(gene)
VSC_DEG = CSC_VSC_DEG %>% mutate(id=rownames(CSC_VSC_DEG)) %>% left_join(gene_id) %>% filter(logFC<0) %>% select(gene)

# Compare between cochlear HC & cochlear SC
CHC_CSC <- y[,c("C.hc.rep1","C.hc.rep2","C.hc.rep3","C.sc.rep1","C.sc.rep2","C.sc.rep3")]
design <- cbind(CHC=c(1,1,1,0,0,0),CSC=c(0,0,0,1,1,1))
fit <- lmFit(CHC_CSC$E,design)
cont.matrix <- makeContrasts(CHCvsCSC=CHC-CSC, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
CHC_CSC_DEG = topTable(fit2, adjust="BH",n=Inf,sort="logFC",resort.by="logFC")

# Compare between vestibular HC & vestibular SC
VHC_VSC <- y[,c("V.hc.rep1","V.hc.rep2","V.hc.rep3","V.sc.rep1","V.sc.rep2","V.sc.rep3")]
design <- cbind(VHC=c(1,1,1,0,0,0),VSC=c(0,0,0,1,1,1))
fit <- lmFit(VHC_VSC$E,design)
cont.matrix <- makeContrasts(VHCvsVSC=VHC-VSC, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
VHC_VSC_DEG = topTable(fit2, adjust="BH",n=Inf,sort="logFC",resort.by="logFC")

CHC_only_DEG = intersect(CHC_DEG,CHC_CSC_DEG %>% mutate(id=rownames(CHC_CSC_DEG)) %>% left_join(gene_id) %>% filter(logFC>0) %>% select(gene))
VHC_only_DEG = intersect(VHC_DEG,VHC_VSC_DEG %>% mutate(id=rownames(VHC_VSC_DEG)) %>% left_join(gene_id) %>% filter(logFC>0) %>% select(gene))
CSC_only_DEG = intersect(CSC_DEG,CHC_CSC_DEG %>% mutate(id=rownames(CHC_CSC_DEG)) %>% left_join(gene_id) %>% filter(logFC<0) %>% select(gene))
VSC_only_DEG = intersect(VSC_DEG,VHC_VSC_DEG %>% mutate(id=rownames(VHC_VSC_DEG)) %>% left_join(gene_id) %>% filter(logFC<0) %>% select(gene))
