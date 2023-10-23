################################################
#
# 08_gene_sets.R
# Author: Maximilian Schönung
# Date: 26.09.2023
#
################################################

# Set the paths -----------------------------------------------------------
setwd("/omics/odcf/analysis/OE0565_projects/ptpn11/")
data.dir <- "/omics/odcf/analysis/OE0565_projects/ptpn11/data/"
plot.dir <- "/omics/odcf/analysis/OE0565_projects/ptpn11/plots/"

subDir <- "02_functional/"
if (!file.exists(paste0(plot.dir,subDir))){
  dir.create(paste0(plot.dir,subDir),recursive = T)
}

# Load the libraries ------------------------------------------------------
library(Seurat)
library(dplyr)
library(ggplot2)

# Load the 10x Data -------------------------------------------------------
seurat <- readRDS(paste0(data.dir,"2023-09-18_bm_ptpn11_sct_annotated.RDS"))

# Downsample the data ----------------------------------------------------
meta <- seurat@meta.data
meta$anno <- Idents(seurat)
meta$index <- 1:nrow(meta)
meta <- meta[meta$anno=="HSC",]
ind1 <- meta$index[sample(which(meta$leuk==0),1000)]
ind2 <- meta$index[sample(which(meta$leuk==1),1000)]
ind.all <- c(ind1,ind2)
seurat.hsc <- seurat[,ind.all]
rm(seurat)

# Load the gene_set --------------------------------------------------------------
geneset <- read.delim("data/gene_sets/izzo_2020_ng_gene_modules.txt",header = T,stringsAsFactors = F)
mono_set <- geneset[geneset$Module==16,"Gene"]
ery_set <- geneset[geneset$Module==11,"Gene"]
lym_set <- geneset[geneset$Module==23,"Gene"]
stem_set <- geneset[geneset$Module==4,"Gene"]

# Add the gene_set --------------------------------------------------------
seurat.hsc <- AddModuleScore(seurat.hsc,features=list(mono_set),name="Mono_Score")
seurat.hsc <- AddModuleScore(seurat.hsc,features=list(ery_set),name="Ery_Score")
seurat.hsc <- AddModuleScore(seurat.hsc,features=list(lym_set),name="Lym_Score")
seurat.hsc <- AddModuleScore(seurat.hsc,features=list(stem_set),name="Stem_Score")
seurat.hsc@meta.data$leuk <- factor(seurat.hsc@meta.data$leuk,labels = c("WT","E76K"))

# Plot the score ----------------------------------------------------------
pdf(paste0(plot.dir,subDir,Sys.Date(),"_gene_sets_seurat_score.pdf"))
ggplot(seurat.hsc@meta.data,aes(Mono_Score1,Stem_Score1,fill=leuk))+
  stat_density_2d(geom = "polygon", aes(alpha = ..level.., fill = leuk))+
  facet_wrap(~leuk)+
  scale_fill_manual(values=c("gray40","firebrick"))+
  theme_classic()+
  geom_vline(xintercept=median(seurat.hsc@meta.data[seurat.hsc@meta.data$leuk=="WT","Mono_Score1"]))+
  geom_hline(yintercept=mean(seurat.hsc@meta.data[seurat.hsc@meta.data$leuk=="WT","Stem_Score1"]))
ggplot(seurat.hsc@meta.data,aes(Ery_Score1,Stem_Score1,fill=leuk))+
  stat_density_2d(geom = "polygon", aes(alpha = ..level.., fill = leuk))+
  facet_wrap(~leuk)+
  scale_fill_manual(values=c("gray40","firebrick"))+
  theme_classic()+
  geom_vline(xintercept=median(seurat.hsc@meta.data[seurat.hsc@meta.data$leuk=="WT","Ery_Score1"]))+
  geom_hline(yintercept=mean(seurat.hsc@meta.data[seurat.hsc@meta.data$leuk=="WT","Stem_Score1"]))
ggplot(seurat.hsc@meta.data,aes(Lym_Score1,Stem_Score1,fill=leuk))+
  stat_density_2d(geom = "polygon", aes(alpha = ..level.., fill = leuk))+
  facet_wrap(~leuk)+
  scale_fill_manual(values=c("gray40","firebrick"))+
  theme_classic()+
  geom_vline(xintercept=median(seurat.hsc@meta.data[seurat.hsc@meta.data$leuk=="WT","Lym_Score1"]))+
  geom_hline(yintercept=mean(seurat.hsc@meta.data[seurat.hsc@meta.data$leuk=="WT","Stem_Score1"]))
dev.off()

#Test for Mono signature
LMM1.stem <- lme4::lmer(Mono_Score1~leuk+(1|orig.ident),data = seurat.hsc@meta.data)
LMM2.stem <- lme4::lmer(Mono_Score1~(1|orig.ident),data = seurat.hsc@meta.data)
anova(LMM2.stem,LMM1.stem)

#Test for stem signature
LMM1.stem <- lme4::lmer(Stem_Score1~leuk+(1|orig.ident),data = seurat.hsc@meta.data)
LMM2.stem <- lme4::lmer(Stem_Score1~(1|orig.ident),data = seurat.hsc@meta.data)
anova(LMM2.stem,LMM1.stem)
