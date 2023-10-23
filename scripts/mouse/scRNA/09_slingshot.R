################################################
#
# 09_singshot
# Author: Maximilian Schönung
# Date: 23.10.2023
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
library(slingshot)
library(tradeSeq)

# Cluster colour ----------------------------------------------------------
col.vec = c("HSC"="black","MPP1"="grey89","MPP2"="grey50",
            "MPP_Cycling1"="ivory3","MPP_Cycling2"="ivory4",
            "MEP"="coral2","Ery1"="orangered2","Ery2"="orangered3","MPP_EryMeg"="darkred",
            "MyP1"="darkolivegreen4","MyP2"="darkolivegreen3","MonoP"="darkolivegreen2",
            "NeutroP"="palegreen3",
            "Mono1"="goldenrod","Mono2"="gold3","DC"="gold",
            "Neutro1"="orange","Neutro2"="darkorange","Neutro3"="darkorange1","Neutro4"="darkorange2","Neutro5"="darkorange3",
            "MPP_Lympho"="lightsteelblue4","LyP"="steelblue1","T-cell"="lightskyblue",
            "B-cell1"="paleturquoise","B-cell2"="paleturquoise1","B-cell3"="paleturquoise2","B-cell4"="paleturquoise3",
            "Basophils"="blueviolet")

# Load the 10x Data -------------------------------------------------------
seurat <- readRDS(paste0(data.dir,"2023-09-18_bm_ptpn11_sct_annotated.RDS"))
seurat.sub <- subset(seurat,idents = c("HSC","MPP1","MPP2","MyP1","MyP2","MonoP","Mono1","Mono2"))
rm(seurat)
sce <- as.SingleCellExperiment(seurat.sub)
reducedDim(sce, "UMAP") <- seurat.sub@reductions$umap@cell.embeddings

# Extract the data for Slingshot ------------------------------------------
set.seed(1234)
pseudotime <- slingshot(sce,reducedDim="UMAP", start.clus="HSC")
pseudotime$PseudotimeRank <- rank(-pseudotime$slingPseudotime_1)
curves <- slingCurves(pseudotime, as.df = TRUE)
sling.weights <- slingCurveWeights(pseudotime)
umap.sub <- as.data.frame(reducedDim(pseudotime,"UMAP"))


pdf(paste0(plot.dir,subDir,Sys.Date(),"_slingshot_myeloid.pdf"))
ggplot(umap.sub, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(colour=pseudotime$ident)) +
  scale_colour_manual(values=col.vec) +
  theme_classic()+
  geom_path(data = curves %>% arrange(Order),
            aes(group = Lineage),colour="red")
ggplot(umap.sub,aes(UMAP_1,UMAP_2,colour=pseudotime$slingPseudotime_1))+
  geom_point(size=.1)+
  theme_classic()+
  viridis::scale_colour_viridis(option="A")
ggplot(umap.sub,aes(UMAP_1,UMAP_2,colour=pseudotime$PseudotimeRank))+
  geom_point(size=.1)+
  theme_classic()+
  viridis::scale_colour_viridis(option="A")
dev.off()

# Plot the lineage progression --------------------------------------------
umap.sub$PseudotimeRank <- pseudotime$PseudotimeRank
umap.sub$Pseudotime <- -pseudotime$slingPseudotime_1
umap.sub$leuk <- factor(pseudotime$leuk,labels = c("WT","E76K"))
umap.sub$ident <- factor(pseudotime$ident)

# downsample
set.seed(1234)
ind1 <- sample(which(umap.sub$leuk=="E76K"),5000)
ind2 <- sample(which(umap.sub$leuk=="WT"),5000)
ind.all <- c(ind1,ind2)
umap.down <- umap.sub[ind.all,]

pdf(paste0(plot.dir,subDir,Sys.Date(),"_diffaxis_myeloid_slingshot.pdf"),height = 4)
ggplot(umap.down,aes(PseudotimeRank,fill=ident))+
  geom_density(aes(x=PseudotimeRank, y=(..count../sum(..count..))*100, fill=ident), alpha=1/2)+
  theme_classic()+
  scale_fill_manual(values = col.vec[names(col.vec)%in%levels(umap.sub$ident)])+
  ylab("Percentage of Cells")+
  xlab("Differentiation Pseudotime Rank (HSPC - Monocytes)")
ggplot(umap.down,aes(Pseudotime,fill=ident))+
  geom_density(aes(x=Pseudotime, y=(..count../sum(..count..))*100, fill=ident), alpha=1/2)+
  theme_classic()+
  scale_fill_manual(values = col.vec[names(col.vec)%in%levels(umap.sub$ident)])+
  ylab("Percentage of Cells")+
  xlab("Differentiation Pseudotime (HSPC - Monocytes)")
ggplot(umap.down)+
  geom_density(data=umap.down[umap.down$leuk=="E76K",],aes(x=PseudotimeRank, y=(..count../sum(..count..))*100, fill=leuk), alpha=1/2)+
  geom_density(data=umap.down[umap.down$leuk=="WT",],aes(x=PseudotimeRank, y=(..count../sum(..count..))*100, fill=leuk), alpha=1/2)+
  theme_classic()+
  scale_fill_manual(values=c("red3","royalblue4"))+
  geom_rug(aes(x=PseudotimeRank,colour=ident))+
  scale_colour_manual(values = col.vec[names(col.vec)%in%levels(umap.sub$ident)])+
  ylab("Percentage of Cells per Genotype")+
  xlab("Differentiation Pseudotime Rank (HSPC - Monocytes)")
ggplot(umap.down)+
  geom_density(data=umap.down[umap.down$leuk=="E76K",],aes(x=Pseudotime, y=(..count../sum(..count..))*100, fill=leuk), alpha=1/2)+
  geom_density(data=umap.down[umap.down$leuk=="WT",],aes(x=Pseudotime, y=(..count../sum(..count..))*100, fill=leuk), alpha=1/2)+
  theme_classic()+
  scale_fill_manual(values=c("red3","royalblue4"))+
  geom_rug(aes(x=Pseudotime,colour=ident))+
  scale_colour_manual(values = col.vec[names(col.vec)%in%levels(umap.sub$ident)])+
  ylab("Percentage of Cells per Genotype")+
  xlab("Differentiation Pseudotime (HSPC - Monocytes)")
dev.off()

