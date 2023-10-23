################################################
#
# 05_dim_qc
# Author: Maximilian Schönung
# Date: 22.06.2023
#
################################################


# Set the paths -----------------------------------------------------------
setwd("/omics/odcf/analysis/OE0565_projects/ptpn11/")
data.dir <- "/omics/odcf/analysis/OE0565_projects/ptpn11/data/"
plot.dir <- "/omics/odcf/analysis/OE0565_projects/ptpn11/plots/"


subDir <- "01_DimRed/"
if (!file.exists(paste0(plot.dir,subDir))){
  dir.create(paste0(plot.dir,subDir),recursive = T)
}

# Load the libraries ------------------------------------------------------
library(Seurat)
library(ggplot2)
library(dplyr)
library(Nebulosa)

# Load the 10x Data & Run PCA -------------------------------------------------------
seurat <- readRDS(paste0(data.dir,"2023-09-18_bm_ptpn11_sct_annotated.RDS"))
seurat@meta.data$celltype <- Idents(seurat)
meta <- seurat@meta.data
meta$celltype <- Idents(seurat)
meta$UMAP1 <- seurat@reductions$umap@cell.embeddings[,1]
meta$UMAP2 <- seurat@reductions$umap@cell.embeddings[,2]

# Again run the QC for generating plots after filtering -------------------
meta$orig.ident <- factor(meta$orig.ident,levels = c("345BM","347BM","349BM","348BM","351BM","352BM"))
meta$leuk <- factor(meta$leuk,labels=c("WT","E76K"))
meta.sub <- meta[,c("orig.ident","leuk","id","nCount_RNA","nFeature_RNA","percent.mt")]
meta.melt <- reshape2::melt(meta.sub,id.vars=c("orig.ident","leuk","id"))

pdf(paste0(plot.dir,subDir,Sys.Date(),"_cluster_qc_filtered.pdf"),width = 10,height=5)
ggplot(meta.melt,aes(orig.ident,value,fill=leuk))+
  geom_violin()+
  theme_classic()+
  facet_wrap(~variable,scales="free")+
  xlab("Identity")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_fill_manual(values=c("royalblue4","red3"))
dev.off()

# Plot the UMAP with different QC metrices --------------------------------
set.seed(1234)
rand <- sample(1:nrow(meta),10000)
pdf(paste0(plot.dir,subDir,Sys.Date(),"_cluster_qc.pdf"))
ggplot(meta,aes(UMAP1,UMAP2,colour=factor(leuk)))+
  geom_point(size=.1,alpha=.2)+
  theme_classic()+
  scale_colour_manual(values=c("royalblue4","red3"))
ggplot(meta[rand,],aes(UMAP1,UMAP2,colour=factor(leuk)))+
  geom_point(size=.1)+
  theme_classic()+
  scale_colour_manual(values=c("royalblue4","red3"))
ggplot(meta,aes(UMAP1,UMAP2,colour=factor(orig.ident)))+
  geom_point(size=.1,alpha=.5)+
  theme_classic()
ggplot(meta,aes(UMAP1,UMAP2,colour=factor(Phase)))+
  geom_point(size=.1,alpha=.5)+
  theme_classic()
ggplot(meta,aes(UMAP1,UMAP2,colour=G2M.Score))+
  geom_point(size=.1,alpha=.2)+
  scale_colour_continuous(type = "viridis")+
  theme_classic()
ggplot(meta,aes(UMAP1,UMAP2,colour=S.Score))+
  geom_point(size=.1,alpha=.2)+
  scale_colour_continuous(type = "viridis")+
  theme_classic()
ggplot(meta,aes(UMAP1,UMAP2,colour=percent.mt))+
  geom_point(size=.1,alpha=.2)+
  scale_colour_continuous(type = "viridis")+
  theme_classic()
ggplot(meta,aes(UMAP1,UMAP2,colour=nCount_RNA))+
  geom_point(size=.1,alpha=.2)+
  scale_colour_continuous(type = "viridis")+
  theme_classic()
ggplot(meta,aes(UMAP1,UMAP2,colour=nFeature_RNA))+
  geom_point(size=.1,alpha=.2)+
  scale_colour_continuous(type = "viridis")+
  theme_classic()
dev.off()

# Percentage of cells from each Genotype in each cluster ------------------
pdf(paste0(plot.dir,subDir,Sys.Date(),"_vln_qc.pdf"),width = 12)
VlnPlot(seurat, features = "nFeature_RNA", group.by = "celltype", pt.size = 0, combine = FALSE)
VlnPlot(seurat, features = "nFeature_RNA", split.by = "leuk", group.by = "celltype", pt.size = 0, combine = FALSE)
VlnPlot(seurat, features = "nCount_RNA", group.by = "celltype", pt.size = 0, combine = FALSE)
VlnPlot(seurat, features = "nCount_RNA", split.by = "leuk", group.by = "celltype", pt.size = 0, combine = FALSE)
VlnPlot(seurat, features = "percent.mt", group.by = "celltype", pt.size = 0, combine = FALSE)
VlnPlot(seurat, features = "percent.mt", split.by = "leuk", group.by = "celltype", pt.size = 0, combine = FALSE)
dev.off()

plot.df <- as.data.frame(table(seurat$celltype,seurat$Phase,seurat$leuk))

pdf(paste0(plot.dir,subDir,Sys.Date(),"_bar_qc.pdf"))
ggplot(plot.df,aes(Var3,Freq,fill=Var2))+
  geom_bar(position = "fill",stat = "identity") +
  scale_y_continuous(labels = scales::percent_format())+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  facet_wrap(~Var1)
dev.off()


# Marker per Cluster ------------------------------------------------------
png(paste0(plot.dir,subDir,Sys.Date(),"_density_marker.png"),width=2400, height=1600)
plot_density(seurat,features=c("Meis1","Mecom","Cd34","Pcna","Flt3","Cd48","Slamf1","Fcgr3","Mpo","Csf1r","Gfi1b","Kit"),size=.1,pal="magma")& NoAxes()& NoLegend()
dev.off()

png(paste0(plot.dir,subDir,Sys.Date(),"_density_marker2.png"),width=2400, height=1600)
plot_density(seurat,features=c("Elane","Fcnb","F13a1","Siglech","Cd7","Itgb7","Naaa","Cebpe","Irf8","Ly6g","Cd3e","Ebf1"),size=.1,pal="magma")& NoAxes()& NoLegend()
dev.off()


pdf(paste0(plot.dir,subDir,Sys.Date(),"_dotplot_HSC_MPP.pdf"),width = 10)
DotPlot(seurat,features=c("Procr","Meis1","Mecom","Cd34","Kit","Cd48","Flt3","Gfi1b"),
        idents=c("HSC","MPP1","MPP2","MPP_Cycling1","MPP_Cycling2","MPP_EryMeg","MPP_Lympho"),
        cluster.idents=T,dot.scale = 9)
dev.off()
