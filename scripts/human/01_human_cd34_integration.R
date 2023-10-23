################################################
#
# 01_human_cd34_integration
# Author: Maximilian Schönung
# Date: 23.10.2023
#
################################################

# Set the paths -----------------------------------------------------------

setwd("/omics/odcf/analysis/OE0565_projects/ptpn11/")
plot.dir <- "/omics/odcf/analysis/OE0565_projects/ptpn11/plots/"
table.dir <- "/omics/odcf/analysis/OE0565_projects/ptpn11/tables/"

subDir <- "04_human/"
if (!file.exists(paste0(plot.dir,subDir))){
  dir.create(paste0(plot.dir,subDir),recursive = T)
}

# Load the libraries ------------------------------------------------------
library(Seurat)
library(ggplot2)
library(RColorBrewer)

# Load the 10x Data -------------------------------------------------------
ref <- readRDS("/omics/odcf/analysis/OE0565_projects/ptpn11/data_human/JMML_MNC-Postnatal_REF_Sampled-Anchor.Seurat.PC1-18.rds")
cd34 <- readRDS("/omics/odcf/analysis/OE0219_projects/JMMLC/scRNA_Data/2.Single_Dataset/JMML_CD34/29/JMML_CD34_PC_29.SeuratV3.rds")
anno <- readRDS("/omics/odcf/analysis/OE0565_projects/ptpn11/data_human/10X_meta.data_complete.rds")

# Find the anchors --------------------------------------------------------
anchors <- FindTransferAnchors(reference = ref, query = cd34, reference.reduction = "pca")
ref2 <- RunUMAP(ref,dims = 1:18,return.model = T)
cd34.query <- MapQuery(anchorset = anchors, reference = ref2, query = cd34,
                           refdata = list(celltype = "predicted.id"), reference.reduction = "pca", reduction.model = "umap")

p1 <- DimPlot(ref2, reduction = "umap", group.by = "predicted.id", label = TRUE, label.size = 3,
              repel = TRUE) + NoLegend() + ggtitle("Reference annotations")
p2 <- DimPlot(cd34.query, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE,
              label.size = 3, repel = TRUE) + NoLegend() + ggtitle("Query transferred labels")

p1+p2

ref.df <- as.data.frame(ref2@reductions$umap@cell.embeddings)
query.df <- as.data.frame(cd34.query@reductions$ref.umap@cell.embeddings)
query.meta <- cd34.query@meta.data
all(rownames(cd34@meta.data)==rownames(query.meta))
query.meta$old_id <- anno[rownames(query.meta),"Celltype"]
query.meta$old_id[query.meta$old_id!="CD34+ HSC"] <- "CD34"
table(query.meta$old_id)

pdf(paste0(plot.dir,subDir,Sys.Date(),"_cd34_integration_umap.pdf"),height = 7,width =7)
ggplot(ref.df,aes(UMAP_1,UMAP_2))+
  geom_point(colour="#4f2973",alpha=.01,size=.1)+
  theme_classic()+
  geom_point(data = query.df,aes(refUMAP_1,refUMAP_2,colour=query.meta$old_id),alpha=.5,size=.1)+
  scale_color_manual(values=c("CD34+ HSC"="#252525","CD34"="#b3b3b3"))+
  theme(legend.position="none")
ggplot(ref.df,aes(UMAP_1,UMAP_2))+
  geom_point(colour="#4f2973",alpha=.05,size=.1)+
  theme_classic()+
  geom_point(data = query.df,aes(refUMAP_1,refUMAP_2,colour=query.meta$old_id),alpha=.5,size=.1)+
  scale_color_manual(values=c("CD34+ HSC"="#252525","CD34"="#b3b3b3"))+
  theme(legend.position="none")
dev.off()