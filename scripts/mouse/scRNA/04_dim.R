################################################
#
# 04_dim
# Author: Maximilian Schönung
# Date: 22.06.2023
#
################################################


# Set the paths -----------------------------------------------------------
setwd("/omics/odcf/analysis/OE0565_projects/ptpn11/")
data.dir <- "/omics/odcf/analysis/OE0565_projects/ptpn11/data/"
plot.dir <- "/omics/odcf/analysis/OE0565_projects/ptpn11/plots/"

# Load the libraries ------------------------------------------------------
library(Seurat)
library(ggplot2)
library(dplyr)

# Load the 10x Data & Run PCA -------------------------------------------------------
seurat <- readRDS(paste0(data.dir,"2023-06-22_bm_ptpn11_doubletFilt_sctransform_integrated.RDS"))
seurat <- RunPCA(seurat,npcs = 200) # run PCA

# Dimensional Reduction ---------------------------------------------------
seurat <- RunUMAP(seurat,n.neighbors = 50,dims=1:25,min.dist=0.2)

subDir <- "01_DimRed/"
if (!file.exists(paste0(plot.dir,subDir))){
  dir.create(paste0(plot.dir,subDir),recursive = T)
}

pdf(paste0(plot.dir,subDir,Sys.Date(),"_umap_tsne.pdf"))
print(DimPlot(seurat, reduction = "umap",group.by="leuk"))
print(DimPlot(seurat, reduction = "umap",group.by="Phase"))
dev.off()


#Cluster the data ------------------------------------------------------
seurat <- FindNeighbors(seurat, dims = 1:25)
seurat <- FindClusters(seurat, resolution = .8)
table(seurat@meta.data$seurat_clusters)

pdf(paste0(plot.dir,subDir,Sys.Date(),"_tsne_umap_clust.pdf"))
print(DimPlot(seurat, reduction = "umap",group.by="seurat_clusters",raster=F))
dev.off()


umap <- data.frame(seurat@reductions$umap@cell.embeddings,"clusts"=seurat@meta.data$seurat_clusters)
summary(umap)

pdf(paste0(plot.dir,subDir,Sys.Date(),"_umap_perCluster.pdf"))
for(i in names(table(umap$clusts))){
  print(ggplot(umap,aes(UMAP_1,UMAP_2))+
          geom_point(size=0.1,alpha=0.3)+
          geom_point(data=umap[umap$clusts==i,],aes(UMAP_1,UMAP_2),colour="red")+
          theme_classic()+ theme(legend.title = element_blank())+
          xlim(-11, 18)+  ylim(-15, 16)+
          ggtitle(paste0("Cluster: ",i)))}
dev.off()


# Normalize Data together to find marker genes ----------------------------
DefaultAssay(seurat) <- "RNA"
seurat <- NormalizeData(seurat)
seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(seurat)
seurat <- ScaleData(seurat, features = all.genes)


# Find Marker ---------------------------------------------------------
markers <- FindAllMarkers(seurat, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.5)
markers.sel <- markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
write.table(markers.sel,paste0(data.dir,Sys.Date(),"_bm_top20_marker_res1.0.txt"),sep="\t",quote=F)
write.table(markers,paste0(data.dir,Sys.Date(),"_bm_all_marker_res1.0.txt"),sep="\t",quote=F)

# Feature Plot -------------------------------------------------------------
pdf(paste0(plot.dir,subDir,Sys.Date(),"_umap_marker.pdf"),height = 20,width = 20)
print(FeaturePlot(seurat, features = c("Kit","Flt3","Il7r")))
print(FeaturePlot(seurat, features = c("Stmn1", "Cd34", "Flt3", "Elane", "Mki67", "Meis1", "Ifitm1", "Il7r", 
                                       "Klf1","Vwf","Pf4","Gata2")))
print(FeaturePlot(seurat, features = c("Mpo", "Csf1r", "S100a6", "S100a8", "Car1", "Ccr2", "Cd74", "Ly6c2", 
                                       "Cd4","Cd8a","Siglech","Fcer1a")))
dev.off()


# Save the object ---------------------------------------------------------
saveRDS(seurat,paste0(data.dir,Sys.Date(),"_bm_ptpn11_SCT_clustered.RDS"))

# Load the marker genes ---------------------------------------------------
clus.new <- openxlsx::read.xlsx("tables/2023-06-22_annotation.xlsx")
colSums(table(clus.new$cluster,clus.new$Layer2))
clus.new.sub <- clus.new[!duplicated(clus.new$cluster),]
duplicated(clus.new.sub)
use.ids <- clus.new.sub$Layer2
names(use.ids) <- clus.new.sub$cluster
seurat <- RenameIdents(seurat,use.ids)

# Save the object ---------------------------------------------------------
saveRDS(seurat,paste0(data.dir,Sys.Date(),"_bm_ptpn11_sct_annotated.RDS"))