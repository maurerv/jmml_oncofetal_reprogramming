################################################
#
# 06_heatmap
# Author: Maximilian Schönung
# Date: 19.10.2023 
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
library(dplyr)
library(ggplot2)
library(viridis)

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
all(names(col.vec)%in%unique(Idents(seurat)))
all(unique(Idents(seurat))%in%names(col.vec))
umap <- data.frame(seurat@reductions$umap@cell.embeddings)
umap$cluster <- Idents(seurat)

pdf(paste0(plot.dir,subDir,Sys.Date(),"_umap_annotated_colours.pdf"),useDingbats = F)
DimPlot(seurat, reduction = "umap", label = TRUE, pt.size = 0.1,label.size = 3,cols=col.vec,shuffle=T) + NoLegend()
ggplot(umap,aes(UMAP_1,UMAP_2,fill=cluster,colour=cluster))+
  geom_point(shape=21,alpha=.5,size=.3)+
  scale_color_manual(values=col.vec)+
  scale_fill_manual(values=col.vec)+
  guides(color="none",fill="none")+
  theme_classic()
dev.off()


seurat.down <- subset(x = seurat, downsample = 25)
rm(seurat)
markers <- openxlsx::read.xlsx("tables/2023-06-22_annotation.xlsx")
markers.sel <- markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
seurat.down2 <- seurat.down
seurat.down2@active.ident <- factor(as.character(seurat.down2@active.ident),
                                    levels=c("HSC","MPP1","MPP2",
                                             "MPP_Cycling1","MPP_Cycling2",
                                             "MPP_EryMeg",
                                             "MyP1","MyP2","MonoP","NeutroP",
                                             "Mono2","Mono1","DC",
                                             "Neutro4","Neutro1","Neutro3","Neutro5","Neutro2",
                                             "MEP","Ery1","Ery2",
                                             "MPP_Lympho","LyP","T-cell",
                                             "B-cell1","B-cell2","B-cell3","B-cell4",
                                             "Basophils"))

# Order the features ------------------------------------------------------
markers.sun <- markers.sel
ind <- which(!duplicated(seurat.down@meta.data$seurat_clusters))
clust.df <- data.frame("clust"=as.character(Idents(seurat.down))[ind],"ind"=seurat.down@meta.data$seurat_clusters[ind])
rownames(clust.df) <- clust.df$clust
ind2 <- clust.df[levels(seurat.down2@active.ident),"ind"]
list.marker <- list()
for(i in ind2){
  list.marker[[i]] <- markers.sel[markers.sel$cluster==i,]
}
marker.plot <- do.call(rbind,list.marker)
marker.plot <- unique(marker.plot$gene)

cols <- viridis(100)[c(1, 50, 100)]
pdf(paste0(plot.dir,subDir,Sys.Date(),"_marker_heatmap.pdf"))
DoHeatmap(seurat.down2,
          marker.plot,
          size=2,group.colors = col.vec)+
  theme(text = element_text(size = 4))+
  scale_fill_viridis()+
  guides(color="none")
dev.off()