################################################
#
# 02_doublet_finder
# Author: Maximilian Schönung
# Date: 19.07.2023
#
################################################
#https://rpubs.com/olneykimberly/Ecoli_brain_processing

# Set the paths -----------------------------------------------------------
setwd("/omics/odcf/analysis/OE0565_projects/ptpn11/")
data.dir <- "/omics/odcf/analysis/OE0565_projects/ptpn11/data/"
plot.dir <- "/omics/odcf/analysis/OE0565_projects/ptpn11/plots/"

# Load the libraries ------------------------------------------------------
library(Seurat)
library(ggplot2)
library(dplyr)
library(DoubletFinder)

# Load the 10x Data -------------------------------------------------------
seurat.filt <- readRDS(paste0(data.dir,"2023-06-19_bm_ptpn11_combined_filt_strict.RDS"))

# Split the object --------------------------------------------------------
# DoubletFinder needs to be performed on a single sample basis
seurat.split <- SplitObject(seurat.filt, split.by = "orig.ident") 

# Doublet Finder function -------------------------------------------------
double_fun <- function(x){
  norm <- NormalizeData(x)
  norm <- FindVariableFeatures(norm, selection.method = "vst", nfeatures = 2000)
  norm <- ScaleData(norm)  
  norm <- RunPCA(norm)
  norm <- RunUMAP(norm, dims = 1:30, reduction = "pca")
  
  #start doublet finder
  sweep.res <- paramSweep_v3(norm, PCs = 1:30, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  
  # find pK
  bcmvn_max <- bcmvn[which.max(bcmvn$BCmetric),]
  pK_value <- bcmvn_max$pK
  pK_value <- as.numeric(levels(pK_value))[pK_value]
  
  #expect around 8% of doublets
  nExp_poi <- round(0.08*nrow(norm@meta.data))
  
  #doublet findet
  norm <- doubletFinder_v3(norm, PCs = 1:30, 
                           pN = 0.25, pK = pK_value, nExp = nExp_poi, 
                           reuse.pANN = FALSE, sct = FALSE)
  colnames(norm@meta.data)[length(colnames(norm@meta.data))] <- "DoubletClassification"
  return(norm)
}

# Apply Doublet Finder Function -------------------------------------------
seurat.doublets <- lapply(seurat.split,double_fun)

subDir <- "00_qc/"
if (!file.exists(paste0(plot.dir,subDir))){
  dir.create(paste0(plot.dir,subDir),recursive = T)
}
pdf(paste0(plot.dir,subDir,Sys.Date(),"_doublet_umap.pdf"))
lapply(seurat.doublets,function(x)print(DimPlot(x,reduction = "umap",group.by = "DoubletClassification")))
dev.off()

# Get the doublet table ---------------------------------------------------
doublets.list <- lapply(seurat.doublets,function(x)data.frame("cellID"=rownames(x@meta.data),"id"=x@meta.data$orig.id,"Doublets"=x@meta.data$DoubletClassification))

doublets.df <- do.call(rbind,doublets.list)
plot.doublets <- as.data.frame(table(doublets.df$id,doublets.df$Doublets))
plot.doublets$Var1 <- factor(plot.doublets$Var1,c("345BM","347BM","349BM","348BM","351BM","352BM"))

pdf(paste0(plot.dir,subDir,Sys.Date(),"_doublet_bar.pdf"))
ggplot(plot.doublets,aes(Var1,Freq,fill=Var2))+
  geom_bar(stat="identity",position = "dodge")+
  theme_classic()+
  xlab("Sample")+
  ylab("Number of Cells")
ggplot(plot.doublets,aes(Var1,log10(Freq),fill=Var2))+
  geom_bar(stat="identity",position = "dodge")+
  theme_classic()+
  xlab("Sample")+
  ylab("Number of Cells [log10]")
dev.off()

saveRDS(doublets.df,paste0(data.dir,Sys.Date(),"_doublet_information.RDS"))


# Remove the doublets and save the seurat object --------------------------
doublets <- doublets.df[doublets.df$Doublets=="Doublet","cellID"]
seurat.filt2 <- seurat.filt[,(!colnames(seurat.filt)%in%doublets)]

saveRDS(seurat.filt2,paste0(data.dir,Sys.Date(),"_bm_ptpn11_doublet_filtered.RDS"))