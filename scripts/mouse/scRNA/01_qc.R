################################################
#
# 01_qc
# Author: Maximilian Schönung
# Date: 19.07.2023
#
################################################


# Set the paths -----------------------------------------------------------
setwd("/omics/odcf/analysis/OE0565_projects/ptpn11/")
data.dir <- "/omics/odcf/analysis/OE0565_projects/ptpn11/data/"
plot.dir <- "/omics/odcf/analysis/OE0565_projects/ptpn11/plots/"

# Load the libraries ------------------------------------------------------
library(Seurat)
library(ggplot2)
library(RColorBrewer)

# Load the 10x Data -------------------------------------------------------
leuk <- readRDS(paste0(data.dir,"2023-06-19_bm_ptpn11_combined_raw.RDS"))

# Add meta data and mito to qc metrix ---------------------------------------------------
leuk@meta.data$leuk <- as.numeric(leuk@meta.data$orig.ident%in%c("345BM","347BM","349BM"))
leuk@meta.data$id <- substr(leuk@meta.data$orig.ident, 0, 3)
leuk@meta.data$sex <- rep(NA,length(leuk@meta.data$id))
leuk@meta.data$sex[leuk@meta.data$id%in%c(351,352)] <- "M"
leuk@meta.data$sex[leuk@meta.data$id%in%c(345,347,348,349)] <- "F"
leuk@meta.data$batch <- rep(NA,length(leuk@meta.data$id))
leuk@meta.data$batch[leuk@meta.data$orig.ident%in%c("345BM","347BM","348BM","349BM","351BM")] <- 1
leuk@meta.data$batch[leuk@meta.data$orig.ident%in%c("352BM","345S","347S","349S")] <- 2
leuk[["percent.mt"]] <- PercentageFeatureSet(leuk, pattern = "^mt-")


# Inspect the meta data ---------------------------------------------------
meta <- leuk@meta.data

meta$orig.ident <- factor(meta$orig.ident)
meta.melt <- reshape2::melt(meta,id.vars=c("orig.ident","leuk","id"))


qc1 <- ggplot(meta,aes(log10(nCount_RNA),log10(nFeature_RNA),color=percent.mt>=5))+
  geom_point()+
  theme_classic()+
  geom_rug(alpha=.2,color="blue")+
  ggtitle("Percent Mitochondria 5% Cutoff")

qc2 <- ggplot(meta,aes(log10(nCount_RNA),log10(nFeature_RNA),color=percent.mt>=10))+
  geom_point()+
  theme_classic()+
  geom_rug(alpha=.2,color="blue")+
  ggtitle("Percent Mitochondria 10% Cutoff")

qc3 <- ggplot(meta,aes(orig.ident,nFeature_RNA,fill=orig.ident))+
  geom_violin()+
  geom_jitter(alpha=0.1)+
  theme_classic()+
  ggtitle("Number Counts per Cell")+
  scale_fill_brewer(palette = "Set1")

qc3v <- ggplot(meta,aes(orig.ident,nFeature_RNA,fill=orig.ident))+
  geom_violin()+
  theme_classic()+
  ggtitle("Number Counts per Cell")+
  scale_fill_brewer(palette = "Set1")

qc3l <- ggplot(meta,aes(orig.ident,nFeature_RNA,fill=orig.ident))+
  geom_violin()+
  geom_jitter(alpha=0.1)+
  theme_classic()+
  ggtitle("Number Counts per Cell")+
  geom_hline(yintercept = 500,linetype="dotted",color="red")+
  geom_hline(yintercept = 5800,linetype="dotted",color="red")+
  scale_fill_brewer(palette = "Set1")

qc4 <- ggplot(meta,aes(orig.ident,nCount_RNA,fill=orig.ident))+
  geom_violin()+
  geom_jitter(alpha=0.1)+
  theme_classic()+
  ggtitle("Number Genes per Cell")+
  scale_fill_brewer(palette = "Set1")


qc4v <- ggplot(meta,aes(orig.ident,nCount_RNA,fill=orig.ident))+
  geom_violin()+
  theme_classic()+
  ggtitle("Number Genes per Cell")+
  scale_fill_brewer(palette = "Set1")

qc5 <- ggplot(meta,aes(orig.ident,percent.mt,fill=orig.ident))+
  geom_violin()+
  geom_jitter(alpha=0.3)+
  theme_classic()+
  ggtitle("Percentage Mito")

qc5l <- ggplot(meta,aes(orig.ident,percent.mt,fill=orig.ident))+
  geom_violin()+
  geom_jitter(alpha=0.3)+
  theme_classic()+
  ggtitle("Percentage Mito; Cutoff 5%")+
  geom_hline(yintercept = 5,linetype="dotted",color="red")

qc5l2 <- ggplot(meta,aes(orig.ident,percent.mt,fill=orig.ident))+
  geom_violin()+
  geom_jitter(alpha=0.3)+
  theme_classic()+
  ggtitle("Percentage Mito; Cutoff 10%")+
  geom_hline(yintercept = 10,linetype="dotted",color="red")



# Create QC Subdirectory --------------------------------------------------
subDir <- "00_qc/"
if (!file.exists(paste0(plot.dir,subDir))){
  dir.create(paste0(plot.dir,subDir))
}

pdf(paste0(plot.dir,subDir,Sys.Date(),"_qc_plot.pdf"))
qc1
qc2
qc3
qc3v
qc3l
qc4
qc4v
qc5
qc5l
qc5l2
dev.off()


table(leuk$percent.mt<5&leuk$nFeature_RNA>500&leuk$nFeature_RNA<5800)
# removes 2821 cells and retains 31 600 (for 5% mito)

# Filter the set ----------------------------------------------------------
filt <- subset(leuk, subset = nFeature_RNA > 500 & nFeature_RNA < 5800 & percent.mt < 5)

# Exclude Genes which are just expressed in less than 11 cells ----------------------
counts <- GetAssayData(object = filt, slot = "counts")
filtered_seurat <- CreateSeuratObject(counts, meta.data = filt@meta.data, min.cells = 10)

# Save the dataset --------------------------------------------------------
saveRDS(filtered_seurat,paste0(data.dir,Sys.Date(),"_bm_ptpn11_combined_filt_strict.RDS"))
