################################################
#
# 03_norm
# Author: Maximilian Schönung
# Date: 06.07.2022
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

# Load the 10x Data -------------------------------------------------------
seurat.filt <- readRDS(paste0(data.dir,"2023-06-19_bm_ptpn11_doublet_filtered.RDS"))

# Normalize ---------------------------------------------------
seurat.n <- NormalizeData(seurat.filt)

# Assign Cell Cycle State -------------------------------------------------
cell_cycle_markers <- read.delim(paste0(data.dir,"mus_musculus_cell_cycle_genes.txt"),stringsAsFactors = F)

# Acquire the S phase genes
s_genes <- cell_cycle_markers %>%
  dplyr::filter(phase == "S") %>%
  pull("gene_name")

# Acquire the G2M phase genes        
g2m_genes <- cell_cycle_markers %>%
  dplyr::filter(phase == "G2/M") %>%
  pull("gene_name")

seurat.n <- CellCycleScoring(seurat.n, 
                             g2m.features = g2m_genes, 
                             s.features = s_genes)


# Add the cell cycle state to the filtered seurat -------------------------
table(rownames(seurat.filt@meta.data)==rownames(seurat.n@meta.data)) #just if all true!
seurat.filt@meta.data <- seurat.n@meta.data

# SCT Integration ----------------------------------------------------
split_seurat <- SplitObject(seurat.filt, split.by = "leuk")
split_seurat <- split_seurat[c("0", "1")]
split_seurat <- lapply(split_seurat,SCTransform)

# Select the most variable features to use for integration
integ_features <- SelectIntegrationFeatures(object.list = split_seurat, nfeatures = 3000) 

# Prepare the SCT list object for integration
split_seurat <- PrepSCTIntegration(object.list = split_seurat, anchor.features = integ_features)

# Find best buddies - can take a while to run
integ_anchors <- FindIntegrationAnchors(object.list = split_seurat, normalization.method = "SCT", anchor.features = integ_features)

# Integrate across conditions
seurat_integrated <- IntegrateData(anchorset = integ_anchors, normalization.method = "SCT")

# save the integrated dataset
saveRDS(seurat_integrated,paste0(data.dir,Sys.Date(),"_bm_ptpn11_doubletFilt_sctransform_integrated.RDS"))