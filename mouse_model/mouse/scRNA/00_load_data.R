################################################
#
# 00_Combine 10x Datasets (just BM)
# Author: Maximilian Schönung
# Date: 19.07.2023 
#
################################################


# Set the paths -----------------------------------------------------------
setwd("/omics/odcf/analysis/OE0565_projects/ptpn11/")
data.dir <- "/omics/odcf/analysis/OE0565_projects/ptpn11/data/"

# Load the libraries ------------------------------------------------------
library(Seurat)

# Load the 10x Data -------------------------------------------------------
bm.345_data <- Read10X(data.dir = "/omics/odcf/project/OE0565/ptpn11/sequencing/10x_scRNA_sequencing/view-by-pid/OE0565_PTPN11_345/bone-marrow01/paired/merged-alignment/RG_hg_GRCm38-2020-A_TV_CELL_RANGER-1.2.0_EC_-_FC_-_PV_cellranger-4.0.0_ID_0/filtered_feature_bc_matrix/")
bm.347_data <- Read10X(data.dir = "/omics/odcf/project/OE0565/ptpn11/sequencing/10x_scRNA_sequencing/view-by-pid/OE0565_PTPN11_347/bone-marrow01/paired/merged-alignment/RG_hg_GRCm38-2020-A_TV_CELL_RANGER-1.2.0_EC_-_FC_-_PV_cellranger-4.0.0_ID_0/filtered_feature_bc_matrix/")
bm.348_data <- Read10X(data.dir = "/omics/odcf/project/OE0565/ptpn11/sequencing/10x_scRNA_sequencing/view-by-pid/OE0565_PTPN11_348/bone-marrow01/paired/merged-alignment/RG_hg_GRCm38-2020-A_TV_CELL_RANGER-1.2.0_EC_-_FC_-_PV_cellranger-4.0.0_ID_0/filtered_feature_bc_matrix/")
bm.349_data <- Read10X(data.dir = "/omics/odcf/project/OE0565/ptpn11/sequencing/10x_scRNA_sequencing/view-by-pid/OE0565_PTPN11_349/bone-marrow01/paired/merged-alignment/RG_hg_GRCm38-2020-A_TV_CELL_RANGER-1.2.0_EC_-_FC_-_PV_cellranger-4.0.0_ID_0/filtered_feature_bc_matrix/")
bm.351_data <- Read10X(data.dir = "/omics/odcf/project/OE0565/ptpn11/sequencing/10x_scRNA_sequencing/view-by-pid/OE0565_PTPN11_351/bone-marrow01/paired/merged-alignment/RG_hg_GRCm38-2020-A_TV_CELL_RANGER-1.2.0_EC_-_FC_-_PV_cellranger-4.0.0_ID_0/filtered_feature_bc_matrix/")
bm.352_data <- Read10X(data.dir = "/omics/odcf/project/OE0565/ptpn11/sequencing/10x_scRNA_sequencing/view-by-pid/OE0565_PTPN11_352/bone-marrow01/paired/merged-alignment/RG_hg_GRCm38-2020-A_TV_CELL_RANGER-1.2.0_EC_-_FC_-_PV_cellranger-4.0.0_ID_0/filtered_feature_bc_matrix/")

# Create Seurat Datasets from the Data ------------------------------------
bm.345 <- CreateSeuratObject(counts = bm.345_data, project = "345BM")
bm.347 <- CreateSeuratObject(counts = bm.347_data, project = "347BM")
bm.348 <- CreateSeuratObject(counts = bm.348_data, project = "348BM")
bm.349 <- CreateSeuratObject(counts = bm.349_data, project = "349BM")
bm.351 <- CreateSeuratObject(counts = bm.351_data, project = "351BM")
bm.352 <- CreateSeuratObject(counts = bm.352_data, project = "352BM")


# Combine the datasets ----------------------------------------------------
leuk <- merge(bm.345,y=c(bm.347,bm.348,bm.349,bm.351,bm.352), 
              add.cell.ids = c("345BM","347BM","348BM","349BM","351BM","352BM"), project = "PTPN11")
saveRDS(leuk,paste0(data.dir,Sys.Date(),"bm_ptpn11_combined_raw.RDS"))
