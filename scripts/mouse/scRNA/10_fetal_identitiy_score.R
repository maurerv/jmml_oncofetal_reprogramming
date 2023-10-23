################################################
#
# 10_fetal_identity_score
# Author: Maximilian Schönung
# Date: 23.10.2023 
#
################################################

# Set the paths -----------------------------------------------------------

setwd("/omics/odcf/analysis/OE0565_projects/ptpn11/")
data.dir <- "/omics/odcf/analysis/OE0565_projects/ptpn11/data/"
plot.dir <- "/omics/odcf/analysis/OE0565_projects/ptpn11/plots/"

subDir <- "03_gsea/"
if (!file.exists(paste0(plot.dir,subDir))){
  dir.create(paste0(plot.dir,subDir),recursive = T)
}

# Load the libraries ------------------------------------------------------
library(Seurat)
library(ggplot2)
library(dplyr)
library(msigdbr)
library(fgsea)
library(RColorBrewer)
library(biomaRt)

# Load the 10x Data -------------------------------------------------------
seurat <- readRDS(paste0(data.dir,"2023-09-18_bm_ptpn11_sct_annotated.RDS"))
seurat.hsc.bm <- subset(seurat,idents = "HSC")
rm(seurat) #make the env lighter
included.genes <- rownames(seurat.hsc.bm@assays$RNA@counts)

# Remove sex chromosomes --------------------------------------------------
ensembl <- useMart('ensembl', dataset = 'mmusculus_gene_ensembl',host="http://nov2020.archive.ensembl.org")
genes <- getBM(
  attributes = c('external_gene_name','chromosome_name'),
  filters = 'chromosome_name',
  values = c("X","Y"),
  mart = ensembl)
tail(genes)
included.genes <- included.genes[!included.genes%in%genes$external_gene_name]
hsc.sub <- subset(seurat.hsc.bm,features=included.genes)

# Change the Identities for Diff Exp --------------------------------------
hsc.sub <- NormalizeData(hsc.sub)
Idents(hsc.sub) <- hsc.sub$leuk
diff.genes <- FindMarkers(hsc.sub, ident.1 = "1", ident.2 = "0",
                           logfc.threshold = -Inf, min.pct = 0.05, max.cells.per.ident = 1000)

sorted <- order(diff.genes$avg_log2FC,decreasing = T)
ranks <- diff.genes$avg_log2FC[sorted]
names(ranks) <- rownames(diff.genes)[sorted]

#get gene sets
fetal <- openxlsx::read.xlsx("tables/1-s2.0-S1934590920303970-mmc4.xlsx")
fgsea_sets <- as.list(as.data.frame(apply(fetal,2,stringr::str_to_title)))
fgsea_sets <- lapply(fgsea_sets,function(x)x[nzchar(x)])
fgsea_sets <- lapply(fgsea_sets,function(x)x[x%in%names(ranks)])

#get gene sets
mark_gset <- read.csv("tables/JMML.REF-HSC.DEG.Top50_mod.csv")
for(i in levels(as.factor(mark_gset$stage))){
  fgsea_sets[[i]] <- mark_gset[mark_gset$stage==i,"gene"]
}
fgsea_sets <- lapply(fgsea_sets,stringr::str_to_title)
fgsea_sets <- lapply(fgsea_sets,function(x)x[x%in%included.genes])

# fGSEA
fgseaRes<- fgsea(fgsea_sets, stats = ranks)
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))
write.table(as.data.frame(fgseaResTidy)[,1:7],paste0("tables/",Sys.Date(),"_mouse_gsea.txt"),sep="\t",quote = F)

# Sig genes ---------------------------------------------------------------
pdf(paste0(plot.dir,subDir,Sys.Date(),"_GSEA_mark_fetal_hsc_hpc.pdf"),width = 10)
ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill= NES < 0)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  ylim(-2.5,2.5)+
  theme_classic()+
  geom_hline(yintercept = 0)

ggplot(fgseaResTidy[fgseaResTidy$pathway%in%c("FCA_Sampled","HCA_CB_Sampled","HCA_BM_Sampled","JUV_Sampled"),],
       aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill= NES < 0)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  ylim(-2.5,2.5)+
  theme_classic()+
  geom_hline(yintercept = 0)

ggplot(fgseaResTidy[!fgseaResTidy$pathway%in%c("FCA_Sampled","HCA_CB_Sampled","HCA_BM_Sampled","JUV_Sampled"),],
       aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill= NES < 0)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  ylim(-2.5,2.5)+
  theme_classic()+
  geom_hline(yintercept = 0)

ggplot(fgseaResTidy[fgseaResTidy$pathway%in%c("Fetal.HSC.Identity.Genes","Adult.HSC.Identity.Genes"),],
       aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill= NES < 0)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  ylim(-2.5,2.5)+
  theme_classic()+
  geom_hline(yintercept = 0)

plotEnrichment(fgsea_sets[["FCA_Sampled"]],
               ranks)+
  theme_classic()+
  ggtitle("FCA_Sampled")

plotEnrichment(fgsea_sets[["HCA_CB_Sampled"]],
               ranks)+
  theme_classic()+
  ggtitle("HCA_CB_Sampled")

plotEnrichment(fgsea_sets[["HCA_BM_Sampled"]],
               ranks)+
  theme_classic()+
  ggtitle("HCA_BM_Sampled")

plotEnrichment(fgsea_sets[["JUV_Sampled"]],
               ranks)+
  theme_classic()+
  ggtitle("JUV_Sampled")

plotEnrichment(fgsea_sets[["Fetal.HSC.Identity.Genes"]],
               ranks)+
  theme_classic()+
  ggtitle("Fetal.HSC.Identity.Genes")

plotEnrichment(fgsea_sets[["Adult.HSC.Identity.Genes"]],
               ranks)+
  theme_classic()+
  ggtitle("Adult.HSC.Identity.Genes")

plotEnrichment(fgsea_sets[["Adult.HPC.Identity.Genes"]],
               ranks)+
  theme_classic()+
  ggtitle("Adult.HPC.Identity.Genes")

plotEnrichment(fgsea_sets[["Fetal.HPC.Identity.Genes"]],
               ranks)+
  theme_classic()+
  ggtitle("Fetal.HPC.Identity.Genes")
dev.off()


# Plot some examples ------------------------------------------------------
pdf(paste0(plot.dir,subDir,Sys.Date(),"_murine_fetal_adult_example_vlnPlts.pdf"),width = 10)
VlnPlot(seurat.hsc.bm,features=c("Cd52"),split.by = "leuk",cols =c("gray40","firebrick"),pt.size = F)
VlnPlot(seurat.hsc.bm,features=c("Lgals1","Plac8","Ran"),split.by = "leuk",cols =c("gray40","firebrick"),pt.size = F)
VlnPlot(seurat.hsc.bm,features=c("H3f3b","Fos","Procr"),split.by = "leuk",cols =c("gray40","firebrick"),pt.size = F)
VlnPlot(seurat.hsc.bm,features=c("Jun","Fos"),split.by = "leuk",cols =c("gray40","firebrick"),pt.size = F)
VlnPlot(seurat.hsc.bm,features=c("Mecom","H3f3b","Nfkbia"),split.by = "leuk",cols =c("gray40","firebrick"),pt.size =F)
VlnPlot(seurat.hsc.bm,features=c("Nfkbia"),split.by = "leuk",cols =c("gray40","firebrick"),pt.size = T)
VlnPlot(seurat.hsc.bm,features=c("Mef2c"),split.by = "leuk",cols =c("gray40","firebrick"),pt.size = T)
dev.off()

pdf(paste0(plot.dir,subDir,Sys.Date(),"_Cd52_VlnPlt_wBoxplot.pdf"),width = 10)
VlnPlot(seurat.hsc.bm,features=c("Cd52"),split.by = "leuk",cols =c("gray40","firebrick"),pt.size = F)+
  geom_boxplot(colour="white",width=.1,position=position_dodge(.9))
dev.off()



# Use Module Scores to verify the GSEA ------------------------------------
seurat <- AddModuleScore(seurat.hsc.bm,features=list(fgsea_sets$Fetal.HSC.Identity.Genes),name="FetalHSC")
seurat <- AddModuleScore(seurat,features=list(fgsea_sets$Adult.HSC.Identity.Genes),name="AdultHSC")
VlnPlot(seurat,features=c("FetalHSC1","AdultHSC1"),split.by = "leuk",cols =c("gray40","firebrick"),pt.size =.1)

#Test for fetal signature
LMM1.stem <- lme4::lmer(FetalHSC1~leuk+(1|orig.ident),data = seurat@meta.data)
LMM2.stem <- lme4::lmer(FetalHSC1~(1|orig.ident),data = seurat@meta.data)
anova(LMM2.stem,LMM1.stem)

#Test for adult signature
LMM1.stem <- lme4::lmer(AdultHSC1~leuk+(1|orig.ident),data = seurat@meta.data)
LMM2.stem <- lme4::lmer(AdultHSC1~(1|orig.ident),data = seurat@meta.data)
anova(LMM2.stem,LMM1.stem)

# Stats for diff expression of markers ------------------------------------
diff.genes[c("Lgals1","Plac8","H3f3b","Fos","Procr","Cd52"),]

