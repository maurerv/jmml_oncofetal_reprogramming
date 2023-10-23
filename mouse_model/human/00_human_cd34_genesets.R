################################################
#
# 00_human_cd34_genesets.R
# Author: Maximilian Schönung
# Date: 12.10.2023 
#
################################################

# Set the paths -----------------------------------------------------------

setwd("/omics/odcf/analysis/OE0565_projects/ptpn11/")
data.dir <- "/omics/odcf/analysis/OE0219_projects/JMMLC/scRNA_Data/2.Single_Dataset/JMML_CD34/29/"
plot.dir <- "/omics/odcf/analysis/OE0565_projects/ptpn11/plots/"
table.dir <- "/omics/odcf/analysis/OE0565_projects/ptpn11/tables/"

subDir <- "04_human/"
if (!file.exists(paste0(plot.dir,subDir))){
  dir.create(paste0(plot.dir,subDir),recursive = T)
}

# Load the libraries ------------------------------------------------------
library(Seurat)
library(ggplot2)
library(dplyr)
library(fgsea)
library(org.Hs.eg.db)
library(RColorBrewer)

# Load the 10x Data -------------------------------------------------------
human <- readRDS(paste0(data.dir,"JMML_CD34_PC_29.SeuratV3.rds"))
DefaultAssay(human) <- "RNA"
human@meta.data$cell_id <- rownames(human@meta.data)
meta <- readRDS("/omics/odcf/analysis/OE0565_projects/ptpn11/data_human/10X_meta.data_reduced.rds")
human@meta.data <- cbind(human@meta.data,meta[rownames(human@meta.data),])
human <- subset(human,subset = predicted.id == "CD34+ HSC")

# Load the genesets -------------------------------------------------------
geneset <- openxlsx::read.xlsx(paste0(table.dir,"Maxi_human_gene_sets_v2.xlsx"))
mat <- as.matrix(GetAssayData(human, slot = "counts",assay = "RNA"))
annots <- select(org.Hs.eg.db, keys=rownames(mat), 
                 columns="SYMBOL", keytype="ENSEMBL")
annots <- annots[!duplicated(annots$ENSEMBL),]
rownames(annots) <- annots$ENSEMBL

# Density Plot Risk ----------------------------------------------------
selected.set <- as.character(na.omit(geneset$Dick_17LSC))
seurat <- AddModuleScore(human,features=list(annots$ENSEMBL[annots$SYMBOL%in%na.omit(selected.set)]),
                         name="Dick_17LSC")
selected.set <- as.character(na.omit(geneset$Bresolin_JMML_AML_like))
seurat <- AddModuleScore(seurat,features=list(annots$ENSEMBL[annots$SYMBOL%in%na.omit(selected.set)]),
                         name="Bresolin_JMML_AML_like")
selected.set <- as.character(na.omit(geneset$`Pellin_CD34_2_3_MEP_Ly-My_MEP`))
seurat <- AddModuleScore(seurat,features=list(annots$ENSEMBL[annots$SYMBOL%in%na.omit(selected.set)]),
                         name="MkEryPriming")
selected.set <- as.character(na.omit(geneset$`Pellin_CD34_2_3_MEP_Ly-My_Ly-My`))
seurat <- AddModuleScore(seurat,features=list(annots$ENSEMBL[annots$SYMBOL%in%na.omit(selected.set)]),
                         name="MyeloidPriming")
sets <- data.frame("Dick_17LSC"=seurat@meta.data$Dick_17LSC1,
                   "Bresolin_JMML_AML_like"=seurat@meta.data$Bresolin_JMML_AML_like1,
                   "MkEryPriming"=seurat@meta.data$MkEryPriming1,
                   "MyeloidPriming"=seurat@meta.data$MyeloidPriming1,
                   "Subgroup"=seurat@meta.data$Epigenotype,
                   "Donor"=factor(seurat@meta.data$orig.ident),
                   "Sex"=factor(seurat@meta.data$Sex),
                   "Mutation"=factor(seurat@meta.data$Genotype))
summary(sets)

## subset to same cell number
set.seed(1234)
sub.list <- list()
for(j in levels(factor(sets$Subgroup))){
  sub.list[[j]] <- sample(which(sets$Subgroup==j),1000)
}
sub.vec <- unlist(sub.list)
sets.plot <- sets[sub.vec,]
sets.plot$Subgroup <- factor(sets.plot$Subgroup,levels=c("LM","IM","HM"))

pdf(paste0(plot.dir,subDir,Sys.Date(),"_gene_modules_risk.pdf"),height = 4)
ggplot(sets.plot,
       aes(Dick_17LSC,Bresolin_JMML_AML_like,fill=Subgroup))+
  stat_density_2d(geom = "polygon", aes(alpha = ..level.., fill = Subgroup))+
  facet_wrap(~Subgroup,nrow=1)+
  scale_fill_manual(values=c(HM ="#c33126", IM = "#fbbb25", LM = "#0058b4"))+
  theme_classic()+
  geom_hline(yintercept = median(sets.plot[sets.plot$Subgroup=="LM","Bresolin_JMML_AML_like"]),linetype="dotted",color="black")+
  geom_vline(xintercept = median(sets.plot[sets.plot$Subgroup=="LM","Dick_17LSC"]),linetype="dotted",color="black")
ggplot(sets.plot,
       aes(MkEryPriming,MyeloidPriming,fill=Subgroup))+
  stat_density_2d(geom = "polygon", aes(alpha = ..level.., fill = Subgroup))+
  facet_wrap(~Subgroup,nrow=1)+
  scale_fill_manual(values=c(HM ="#c33126", IM = "#fbbb25", LM = "#0058b4"))+
  theme_classic()+
  geom_hline(yintercept = median(sets.plot[sets.plot$Subgroup=="LM","MyeloidPriming"]),linetype="dotted",color="black")+
  geom_vline(xintercept = median(sets.plot[sets.plot$Subgroup=="LM","MkEryPriming"]),linetype="dotted",color="black")
dev.off()


# LMM and ANOVA -----------------------------------------------------------
#Test for fetal signature
LMM1.stem <- lme4::lmer(Dick_17LSC~Subgroup+(1|Donor),data = sets.plot)
LMM2.stem <- lme4::lmer(Dick_17LSC~(1|Donor),data = sets.plot)
anova(LMM2.stem,LMM1.stem)

LMM1.stem <- lme4::lmer(Bresolin_JMML_AML_like~Subgroup+(1|Donor),data = sets.plot)
LMM2.stem <- lme4::lmer(Bresolin_JMML_AML_like~(1|Donor),data = sets.plot)
anova(LMM2.stem,LMM1.stem)

LMM1.stem <- lme4::lmer(MyeloidPriming~Subgroup+(1|Donor),data = sets.plot)
LMM2.stem <- lme4::lmer(MyeloidPriming~(1|Donor),data = sets.plot)
anova(LMM2.stem,LMM1.stem)

LMM1.stem <- lme4::lmer(MkEryPriming~Subgroup+(1|Donor),data = sets.plot)
LMM2.stem <- lme4::lmer(MkEryPriming~(1|Donor),data = sets.plot)
anova(LMM2.stem,LMM1.stem)



