#cell type prediction with bock celltype predictor
#test only with bock data
# rm -r ~/icgc/jmmlc_pbat/data/odcf_md/analysis/200904_Celltype_Prediction_Bock/predict_Celltypes_testBock
# Rscript ~/icgc/jmmlc_pbat/data/Bock_Celltype_Predictor/cell_type_predictor_code/cellTypePredictor.R \
#     --features ~/icgc/jmmlc_pbat/data/Bock_Celltype_Predictor/cell_type_predictor_input/methMatrix.tsv \
#     --t ~/icgc/jmmlc_pbat/data/Bock_Celltype_Predictor/cell_type_predictor_input/sampleAnnot.tsv \
#     --out ~/icgc/jmmlc_pbat_normals/data/odcf_md/analysis/200904_Celltype_Prediction_Bock/predict_Celltypes_testBock

#Libraries
library(DSS)
library(bsseq)
library(doParallel)
library(ChIPseeker)
library(foreach)
library(SummarizedExperiment)
library(doMC)
library(rtracklayer)
library(HDF5Array)
library(org.Mm.eg.db)
library(ggpubr)
library(randomcoloR)
library(RColorBrewer)
library(VennDiagram)
library(ChIPpeakAnno)
library(dendextend)
library(pheatmap)
library(data.table)
#Directories
input.dir <- "~/icgc/jmmlc_pbat_normals/data/odcf_md/analysis/"
analysis.dir <- "~/icgc/jmmlc_pbat_normals/data/odcf_md/analysis/Celltype_Prediction_Bock"
dir.create(analysis.dir)

#with own normal reference and JmmL
#load data
bsseq_all <- readRDS(file.path(input.dir ,"bsseq", "bsseq_HSC_combTotal_snpRemoved.rds"))
pheno <- pData(bsseq_all)

#read cell type region tation to extract methylation values
library(data.table)
input_regions <- fread("~/icgc/jmmlc_pbat/data/cell_type_predictor_input/regionAnnot.tsv")
input_regions <- makeGRangesFromDataFrame(as.data.frame(input_regions))
# liftOver
#get lift
ch <- import.chain("/C010-datasets/Internal/Joschka/data/hg38ToHg19.over.chain")
#perform liftOver
input_regions_19 <- liftOver(input_regions, ch)
#input_regions_19_red <- lapply(input_regions_19, reduce)
#select regions with a unique liftover
idx <- elementNROWS(input_regions_19)
idx <- which(idx == 1)
input_regions_19 <- unlist(input_regions_19[idx,])
length(input_regions_19)
length(idx)

#extract methylation levels
meth_predictorInput <- bsseq::getMeth(bsseq_all, regions= input_regions_19, type = "raw", what=c("perRegion"))

#prepare pheno
pheno$group <-as.character(pheno$Epigenotype)
pheno[pheno$Sample_Type=="normal", ]$group <-paste0(as.character(pheno[pheno$Sample_Type=="normal", ]$Celltype), "_", as.character(pheno[pheno$Sample_Type=="normal", ]$Tissue))
pheno_new <- data.frame(sampleName=rownames(pheno), class=pheno$group,
 donorId=pheno$Patient, nominalCellNumber=100, cellSourceCurated="JMML",
    train=ifelse(pheno$Sample_Type=="normal", TRUE, FALSE)
    )
#export 
write.table(meth_predictorInput, file.path(analysis.dir, "meth_HSC_comb_snpRemoved_repMerged_sub_cbHSC.tsv"), 
    quote=FALSE, sep="\t", row.names=FALSE)
write.table(pheno_new, file.path(analysis.dir, "sample__HSC_comb_snpRemoved_repMerged_sub_cbHSC.tsv"), 
    quote=FALSE, sep="\t", row.names=FALSE)   

#run commandline cell type tation
#rm -r ~/icgc/jmmlc_pbat_normals/data/odcf_md/analysis/Celltype_Prediction_Bock/predict_Celltypes_JMML_HSCUlm
system(command = "Rscript ~/icgc/jmmlc_pbat/data/Bock_Celltype_Predictor/cell_type_predictor_code/cellTypePredictor.R \\
--features ~/icgc/jmmlc_pbat_normals/data/odcf_md/analysis/Celltype_Prediction_Bock/meth_HSC_comb_snpRemoved_repMerged_sub_cbHSC.tsv \\
--annot ~/icgc/jmmlc_pbat_normals/data/odcf_md/analysis/Celltype_Prediction_Bock/sample__HSC_comb_snpRemoved_repMerged_sub_cbHSC.tsv \\
--out ~/icgc/jmmlc_pbat_normals/data/odcf_md/analysis/Celltype_Prediction_Bock/predict_Celltypes_JMML_HSCUlm"
)



#with Bock normal reference and JmmL
analysis.dir <- "~/icgc/jmmlc_pbat_normals/data/odcf_md/analysis/Celltype_Prediction_BockTrain_JMML"
dir.create(analysis.dir)
#read cell type region tation to extract methylation values
input_regions <- fread("~/icgc/jmmlc_pbat/data/cell_type_predictor_input/regionAnnot.tsv")
input_regions <- makeGRangesFromDataFrame(as.data.frame(input_regions))
# liftOver
#get lift
ch <- import.chain("/C010-datasets/Internal/Joschka/data/hg38ToHg19.over.chain")
#perform liftOver
input_regions_19 <- liftOver(input_regions, ch)
#input_regions_19_red <- lapply(input_regions_19, reduce)
#select regions with a unique liftover
idx <- elementNROWS(input_regions_19)
idx <- which(idx == 1)
input_regions_19 <- unlist(input_regions_19[idx,])
length(input_regions_19)
length(idx)

#extract methylation levels
meth_predictorInput_jmml <- bsseq::getMeth(bsseq_all, regions= input_regions_19, type = "raw", what=c("perRegion"))

#read bock normal references
#methylation
meth_bock <- fread("~/icgc/jmmlc_pbat/data/Bock_Celltype_Predictor/cell_type_predictor_input/methMatrix.tsv")
meth_bock<- as.data.frame(meth_bock)
#subset regions that were used for liftover
meth_bock_sub <- meth_bock[idx,]
annno_bock <- fread("~/icgc/jmmlc_pbat/data/Bock_Celltype_Predictor/cell_type_predictor_input/sampleAnnot.tsv")
annno_bock<- as.data.frame(annno_bock)

#subset only train data
annno_bock_sub <- annno_bock[annno_bock$train ==TRUE,]
dim(annno_bock)
dim(annno_bock_sub)
meth_bock_sub_train <- meth_bock_sub[, colnames(meth_bock_sub) %in% annno_bock_sub$sampleName]
dim(meth_bock_sub_train)

#get sample tation with matching columns 
pheno$group <-as.character(pheno$Epigenotype)
pheno[pheno$Sample_Type=="normal", ]$group <-paste0(as.character(pheno[pheno$Sample_Type=="normal", ]$Celltype), "_", as.character(pheno[pheno$Sample_Type=="normal", ]$Tissue))
pheno_new <- data.frame(sampleName=rownames(pheno), class=pheno$group,
 donorId=pheno$Patient, nominalCellNumber=100, cellSourceCurated="JMML",
    train=ifelse(pheno$Sample_Type=="normal", TRUE, FALSE)
    )
pheno_new_tumor <- pheno_new[pheno_new$train==FALSE,]
meth_predictorInput_jmml_tumor <- meth_predictorInput_jmml[,as.character(pheno_new_tumor$sampleName)]

#merge bock with own data
meth_predictorInput <- cbind(meth_bock_sub_train,meth_predictorInput_jmml_tumor )
pheno_predictorInput <- rbind(annno_bock_sub[, colnames(pheno_new_tumor)], pheno_new_tumor)
write.table(pheno_predictorInput, file.path(analysis.dir, "sample_anno.tsv"),
    quote=FALSE, sep="\t", row.names=FALSE)
    
write.table(meth_predictorInput, file.path(analysis.dir, "meth.tsv"), 
    quote=FALSE, sep="\t", row.names=FALSE)
  
#run commandline cell type tation
#rm -r ~/icgc/jmmlc_pbat_normals/data/odcf_md/analysis/Celltype_Prediction_BockTrain_JMML/predictCelltype_TrainBock_ClassJMML/
system(command = "Rscript ~/icgc/jmmlc_pbat/data/Bock_Celltype_Predictor/cell_type_predictor_code/cellTypePredictor.R \\
--features ~/icgc/jmmlc_pbat_normals/data/odcf_md/analysis/Celltype_Prediction_BockTrain_JMML/meth.tsv \\
--annot ~/icgc/jmmlc_pbat_normals/data/odcf_md/analysis/Celltype_Prediction_BockTrain_JMML/sample_anno.tsv \\
--out ~/icgc/jmmlc_pbat_normals/data/odcf_md/analysis/Celltype_Prediction_BockTrain_JMML/predictCelltype_TrainBock_ClassJMML/"
)
# rm -r ~/icgc/jmmlc_pbat_normals/data/odcf_md/analysis/Celltype_Prediction_BockTrain_JMML/predictCelltype_TrainBock_ClassJMML/
# Rscript ~/icgc/jmmlc_pbat/data/Bock_Celltype_Predictor/cell_type_predictor_code/cellTypePredictor.R \
# --features ~/icgc/jmmlc_pbat_normals/data/odcf_md/analysis/Celltype_Prediction_BockTrain_JMML/meth.tsv \
# --annot ~/icgc/jmmlc_pbat_normals/data/odcf_md/analysis/Celltype_Prediction_BockTrain_JMML/sample_anno.tsv \
# --out ~/icgc/jmmlc_pbat_normals/data/odcf_md/analysis/Celltype_Prediction_BockTrain_JMML/predictCelltype_TrainBock_ClassJMML/



#with Bock normal reference all and JmmL
analysis.dir <- "~/icgc/jmmlc_pbat_normals/data/odcf_md/analysis/Celltype_Prediction_BockTrain_JMML_All"

dir.create(analysis.dir)
#read cell type region tation to extract methylation values
input_regions <- fread("~/icgc/jmmlc_pbat/data/cell_type_predictor_input/regionAnnot.tsv")
input_regions <- makeGRangesFromDataFrame(as.data.frame(input_regions))
# liftOver
#get lift
ch <- import.chain("/C010-datasets/Internal/Joschka/data/hg38ToHg19.over.chain")
#perform liftOver
input_regions_19 <- liftOver(input_regions, ch)
#input_regions_19_red <- lapply(input_regions_19, reduce)
#select regions with a unique liftover
idx <- elementNROWS(input_regions_19)
idx <- which(idx == 1)
input_regions_19 <- unlist(input_regions_19[idx,])
length(input_regions_19)
length(idx)

#extract methylation levels
meth_predictorInput_jmml <- bsseq::getMeth(bsseq_all, regions= input_regions_19, type = "raw", what=c("perRegion"))

#read bock normal references
#methylation
meth_bock <- fread("~/icgc/jmmlc_pbat/data/Bock_Celltype_Predictor/cell_type_predictor_input/methMatrix.tsv")
meth_bock<- as.data.frame(meth_bock)
#subset regions that were used for liftover
meth_bock_sub <- meth_bock[idx,]
annno_bock <- fread("~/icgc/jmmlc_pbat/data/Bock_Celltype_Predictor/cell_type_predictor_input/sampleAnnot.tsv")
annno_bock<- as.data.frame(annno_bock)
annno_bock$class <- paste0(annno_bock$class, "_",annno_bock$cellSourceCurated )

#subset  train data
annno_bock_sub <- annno_bock
annno_bock_sub$train<- TRUE
dim(annno_bock)
dim(annno_bock_sub)
#subset classes with less than 5 samples
temp <- names(table(annno_bock_sub$class)[table(annno_bock_sub$class)> 5])
annno_bock_sub <- annno_bock_sub[annno_bock_sub$class %in% temp,]
dim(annno_bock_sub)
meth_bock_sub <- meth_bock[, colnames(meth_bock) %in% annno_bock_sub$sampleName]

#get sample tation with matching columns 
pheno$group <-as.character(pheno$Epigenotype)
pheno[pheno$Sample_Type=="normal", ]$group <-paste0(as.character(pheno[pheno$Sample_Type=="normal", ]$Celltype), "_", as.character(pheno[pheno$Sample_Type=="normal", ]$Tissue))
pheno_new <- data.frame(sampleName=rownames(pheno), class=pheno$group,
 donorId=pheno$Patient, nominalCellNumber=100, cellSourceCurated="JMML",
    train=ifelse(pheno$Sample_Type=="normal", TRUE, FALSE)
    )
pheno_new_tumor <- pheno_new[pheno_new$train==FALSE,]
meth_predictorInput_jmml_tumor <- meth_predictorInput_jmml[,as.character(pheno_new_tumor$sampleName)]

#merge bock with own data
meth_predictorInput <- cbind(meth_bock_sub_train,meth_predictorInput_jmml_tumor )
pheno_predictorInput <- rbind(annno_bock_sub[, colnames(pheno_new_tumor)], pheno_new_tumor)

 write.table(pheno_predictorInput, file.path(analysis.dir, "sample_anno.tsv"), 
    quote=FALSE, sep="\t", row.names=FALSE)
    
write.table(meth_predictorInput, file.path(analysis.dir, "meth.tsv"), 
    quote=FALSE, sep="\t", row.names=FALSE)
  
#run commandline cell type tation
#rm -r ~/icgc/jmmlc_pbat_normals/data/odcf_md/analysis/Celltype_Prediction_BockTrain_JMML_All/predictCelltype_TrainBock_ClassJMML/
system(command = "Rscript ~/icgc/jmmlc_pbat/data/Bock_Celltype_Predictor/cell_type_predictor_code/cellTypePredictor.R \\
--features ~/icgc/jmmlc_pbat_normals/data/odcf_md/analysis/Celltype_Prediction_BockTrain_JMML_All/meth.tsv \\
--annot ~/icgc/jmmlc_pbat_normals/data/odcf_md/analysis/Celltype_Prediction_BockTrain_JMML_All/sample_anno.tsv \\
--out ~/icgc/jmmlc_pbat_normals/data/odcf_md/analysis/Celltype_Prediction_BockTrain_JMML_All/predictCelltype_TrainBock_ClassJMML/"
)




#with Bock normal reference alll and JmmL
analysis.dir <- "~/icgc/jmmlc_pbat_normals/data/odcf_md/analysis/Celltype_Prediction_BockTrain_withOrigin_JMML"
dir.create(analysis.dir)
#read cell type region tation to extract methylation values
input_regions <- fread("~/icgc/jmmlc_pbat/data/cell_type_predictor_input/regionAnnot.tsv")
input_regions <- makeGRangesFromDataFrame(as.data.frame(input_regions))
#extract methylation levels
meth_predictorInput_jmml <- bsseq::getMeth(bsseq_all, regions= input_regions, type = "raw", what=c("perRegion"))

#read bock normal references
meth_bock <- fread("~/icgc/jmmlc_pbat/data/Bock_Celltype_Predictor/cell_type_predictor_input/methMatrix.tsv")
meth_bock<- as.data.frame(meth_bock)
annno_bock <- fread("~/icgc/jmmlc_pbat/data/Bock_Celltype_Predictor/cell_type_predictor_input/sampleAnnot.tsv")
annno_bock<- as.data.frame(annno_bock)
annno_bock$class <- paste0(annno_bock$class, "_",annno_bock$cellSourceCurated )
#subset only train data
annno_bock_sub <- annno_bock
#subset classes with less than 5 samples
temp <- names(table(annno_bock_sub$class)[table(annno_bock_sub$class)> 5])
annno_bock_sub <- annno_bock_sub[annno_bock_sub$class %in% temp,]
#make all train data
annno_bock_sub$train <- TRUE
meth_bock_sub <- meth_bock[, colnames(meth_bock) %in% as.character(annno_bock_sub$sampleName)]
dim(annno_bock)
dim(annno_bock_sub)

#get sample tation with matching columns 
pheno <- pData(bsseq_all)
pheno_new <- data.frame(sampleName=rownames(pheno), class=paste0(pheno$Epigenotype), 
donorId=pheno$Patient, nominalCellNumber=100, cellSourceCurated="JMML",
    train=ifelse(pheno$Sample_Type=="normal", TRUE, FALSE)
    )
    
pheno_new_tumor <- pheno_new[pheno_new$train==FALSE,]
meth_predictorInput_jmml_tumor <- meth_predictorInput_jmml[,as.character(pheno_new_tumor$sampleName)]

#merge bock with own data
meth_predictorInput <- cbind(meth_bock_sub,meth_predictorInput_jmml_tumor )
pheno_predictorInput <- rbind(annno_bock_sub[, colnames(pheno_new_tumor)], pheno_new_tumor)
####
table(colnames(meth_predictorInput) %in% pheno_predictorInput$sampleName)

 write.table(pheno_predictorInput, file.path(analysis.dir, "sample_anno.tsv"), 
    quote=FALSE, sep="\t", row.names=FALSE)
    
write.table(meth_predictorInput, file.path(analysis.dir, "meth.tsv"), 
    quote=FALSE, sep="\t", row.names=FALSE)
  
#run commandline cell type tation
# rm -r ~/icgc/jmmlc_pbat_normals/data/odcf_md/analysis/Celltype_Prediction_BockTrain_withOrigin_JMML/predictCelltype_TrainBock_ClassJMML/
system(command = "Rscript ~/icgc/jmmlc_pbat/data/Bock_Celltype_Predictor/cell_type_predictor_code/cellTypePredictor.R \\
    --features ~/icgc/jmmlc_pbat_normals/data/odcf_md/analysis/Celltype_Prediction_BockTrain_withOrigin_JMML/meth.tsv \\
    --annot ~/icgc/jmmlc_pbat_normals/data/odcf_md/analysis/Celltype_Prediction_BockTrain_withOrigin_JMML/sample_anno.tsv \\
    --out ~/icgc/jmmlc_pbat_normals/data/odcf_md/analysis/Celltype_Prediction_BockTrain_withOrigin_JMML/predictCelltype_TrainBock_ClassJMML/"




#with Bock normal reference alll and JmmL and our normal references
analysis.dir <- "~/icgc/jmmlc_pbat_normals/data/odcf_md/analysis/Celltype_Prediction_BockTrain_withOrigin_JMML_ourNormal"
dir.create(analysis.dir)
#read cell type region tation to extract methylation values
input_regions <- fread("~/icgc/jmmlc_pbat/data/cell_type_predictor_input/regionAnnot.tsv")
input_regions <- makeGRangesFromDataFrame(as.data.frame(input_regions))
#extract methylation levels
meth_predictorInput_jmml <- bsseq::getMeth(bsseq_all, regions= input_regions, type = "raw", what=c("perRegion"))

#read bock normal references
meth_bock <- fread("~/icgc/jmmlc_pbat/data/Bock_Celltype_Predictor/cell_type_predictor_input/methMatrix.tsv")
meth_bock<- as.data.frame(meth_bock)
annno_bock <- fread("~/icgc/jmmlc_pbat/data/Bock_Celltype_Predictor/cell_type_predictor_input/sampleAnnot.tsv")
annno_bock<- as.data.frame(annno_bock)
annno_bock$class <- paste0(annno_bock$class, "_",annno_bock$cellSourceCurated )
#subset only train data
annno_bock_sub <- annno_bock
#subset classes with less than 5 samples
temp <- names(table(annno_bock_sub$class)[table(annno_bock_sub$class)> 5])
annno_bock_sub <- annno_bock_sub[annno_bock_sub$class %in% temp,]
#make all train data
annno_bock_sub$train <- TRUE
meth_bock_sub <- meth_bock[, colnames(meth_bock) %in% as.character(annno_bock_sub$sampleName)]
dim(annno_bock)
dim(annno_bock_sub)


#get sample tation with matching columns 
pheno <- pData(bsseq_all)
pheno_new <- data.frame(sampleName=rownames(pheno), class=ifelse(pheno$Sample_Type=="normal", paste0("C010_",as.character(pheno$Celltype)), as.character(pheno$Epigenotype)), donorId=pheno$Patient, nominalCellNumber=100, cellSourceCurated="JMML",
    train= FALSE
    )
    
pheno_new_tumor <- pheno_new
meth_predictorInput_jmml_tumor <- meth_predictorInput_jmml[,as.character(pheno_new_tumor$sampleName)]

#merge bock with own data
meth_predictorInput <- cbind(meth_bock_sub,meth_predictorInput_jmml_tumor )
pheno_predictorInput <- rbind(annno_bock_sub[, colnames(pheno_new_tumor)], pheno_new_tumor)

table(colnames(meth_predictorInput) %in% pheno_predictorInput$sampleName)

 write.table(pheno_predictorInput, file.path(analysis.dir, "sample_anno.tsv"), 
    quote=FALSE, sep="\t", row.names=FALSE)
    
write.table(meth_predictorInput, file.path(analysis.dir, "meth.tsv"), 
    quote=FALSE, sep="\t", row.names=FALSE)
  
#run commandline cell type tation
# rm -r ~/icgc/jmmlc_pbat_normals/data/odcf_md/analysis/Celltype_Prediction_BockTrain_withOrigin_JMML_ourNormal/predictCelltype_TrainBock_ClassJMML/
system(command = "Rscript ~/icgc/jmmlc_pbat/data/Bock_Celltype_Predictor/cell_type_predictor_code/cellTypePredictor.R \\
    --features ~/icgc/jmmlc_pbat_normals/data/odcf_md/analysis/Celltype_Prediction_BockTrain_withOrigin_JMML_ourNormal/meth.tsv \\
    --annot ~/icgc/jmmlc_pbat_normals/data/odcf_md/analysis/Celltype_Prediction_BockTrain_withOrigin_JMML_ourNormal/sample_anno.tsv \\
    --out ~/icgc/jmmlc_pbat_normals/data/odcf_md/analysis/Celltype_Prediction_BockTrain_withOrigin_JMML_ourNormal/predictCelltype_TrainBock_ClassJMML/"

Rscript ~/icgc/jmmlc_pbat/data/Bock_Celltype_Predictor/cell_type_predictor_code/cellTypePredictor.R \
--features ~/icgc/jmmlc_pbat_normals/data/odcf_md/analysis/Celltype_Prediction_BockTrain_withOrigin_JMML_ourNormal/meth.tsv \
--annot ~/icgc/jmmlc_pbat_normals/data/odcf_md/analysis/Celltype_Prediction_BockTrain_withOrigin_JMML_ourNormal/sample_anno.tsv \
--out ~/icgc/jmmlc_pbat_normals/data/odcf_md/analysis/Celltype_Prediction_BockTrain_withOrigin_JMML_ourNormal/predictCelltype_TrainBock_ClassJMML/



