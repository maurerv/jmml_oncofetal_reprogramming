source('LogisticRegression.Core.R')

matrix <- readRDS('/icgc/dkfzlsdf/analysis/OE0219_projects/JMMLC/scRNA_Result/LogisticRegression/REF.Matrix.rds')
celltypes <- readRDS('/icgc/dkfzlsdf/analysis/OE0219_projects/JMMLC/scRNA_Result/LogisticRegression/REF.Celltype.rds')
jmml_matrix <- readRDS('/icgc/dkfzlsdf/analysis/OE0219_projects/JMMLC/scRNA_Result/LogisticRegression/JMML.Matrix.rds')

gm = trainModel(matrix, as.character(celltypes), nParallel = 24)
saveRDS(gm,'REF.TrainModel.rds')
print("===END trainModel===")

pp = predictSimilarity(gm, jmml_matrix, logits=FALSE, minGeneMatch = 0.5)
saveRDS(pp,'JMML.PredictSimilarity.rds')
print("===END trainModel===")