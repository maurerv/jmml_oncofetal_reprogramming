source('LogisticRegression.Core.R')

matrix <- readRDS('HealthyReference.Matrix.rds')
celltypes <- readRDS('HealthyReference.Celltype-Annotation.rds')
jmml_matrix <- readRDS('JMML.Matrix.rds')

REF.TrainModel = trainModel(matrix, as.character(celltypes), nParallel = 24)
print("===END trainModel===")

JMML.PredictSimilarity = predictSimilarity(REF.TrainModel, jmml_matrix, logits=FALSE, minGeneMatch = 0.5)
print("===END trainModel===")