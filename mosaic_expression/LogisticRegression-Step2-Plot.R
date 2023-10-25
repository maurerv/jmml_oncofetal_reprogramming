library(ComplexHeatmap)
library(circlize)
df <- readRDS('JMML.PredictSimilarity.rds')
meta <- readRDS('JMML.Meta.rds')

# Set colors
samplecolor <- c(P1 = "#0058b4", P2 = "#2188c9", 
                       P3 = "#fbbb25", P4 = "#fca349", 
                       P5 = "#ff6b36", P6 = "#e34e2e", 
                       P7 = "#c33126", P8 = "#a41220")
samplecolorlist <- list(Sample = samplecolor)
col_fun <- colorRamp2(c(0, 0.5, 1), c('#e7e1ef', '#ffffff', '#dd1c77'))

# Getting prediction scores of JMML_CD34 
index.CD34 <- meta$Cell == 'CD34'
meta.CD34 <- meta[index.CD34,]
predict.CD34 <- df[index.CD34, ]

# Averaging scores by donor
classes <- meta.CD34$Donor
predict.CD34.Donor <-
  apply(predict.CD34, 2, function(e)
    sapply(split(e, classes), mean))

# Annotation for heatmap
sample.ann <- data.frame(Sample = rownames(predict.CD34.Donor))
rownames(sample.ann) <- rownames(predict.CD34.Donor)
predict.CD34.Donor <- t(predict.CD34.Donor)

# Cleaning NA scores
a <- rowSums(predict.CD34.Donor)
b <- is.na(a)
predict.CD34.Donor <- predict.CD34.Donor[!b, ]

# Ploting
Heatmap(
  predict.CD34.Donor,
  col = col_fun,
  top_annotation = HeatmapAnnotation(df = sample.ann, col = samplecolorlist),
  show_column_names = T,
  show_row_names = T,
  cluster_rows = T,
  cluster_columns = T
)
