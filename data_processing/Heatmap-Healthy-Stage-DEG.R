library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(plyr)
library(ggplot2)
ht_opt$message = FALSE


cpl_patient_col_3 <- c(FET = "#1f0b0d", NEO = "#79675c", 
                       JUV = "#6f8987", ADU ="#c3c6d8", 
                       P1 = "#0058b4", P2 = "#2188c9", 
                       P3 = "#fbbb25", P4 = "#fca349", 
                       P5 = "#ff6b36", P6 = "#e34e2e", 
                       P7 = "#c33126", P8 = "#a41220")
samplecolorlist <- list(Patient = cpl_patient_col_3)

df.sc <- readRDS("JMML-HSC_REF-HSC.DEG.Avg.rds")

sample.ann <- data.frame(Patient = colnames(df.sc))
rownames(sample.ann) <- colnames(df.sc)

gene.ann <- readRDS("JMML-HSC_REF-HSC.DEG.Annotation.rds")

genecolor <- cpl_patient_col_3[1:4]
genecolor <- list(Gene = genecolor)

gl <- read.csv('JMML-HSC_REF-HSC.DEG.Label.csv')
ind <- match(gl$gene,rownames(df.sc))

df.sc <- df.sc[,names(cpl_patient_col_3)]
df.sc <- df.sc[,c(5:12,1:4)]
df.sc.oo <- df.sc

df.sc <- df.sc.oo
for (i in 1:nrow(df.sc)) {
  a <- max(df.sc[i,9:12])
  b <- min(df.sc[i,9:12])
  df.sc[i,][df.sc[i,] > a] <- a
  df.sc[i,][df.sc[i,] < b] <- b
}
df.sc <- t(scale(t(df.sc)))
df.sc[df.sc > 3] <- 3
df.sc[df.sc < -3] <- -3

scico_map = circlize::colorRamp2(c(-4,-3,0,3,4),
                                 c('#053061',"#2166ac","white", "#b2182b",'#67001f'))

for (i in 1:nrow(df.sc)) {
  a <- max(df.sc.oo[i,9:12])
  b <- min(df.sc.oo[i,9:12])
  df.sc[i,][df.sc.oo[i,] > a] <- 4
  df.sc[i,][df.sc.oo[i,] < b] <- -4
}

Heatmap(df.sc,
             col = scico_map,
             row_split = gene.ann$Gene,
             column_split = factor(c(rep('JMML',8),
                                     c('FET', 'NEO', 'JUV', 'ADU')),
                                   levels = c('JMML','FET', 'NEO', 'JUV', 'ADU')),
             cluster_column_slices=F,column_title = NULL,
             right_annotation = rowAnnotation(foo = anno_mark(at = ind, labels = inda)),
             left_annotation = rowAnnotation(df = gene.ann,col=genecolor),
             top_annotation = HeatmapAnnotation(df = sample.ann, 
                                                col = samplecolorlist),
             show_column_names = T,show_row_names = F,
             cluster_row_slices=F,
             cluster_rows = T,cluster_columns = T)

