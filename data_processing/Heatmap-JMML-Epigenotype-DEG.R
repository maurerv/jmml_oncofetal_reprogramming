library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(plyr)
library(ggplot2)
ht_opt$message = FALSE
col = circlize::colorRamp2(c(-3,0,3),
                                 c("#2166ac","white", "#b2182b"))
df.sc <- readRDS('JMML_HSC-Epigenotype.DEG.Avg.Scaled.rds')
gene.ann <- readRDS('JMML_HSC-Epigenotype.DEG.Annotation.rds')
jmml.color <- readRDS('JMML_HSC-Epigenotype.DEG.Color.rds')
jmml.ann <- readRDS('JMML_HSC-Epigenotype.Annotation.rds')

genecolor <- c('#b2182b','#2166ac')
names(genecolor) <- levels(gene.ann$Gene)
genecolor <- list(Gene = genecolor)

mg <- read.table('Mark.HMnonHM.txt')
ind <- match(mg,rownames(df.sc))
ha = rowAnnotation(foo = anno_mark(at = ind, labels = as.character(mg$V1)))

Heatmap(df.sc,
             col = col,
             right_annotation = ha,
             left_annotation = rowAnnotation(df = gene.ann,col=genecolor),
             top_annotation = HeatmapAnnotation(df = jmml.ann, 
                                                col = jmml.color),
             show_column_names = T,show_row_names = F,
             cluster_rows = F,cluster_columns = F)