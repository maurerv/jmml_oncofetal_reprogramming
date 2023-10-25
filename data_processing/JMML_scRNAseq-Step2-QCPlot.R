library(ggpubr)
meta <- readRDS("JMML-REF.MetaData.rds")
m1 <- meta[meta$Project == 'JMML',]
m1$Group <- paste(m1$Project,m1$Cell,sep = '_')
pdf('QCPlot.nCount_RNA.JMML.pdf',7,7)
ma <- m1[m1$Group == 'JMML_CD34',]
gghistogram(ma,'nCount_RNA',facet.by = 'Pseudonym',bins = 50)+
  ggtitle('JMML_CD34')
ma <- m1[m1$Group == 'JMML_MNC',]
gghistogram(ma,'nCount_RNA',facet.by = 'Pseudonym',bins = 50)+
  ggtitle('JMML_MNC')
dev.off()
pdf('QCPlot.nFeature_RNA.JMML.pdf',7,7)
ma <- m1[m1$Group == 'JMML_CD34',]
gghistogram(ma,'nFeature_RNA',facet.by = 'Pseudonym',bins = 50)+
  ggtitle('JMML_CD34')
ma <- m1[m1$Group == 'JMML_MNC',]
gghistogram(ma,'nFeature_RNA',facet.by = 'Pseudonym',bins = 50)+
  ggtitle('JMML_MNC')
dev.off()

