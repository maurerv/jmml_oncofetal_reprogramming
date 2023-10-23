################################################
#
# 07_cluster_contribution
# Author: Maximilian Schönung
# Date: 23.10.2023
#
################################################

# Set the paths -----------------------------------------------------------
setwd("/omics/odcf/analysis/OE0565_projects/ptpn11/")
data.dir <- "/omics/odcf/analysis/OE0565_projects/ptpn11/data/"
plot.dir <- "/omics/odcf/analysis/OE0565_projects/ptpn11/plots/"

subDir <- "02_functional/"
if (!file.exists(paste0(plot.dir,subDir))){
  dir.create(paste0(plot.dir,subDir),recursive = T)
}

# Load the libraries ------------------------------------------------------
library(Seurat)
library(dplyr)
library(ggplot2)

# Cluster colour ----------------------------------------------------------
col.vec = c("HSC"="black","MPP1"="grey89","MPP2"="grey50",
            "MPP_Cycling1"="ivory3","MPP_Cycling2"="ivory4",
            "MEP"="coral2","Ery1"="orangered2","Ery2"="orangered3","MPP_EryMeg"="darkred",
            "MyP1"="darkolivegreen4","MyP2"="darkolivegreen3","MonoP"="darkolivegreen2",
            "NeutroP"="palegreen3",
            "Mono1"="goldenrod","Mono2"="gold3","DC"="gold",
            "Neutro1"="orange","Neutro2"="darkorange","Neutro3"="darkorange1","Neutro4"="darkorange2","Neutro5"="darkorange3",
            "MPP_Lympho"="lightsteelblue4","LyP"="steelblue1","T-cell"="lightskyblue",
            "B-cell1"="paleturquoise","B-cell2"="paleturquoise1","B-cell3"="paleturquoise2","B-cell4"="paleturquoise3",
            "Basophils"="blueviolet")

# Load the 10x Data -------------------------------------------------------
seurat <- readRDS(paste0(data.dir,"2023-09-18_bm_ptpn11_sct_annotated.RDS"))
all(names(col.vec)%in%unique(Idents(seurat)))
all(unique(Idents(seurat))%in%names(col.vec))


# Circular foldChange Plot ------------------------------------------------
meta2 <- seurat@meta.data
meta2$anno <- as.character(seurat@active.ident)
table(meta2$leuk)
set.seed(12345)
ind1 <- sample(which(meta2$leuk==1),10000)
ind2 <- sample(which(meta2$leuk==0),10000)
meta2 <- meta2[c(ind1,ind2),]
df2 <- as.data.frame.matrix(table(meta2$orig.ident,meta2$anno))
df.t2 <- t(df2)
df.leuk2 <- data.frame("E76K_BM"=rowSums(df.t2[,c("345BM","347BM","349BM")]),
                       "WT_BM"=rowSums(df.t2[,c("348BM","351BM","352BM")]))
df.leuk2$fc_bm <- log2(df.leuk2[,1]/df.leuk2[,2])
df.leuk2$celltype <- factor(rownames(df.leuk2),
                            levels=c("HSC","MPP1","MPP2",
                                     "MPP_Cycling1","MPP_Cycling2",
                                     "MPP_EryMeg",
                                     "MyP1","MyP2","MonoP","NeutroP",
                                     "Mono1","Mono2","DC",
                                     "Neutro1","Neutro2","Neutro3","Neutro4","Neutro5",
                                     "MEP","Ery1","Ery2",
                                     "MPP_Lympho","LyP","T-cell",
                                     "B-cell1","B-cell2","B-cell3","B-cell4",
                                     "Basophils"))

list.test <- list()
for(i in names(table(meta2$anno))){
  list.test[[i]] <- fisher.test(table(meta2$anno==i,meta2$leuk))$p.value
}
p.df <- data.frame(do.call(rbind,list.test))
p.df$adj <- p.adjust(p.df[,1],method = "BH")
p.df$sig <- p.df$adj<0.05
p.df$round <- round(p.df$adj,4)
write.table(p.df,"tables/2023-09-18_cluster_freq_pvalue.txt",sep="\t",quote=F,row.names=T)

all(df.leuk2$celltype==rownames(p.df))
df.leuk2$p.adj <- p.df$round
df.leuk2$p_lab <- ifelse(df.leuk2$p.adj==0, "p<0.0001", paste0("p=",format(df.leuk2$p.adj,scientific = F)))
df.leuk3 <- df.leuk2
df.leuk3$fac <- as.numeric(df.leuk3$celltype)
angle = (max(df.leuk3$fac) + 0.45 - df.leuk3$fac) / max(df.leuk3$fac) * 2*pi * 180/pi + 90
angle = angle %% 360
flip = angle > 90 & angle < 270
angle = ifelse(flip, angle + 180, angle)
hjust.val <- ifelse(flip, 1, 0)



pdf(paste0(plot.dir,subDir,Sys.Date(),"_cluster_contribution_circular_10k_pval.pdf"))
ggplot(df.leuk2,aes(celltype,fc_bm,group=1,fill=celltype))+
  geom_hline(yintercept = 1,linetype="dashed")+
  geom_hline(yintercept = -1,linetype="dashed")+
  geom_hline(yintercept = -2,linetype="dashed")+
  geom_col()+
  geom_hline(yintercept = 0)+
  coord_polar()+
  ylim(-5,3)+
  theme_minimal()+
  scale_fill_manual(values = col.vec)+
  ggtitle("log2FC(E76K/WT)")+
  geom_text(aes(label = celltype, x = celltype, y = 1.8, angle = angle, vjust = sign(angle) * .5,hjust=hjust.val),fontface=ifelse(df.leuk2$p.adj<0.01, 2, 1)) +
  geom_text(aes(label = p_lab, x = celltype, y = 1.8, angle = angle, vjust = sign(angle) * 3,hjust=hjust.val),size=2,fontface=ifelse(df.leuk2$p.adj<0.01, 2, 1)) +
  theme(
    legend.position="none",
    panel.background = element_rect(fill = "white", color = "white"),
    panel.grid = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank()
  )
dev.off()

