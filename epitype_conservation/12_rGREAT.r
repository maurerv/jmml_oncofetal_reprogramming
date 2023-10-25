
#Libraries
library(rGREAT)
library(forcats)
library(ggpubr)
#Directories
input.dir <- "/omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md/analysis/"
analysis.dir <-  "/omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md/analysis/DMR_sub"
analysis.dir <-  "/omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md/analysis/200519_DMR_HM_vs_LM"
data.dir <- "/C010-datasets/Internal/mitotic_clock/data/"

#load data
datasets.dir <- "c010-datasets/Internal/COPD/enrichment_databases/"
output.dir <- "/omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md/analysis/bsseq"
analysis.dir <-  "/omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md/analysis/200612_DMR_model_sub_repMerged"

#load data
bsseq_all <- readRDS(file.path(output.dir , "bsseq_all_snpfil_sub_cov_repMerged.rds"))
dmrs_final<- readRDS(file.path(analysis.dir, "dmrs_gr_sub_MethDiff.rds"))
dmrs_red<- readRDS(file.path(analysis.dir, "dmrs_gr_sub_MethDiff_anno_reduced.rds"))

dmrs_final <- lapply(dmrs_final, function(x){
    x$direction = ifelse(x$diff.Methy>0, "hyper", "hypo")
    x
})

#plotting function
bubblPlotEnr<- function(data, n, set){
enr <- data[order(data$Hyper_Adjp_BH, decreasing=FALSE),]
datasetsub<-enr[1:n,]
datasetsub$LogP<- abs(log10(datasetsub$Hyper_Adjp_BH))
datasetsub$order <- 1:nrow(datasetsub)
datasetsub$significant<- ifelse(datasetsub$LogP< -log10(0.05), "No", "Yes" )
datasetsub$set <- set
datasetsub<- mutate(datasetsub, name = fct_reorder(name, LogP))

label_func <- function(x){
    breaks <- x
    breaks[breaks>=300] <- ">=300"
    breaks
}
ggplot(data = datasetsub, aes(y=name, x=set))+coord_fixed()+
    geom_point(aes(size=LogP, fill=Hyper_Fold_Enrichment), pch=21)+
    scale_fill_gradient2( midpoint = 0, low="white", high="darkred", name = "Hyper fold enrichment")+
    scale_size(name="Q-value\n(-log10)", labels = label_func) +
    theme(text =element_text(size=14, color="black", family = "sans"),
          axis.ticks = element_blank(), axis.line = element_blank(), 
          axis.text.x=element_text(size=12, angle = 90, vjust = 0, color="black", family="sans"),
          axis.text.y=element_text(size=12, family="sans", color="black"))+
    scale_x_discrete(name=NULL)+
    theme(legend.text=element_text(size=12, family="sans"), 
          legend.title=element_text(size=12, family= "sans"),
          legend.background = element_rect(fill="white", color="white"),
          panel.background =  element_rect(fill="white"), panel.grid.major = element_line(color="lightgrey"),
          legend.key = element_rect(fill="white"))+rremove("ylab")
}



bubblPlotEnr_comb <- function(data, n){
enr <- data[order(data$Hyper_Adjp_BH, decreasing=FALSE),]
name_OI <- enr[1:n,]$name
datasetsub<-enr[enr$name %in% name_OI,]
datasetsub$LogP<- abs(log10(datasetsub$Hyper_Adjp_BH))
datasetsub$order <- 1:nrow(datasetsub)
datasetsub$significant<- ifelse(datasetsub$LogP< -log10(0.05), "No", "Yes" )
datasetsub<- mutate(datasetsub, name = fct_reorder(name, LogP))

label_func <- function(x){
    breaks <- x
    breaks[breaks>=300] <- ">=300"
    breaks
}
ggplot(data = datasetsub, aes(y=name, x=set))+coord_fixed()+
    geom_point(aes(size=LogP, fill=Hyper_Fold_Enrichment), pch=21)+
    scale_fill_gradient2( midpoint = 0, low="white", high="darkred", name = "Hyper fold enrichment")+
    scale_size(name="Q-value\n(-log10)", labels = label_func) +
    theme(text =element_text(size=14, color="black", family = "sans"),
          axis.ticks = element_blank(), axis.line = element_blank(), 
          axis.text.x=element_text(size=12, angle = 90, vjust = 0, color="black", family="sans"),
          axis.text.y=element_text(size=12, family="sans", color="black"))+
    scale_x_discrete(name=NULL)+
    theme(legend.text=element_text(size=12, family="sans"), 
          legend.title=element_text(size=12, family= "sans"),
          legend.background = element_rect(fill="white", color="white"),
          panel.background =  element_rect(fill="white"), panel.grid.major = element_line(color="lightgrey"),
          legend.key = element_rect(fill="white"))+rremove("ylab")
}

#rgreat enrichment
#rules for rGreat
Dist<- 10000
rule<-"basalPlusExt"


for (i in names(dmrs_final)){
    dir.create(file.path(analysis.dir, "rGreat"))
job_all <- submitGreatJob(dmrs_final[[i]], species= "hg19", includeCuratedRegDoms = TRUE,rule= c("basalPlusExt"),adv_upstream= 5.0,adv_downstream= 1.0,adv_span= 1000.0,adv_twoDistance= 1000.0,
adv_oneDistance= 1000.0, request_interval = 300, max_tries = 10, version = "default")
job_hypo <- submitGreatJob(dmrs_final[[i]][dmrs_final[[i]]$direction=="hypo",], species= "hg19", includeCuratedRegDoms = TRUE,rule= c("basalPlusExt"),adv_upstream= 5.0,adv_downstream= 1.0,adv_span= 1000.0,adv_twoDistance= 1000.0,
adv_oneDistance= 1000.0, request_interval = 300, max_tries = 10, version = "default")
job_hyper <- submitGreatJob(dmrs_final[[i]][dmrs_final[[i]]$direction=="hyper",],  species= "hg19", includeCuratedRegDoms = TRUE,rule= c("basalPlusExt"),adv_upstream= 5.0,adv_downstream= 1.0,adv_span= 1000.0,adv_twoDistance= 1000.0,
adv_oneDistance= 1000.0, request_interval = 300, max_tries = 10, version = "default")

#all DMRs enrichments
enrichments <- list()
dir.create(file.path(analysis.dir, "rGreat", "allDMRs"))
for(j in availableCategories(job_all)){
    enrichments[[j]] <- getEnrichmentTables(job_all, category=j)
    print(j)
    for (k in names(enrichments[[j]])){
        table <- enrichments[[j]][[k]]
    
        pdf(file.path(analysis.dir, "rGreat","allDMRs",paste0("allDMRs_bubble","_",j,"_",k, ".pdf")), height=7, width=14)
        print(bubblPlotEnr(table, 30, "all DMRs"))
        dev.off()
        write.table(table, file.path(analysis.dir, "rGreat","allDMRs",paste0("allDMRs","_",j,"_",k, ".txt")))
    }
}
#hypo DMRs enrichments
enrichments <- list()
dir.create(file.path(analysis.dir, "rGreat", "hypoDMRs"))
for(j in availableCategories(job_hypo)){
    enrichments[[j]] <- getEnrichmentTables(job_hypo, category=j)
    print(j)
    for (k in names(enrichments[[j]])){
        table <- enrichments[[j]][[k]]
    
        pdf(file.path(analysis.dir, "rGreat","hypoDMRs",paste0("hypoDMRs_bubble","_",j,"_",k, ".pdf")), height=7, width=14)
        print(bubblPlotEnr(table, 30, "hypo DMRs"))
        dev.off()
        write.table(table, file.path(analysis.dir, "rGreat","hypoDMRs",paste0("hypoDMRs","_",j,"_",k, ".txt")))
    }
}

#hyper DMRs enrichments
enrichments <- list()
dir.create(file.path(analysis.dir, "rGreat",  "hyperDMRs"))
for(j in availableCategories(job_hyper)){
    enrichments[[j]] <- getEnrichmentTables(job_hyper, category=j)
    print(j)
    for (k in names(enrichments[[j]])){
        table <- enrichments[[j]][[k]]
    
        pdf(file.path(analysis.dir, "rGreat", "hyperDMRs",paste0("hyperDMRs_bubble","_",j,"_",k, ".pdf")), height=7, width=14)
        print(bubblPlotEnr(table, 30, "hyper DMRs"))
        dev.off()
        write.table(table, file.path(analysis.dir, "rGreat", "hyperDMRs",paste0("hyperDMRs","_",j,"_",k, ".txt")))
    }
}

#hyper hypo combined DMRs enrichments
enrichments_hyper <- list()
enrichments_hypo <- list()
dir.create(file.path(analysis.dir, "rGreat",  "hyperhypoDMRs"))
for(j in availableCategories(job_hyper)){
    enrichments_hyper[[j]] <- getEnrichmentTables(job_hyper, category=j)
    enrichments_hypo[[j]] <- getEnrichmentTables(job_hypo, category=j)

    print(j)
    for (k in names(enrichments_hyper[[j]])){
        table_hyper <- enrichments_hyper[[j]][[k]]
        table_hyper$set <- "hyper"
        table_hypo <- enrichments_hypo[[j]][[k]]
        table_hypo$set <- "hypo"
        table <- rbind(table_hyper,table_hypo)

        pdf(file.path(analysis.dir, "rGreat", "hyperhypoDMRs",paste0("hyperhypoDMRs_bubble","_",j,"_",k, ".pdf")), height=7, width=14)
        print(bubblPlotEnr_comb(table, 30))
        dev.off()
        write.table(table, file.path(analysis.dir, "rGreat", "hyperhypoDMRs",paste0("hyperhypoDMRs","_",j,"_",k, ".txt")))
    }
}
}
