library(here)
# Directories -------------------------------------------------------------
maindir= here()
resdir <- file.path(maindir,"output","scRNAseq","PDX_BC976")
QCdir <- file.path(maindir,"output","scRNAseq","QC")

output <- file.path(resdir, "Unsupervised", "TFmotif") ; if(!file.exists(output)){dir.create(output)}

resdir <- file.path(resdir, "Unsupervised") ; if(!file.exists(resdir)){dir.create(resdir)}
resdir_UMAP <- file.path(resdir,"UMAP") ; if(!file.exists(resdir_UMAP)){dir.create(resdir_UMAP)}
resdir_boxplots <- file.path(resdir,"boxplots") ; if(!file.exists(resdir_boxplots)){dir.create(resdir_boxplots)}
resdir_heatmaps = file.path(resdir,"Heatmaps"); if(!dir.exists(resdir_heatmaps)) dir.create(resdir_heatmaps)
resdir_sub = file.path(resdir,"sub500"); if(!dir.exists(resdir_sub)) dir.create(resdir_sub)

RDatadir <- file.path(resdir,"RData") ; if(!file.exists(RDatadir)){dir.create(RDatadir)}

source(file.path(maindir,"Scripts","global_var.R"))

annotCol <- c("sample_id","total_features","rRNA","louvain_partition","cons_BC_lenti","CDH2","TWIST1","TGFB1")
# Select initial population, 4 'persister' states (early) and 3 'resistant' states (late)
control <- c("#E0E0E0","#BDBDBD","#757575","#ffcbcb")
persister_color <- c("#DCEDC8","#9CCC65","#4CAF50","#009688")
res_color <- c("#FFEB3B","#FFC107","#FF9800","#FF5722")
color_PDX <- c(control[c(1)],persister_color[c(1)])


#Common palette for all MM468 experiments
control <- c("#E0E0E0","#BDBDBD","#757575","#ffcbcb")
persister_color <- c("#DCEDC8","#9CCC65","#4CAF50","#009688")
res_color <- c("#FFEB3B","#FFC107","#FF9800","#FF5722")
color_PDX <- c(control[c(1,2)],persister_color[c(1,3)],res_color[c(2)])


samples_PDX = c("HBCx95_chemonaive", "HBCx95_recurrent", "HBCx95_residual",
                "HBCx95_persister",   "HBCx95_resistant")
corres <- data.frame(PDX=samples_PDX,color=color_PDX)

# If re-computing from scratch - set RECOMPUTE to TRUE else if you just want to
# plot the UMAPS, set RECOMPUTE to FALSE 
RECOMPUTE = FALSE

# library(tftargets)


# Overexpressed genes
load(file.path(maindir,"output","scRNAseq","PDX_BC976","Unsupervised","RData","LogCounts.RData"))
LogCounts = LogCounts[,grep("persister",colnames(LogCounts))]
dim(LogCounts)
BinCounts = LogCounts
BinCounts[BinCounts>0] = 1
rowSums(BinCounts)

persister_DA = readxl::read_xlsx(file.path(maindir,"output","scRNAseq","PDX_BC976","Supervised",
                                              "Tables","Differential_analysis_Limma_logFC_Pers_vs_UNT_1.58.xlsx"))
persister_DA$percent_expressed = rowSums(BinCounts[match(persister_DA$Symbol,rownames(BinCounts)),])/ncol(BinCounts)
persister_genes = persister_DA$Symbol
hist(persister_DA$percent_expressed, breaks = 50)
other_genes = readxl::read_xlsx(file.path(maindir,"output","scRNAseq","PDX_BC976","Supervised",
                                          "Tables","Differential_analysis_Limma_logFC_Pers_vs_UNT_1.58.xlsx"), sheet = 3)$Symbol
other_genes = other_genes[which(!other_genes %in% persister_genes)]


# Run ChEA3 enrichment
library(httr)
library(curl)
library(jsonlite)

genes = persister_genes

set.seed(47)
system.time({
    list_TF_enrichment = enrich_for_TF_ChEA3(persister_genes, n_random = 1000, all_genes = c(persister_genes,other_genes))
})



save(list_TF_enrichment, file = file.path(output,"list_TF_enrichment.RData"))
load(file = file.path(output,"list_TF_enrichment.RData"))

scores = as.data.frame(sapply(list_TF_enrichment, function(x) x$Score[1:20]))
scores = scores %>% tidyr::gather("Run", "Score") 
scores = scores %>% mutate("Class" =  gsub("random_genes_.*","random_genes",Run))
scores$Score = 1/as.numeric(scores$Score)

png(file.path(output, "TF_ChEA3_enrichment_score.png"), height = 2000, width = 2000, res = 300)
p  = scores %>% ggplot(aes(x = Class, y = Score)) + ylab("1 / ChEA3 Mean Rank Score ") + 
    geom_violin(aes(fill = Class)) + ylim(c(quantile(scores$Score,0.01),quantile(scores$Score,0.99))) +
    ggpubr::stat_compare_means(paired = F, method = "t.test", label.x.npc = 0.4) +
    theme_classic()
print(ggpubr::add_summary(p, fun = "mean_sd"))
dev.off()

scores_per_TF = sapply(list_TF_enrichment[2:101], function(x) {
    vec = x$Score
    names(vec) = x$TF
    as.numeric(vec[order(names(vec))])
    })
rownames(scores_per_TF) = sort(list_TF_enrichment[[1]]$TF)

# Test whether the inverse score of the top 100 TFs is greater than the score of the top 100 TF enriched in 
# random gene sets
t.test(1/as.numeric(list_TF_enrichment$Genes_of_interest$Score[1:100]),1/rowMeans(scores_per_TF), alternative = "greater")

scores_per_TF_df = data.frame(TF = sort(list_TF_enrichment[[1]]$TF),
                              mean_score =  rowMeans(scores_per_TF),
                              median_score = rowMedians(scores_per_TF))

i = 0
pdf(file.path(output, "TF_ChEA3_enrichment_score_TF_by_TF.pdf"))
for(TF in list_TF_enrichment$Genes_of_interest$TF[1:200]){
    i = i +1
    score_TF_in_persister = as.numeric(list_TF_enrichment$Genes_of_interest$Score[i])
    
    distrib = data.frame(MeanRank = scores_per_TF[TF,])
    n_random_below = length(which(distrib$MeanRank < score_TF_in_persister))
    p = distrib %>% ggplot(aes(MeanRank)) + xlab("Mean Rank (out of 1632 TFs)") + 
        geom_density(bw=50, color = "#009688", fill = alpha("#009688",0.3)) + theme_classic() + xlim(c(-100, 100 + max(c(score_TF_in_persister,distrib$MeanRank)))) +
        geom_vline(xintercept = score_TF_in_persister, color = "red", lty = 2) + 
        theme(axis.text = element_text(size = 12)) + 
        annotate("text", label = TF, x = score_TF_in_persister + 100,  y = 0.001, col = "red", lwd = 5) +
        ggtitle(paste0("Distribution of score for ", TF, " in 100x random gene sets.\nNumber of random subset with lower Mean Rank = ",n_random_below))
    print(p)
}
dev.off()


png(file.path(output, "top_TF_ChEA3_enrichment_by_score_in_random.png"), height = 2000, width = 2000, res = 300)
scores_per_TF_df %>% arrange(mean_score) %>% head(20) %>%
    mutate(TF = fct_reorder(TF, 1/mean_score, .desc = TRUE)) %>%
    ggplot(aes(x = TF, y = 1/mean_score)) + ylab("1 / ChEA3 average(Mean Rank Score)") + 
    geom_bar(fill = "#009688", stat="identity") +theme_classic() + 
    theme(axis.text.x = element_text(size = 15, angle=90),
          axis.text.y = element_text(size = 14)) + xlab("")
dev.off()

# Load ChEA3 enrichment
ChEA3_TF_enrichment_persister = list_TF_enrichment$Genes_of_interest
ChEA3_TF_enrichment_persister$n_targets_in_persisters = sapply(ChEA3_TF_enrichment_persister$Overlapping_Genes,
                                                               function(x) length(unlist(str_split(x, pattern = ","))))
ChEA3_TF_enrichment_persister$Score = as.numeric(ChEA3_TF_enrichment_persister$Score)

scRNA = readxl::read_xlsx(file.path(maindir, "output","scRNAseq","PDX_BC976","Supervised","Tables","Differential_analysis_Limma_logFC_Pers_vs_UNT_1.58.xlsx"))

pdf(file.path(output,"Percentage_of_targeted_persister_genes.pdf"))
ChEA3_TF_enrichment_persister %>% mutate(ratio_n_targets_in_persister = 100* n_targets_in_persisters / length(persister_genes)) %>%
    head(100) %>%
    mutate(TF = fct_reorder(TF, ratio_n_targets_in_persister, .desc = TRUE)) %>%
    ggplot(aes(x=TF, y= ratio_n_targets_in_persister)) + geom_bar(fill = "#009688", stat="identity") +
    theme_classic() + theme(axis.text.x = element_text(size = 15, angle=90),
                            axis.text.y = element_text(size = 14)) + xlab("") + ylim(c(0,50)) + ylab("% Targeted Persister Genes")
dev.off()

pdf(file.path(output,"Percentage_of_targeted_persister_genes_overexpressed.pdf"))
ChEA3_TF_enrichment_persister %>% mutate(ratio_n_targets_in_persister = 100* n_targets_in_persisters / length(persister_genes)) %>% head(100) %>% 
    filter(TF %in% scRNA$Symbol) %>%
    mutate(TF = fct_reorder(TF, ratio_n_targets_in_persister, .desc = TRUE)) %>%
    ggplot(aes(x=TF, y= ratio_n_targets_in_persister)) + geom_bar(fill = "#009688", stat="identity") +
    theme_classic() + theme(axis.text.x = element_text(size = 15, angle=90),
                            axis.text.y = element_text(size = 14)) + xlab("") + ylim(c(0,50)) + ylab("% Targeted Persister Genes")
dev.off()

pdf(file.path(output,"Score_of_enriched_TFs.pdf"))
ChEA3_TF_enrichment_persister %>% head(100) %>%
    mutate(TF = fct_reorder(TF, 1/Score, .desc = T)) %>%
    ggplot(aes(x=TF, y= 1/Score)) + geom_bar(fill = "#009688", stat="identity") +
    theme_classic() + theme(axis.text.x = element_text(size = 15, angle=90),
                            axis.text.y = element_text(size = 14)) + xlab("") + ylab("Inverse of ChEA3 Score")
dev.off()


pdf(file.path(output,"Score_of_enriched_TFs_overexpressed.pdf"))
ChEA3_TF_enrichment_persister %>% head(100) %>% 
    filter(TF %in% scRNA$Symbol) %>%
    mutate(TF = fct_reorder(TF, 1/Score, .desc = T)) %>%
    ggplot(aes(x=TF, y= 1/Score)) + geom_bar(fill = "#009688", stat="identity") +
    theme_classic() + theme(axis.text.x = element_text(size = 15, angle=90),
                            axis.text.y = element_text(size = 14)) + xlab("") + ylab("Inverse of ChEA3 Score")
dev.off()
