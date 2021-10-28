library(here)
# Directories -------------------------------------------------------------
maindir= here()
resdir <- file.path(maindir,"output","scRNAseq","HCC38","Persister")
source(file.path(maindir,"Scripts","global_var.R"))
source(file.path(maindir,"Scripts","functions.R"))

output <- file.path(resdir, "TFmotif") ; if(!file.exists(output)){dir.create(output)}

# Overexpressed genes
load(file.path(maindir,"output","scRNAseq","HCC38","Persister","Unsupervised","RData","LogCounts.RData"))
LogCounts = LogCounts[,grep("persister",colnames(LogCounts))]
dim(LogCounts)
BinCounts = LogCounts
BinCounts[BinCounts>0] = 1
rowSums(BinCounts)
persister_DA = readxl::read_xlsx(file.path(maindir,"output","scRNAseq","HCC38","Persister","Supervised",
                                           "Tables","Differential_analysis_Limma_logFC_1.58.xlsx"))
persister_DA$percent_expressed = rowSums(BinCounts[match(persister_DA$Symbol,rownames(BinCounts)),])/ncol(BinCounts)

persister_genes = persister_DA$Symbol
other_genes = readxl::read_xlsx(file.path(maindir,"output","scRNAseq","HCC38","Persister","Supervised",
                                          "Tables","Differential_analysis_Limma_logFC_1.58.xlsx"), sheet = 3)$Symbol
other_genes = other_genes[which(!other_genes %in% persister_genes)]

# Run ChEA3 enrichment
library(httr)
library(jsonlite)

set.seed(47)
system.time({
    list_TF_enrichment = enrich_for_TF_ChEA3(persister_genes, n_random = 100, all_genes = c(persister_genes,other_genes))
})

save(list_TF_enrichment, file = file.path(output,"list_TF_enrichment.RData"))
load(file = file.path(output,"list_TF_enrichment.RData"))

scores_per_TF = sapply(list_TF_enrichment[2:101], function(x) {
    vec = x$Score
    names(vec) = x$TF
    as.numeric(vec[order(names(vec))])
})
rownames(scores_per_TF) = sort(list_TF_enrichment[[1]]$TF)

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


# All persister genes:
pdf(file.path(output, "TF_ChEA3_enrichment_score.pdf"))
for(top in 1:20){
    scores = as.data.frame(sapply(list_TF_enrichment, function(x) x$Score[1:top]))
    scores = scores %>% tidyr::gather("Run", "Score") 
    scores = scores %>% mutate("Class" =  gsub("random_genes_.*","random_genes",Run))
    scores$Score = 1/as.numeric(scores$Score)
    
    p  = scores %>% ggplot(aes(x = Class, y = Score)) + ylab("1 / ChEA3 Mean Rank Score ") + 
        geom_violin(aes(fill = Class)) +
        ggpubr::stat_compare_means(paired = F, method = "t.test", label.x.npc = 0.4) +
        theme_classic() + ggtitle(paste0('Top ',top,' TFs'))
    print(ggpubr::add_summary(p, fun = "mean_sd"))
}
dev.off()


n_targets_explained = data.frame("TF" = paste0("TF",1:20))
n_targets_list = list()
n=0
for(i in 2:101){
    n = n+1
    ChEA3_TF_enrichment_random = list_TF_enrichment[[i]]
    n_targets_list[[n]] = sort(decreasing = T,unname(sapply(ChEA3_TF_enrichment_random$Overlapping_Genes,
                                                            function(x) length(unlist(str_split(x, pattern = ",")))))[1:20])
}
n_targets_explained$average_targets_explained = rowMeans(as.data.frame(n_targets_list))

pdf(file.path(output,"Percentage_of_targeted_persister_genes_random.pdf"))
print(n_targets_explained %>% mutate(ratio_n_targets_in_random = 100* average_targets_explained / length(persister_genes)) %>% head(20) %>%
          mutate(TF = fct_reorder(TF, ratio_n_targets_in_random, .desc = TRUE)) %>%
          ggplot(aes(x=TF, y= ratio_n_targets_in_random)) + geom_bar(fill = "#009688", stat="identity") +
          theme_classic() + theme(axis.text.x = element_text(size = 15, angle=90),
                                  axis.text.y = element_text(size = 14)) + xlab("") + ylim(c(0,60)) + ylab("% Targeted Persister Genes")
)
dev.off()

ChEA3_TF_enrichment_persister = list_TF_enrichment$Genes_of_interest
ChEA3_TF_enrichment_persister$n_targets_in_persisters = sapply(ChEA3_TF_enrichment_persister$Overlapping_Genes,
                                                               function(x) length(unlist(str_split(x, pattern = ","))))
ChEA3_TF_enrichment_persister$Score = as.numeric(ChEA3_TF_enrichment_persister$Score)

pdf(file.path(output,"Percentage_of_targeted_persister_genes.pdf"))
ChEA3_TF_enrichment_persister %>% mutate(ratio_n_targets_in_persister = 100* n_targets_in_persisters / 168) %>% head(20) %>%
    mutate(TF = fct_reorder(TF, ratio_n_targets_in_persister, .desc = TRUE)) %>%
    ggplot(aes(x=TF, y= ratio_n_targets_in_persister)) + geom_bar(fill = "#009688", stat="identity") +
    theme_classic() + theme(axis.text.x = element_text(size = 15, angle=90),
                            axis.text.y = element_text(size = 14)) + xlab("") + ylab("% Targeted Persister Genes")
dev.off()

pdf(file.path(output,"Score_of_enriched_TFs.pdf"))
ChEA3_TF_enrichment_persister %>% head(20) %>%
    mutate(TF = fct_reorder(TF, 1/Score, .desc = T)) %>%
    ggplot(aes(x=TF, y= 1/Score)) + geom_bar(fill = "#009688", stat="identity") +
    theme_classic() + theme(axis.text.x = element_text(size = 15, angle=90),
                            axis.text.y = element_text(size = 14)) + xlab("") + ylab("Inverse of ChEA3 Score")
dev.off()

n = ChEA3_TF_enrichment_persister$n_targets_in_persisters[ChEA3_TF_enrichment_persister$TF=="FOSL1"]
tot = length(persister_genes)

pdf(file.path(output,"pie_FOSL1.pdf"))
pie(c(n, tot-n), labels =c(n, tot-n))
dev.off()

