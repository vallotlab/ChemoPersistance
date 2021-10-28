library(here)
# Directories -------------------------------------------------------------
maindir= here()
resdir <- file.path(maindir,"output","scRNAseq","MM468","UNC","Unsupervised")
source(file.path(maindir,"Scripts","global_var.R"))
source(file.path(maindir,"Scripts","functions.R"))

output <- file.path(resdir, "TF_motif") ; if(!file.exists(output)){dir.create(output)}

# Overexpressed genes in HBCx95
load(file.path(maindir,"output","scRNAseq","MM468", "UNC_5FU","Unsupervised","RData","UNC_5FU_LogCounts.RData"))
LogCounts = UNC_5FU_LogCounts[,grep("UNC",colnames(UNC_5FU_LogCounts))]
dim(LogCounts)
BinCounts = LogCounts
BinCounts[BinCounts>0] = 1

persister_DA = readxl::read_xlsx(file.path(maindir,"output","scRNAseq","MM468","UNC","Supervised",
                                                  "Tables","Differential_analysis_Limma_logFC_1.58.xlsx"),3)
persister_DA$percent_expressed = rowSums(BinCounts[match(persister_DA$Symbol,rownames(BinCounts)),])/ncol(BinCounts)
persister_DA = persister_DA %>% filter(log2FC.UNC > 1.58 & qval.UNC < 0.01)

persister_genes = persister_DA$Symbol

other_genes = readxl::read_xlsx(file.path(maindir,"output","scRNAseq","PDX","Supervised",
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

# Random Background - Rank of TF enrichment
scores_per_TF = sapply(list_TF_enrichment[2:101], function(x) {
    vec = x$Score
    names(vec) = x$TF
    as.numeric(vec[order(names(vec))])
})
rownames(scores_per_TF) = sort(list_TF_enrichment[[1]]$TF)

scores_per_TF_df = data.frame(TF = sort(list_TF_enrichment[[1]]$TF),
                              mean_score =  rowMeans(scores_per_TF),
                              median_score = rowMedians(scores_per_TF))

scores_per_TF_df$rank = order(scores_per_TF_df$mean_score, decreasing = F)
scores_per_TF_df$rank[scores_per_TF_df$TF=="ATF3"]

i = 0
pdf(file.path(output, "TF_ChEA3_enrichment_score_TF_by_TF.pdf"))
for(TF in list_TF_enrichment$Genes_of_interest$TF[1:40]){
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

i = 0
pdf(file.path(output, "TF_ChEA3_enrichment_score_TF_by_TF_with_high_rank.pdf"))
high_rank = list_TF_enrichment$Genes_of_interest[which(as.numeric(list_TF_enrichment$Genes_of_interest$Score) > 750),]
for(TF in high_rank$TF[1:10]){
    i = i +1
    score_TF_in_persister = as.numeric(high_rank$Score[i])
    
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


png(file.path(output, "top_TF_ChEA3_enrichment_by_score_in_random.png"), height = 2000, width = 2500, res = 300)
scores_per_TF_df %>% arrange(mean_score) %>% head(50) %>%
    mutate(TF = fct_reorder(TF, 1/mean_score, .desc = TRUE)) %>%
    ggplot(aes(x = TF, y = 1/mean_score)) + ylab("1 / ChEA3 average(Mean Rank Score)") + 
    geom_bar(fill = "#009688", stat="identity") +theme_classic() + 
    theme(axis.text.x = element_text(size = 12, angle=90),
          axis.text.y = element_text(size = 12)) + xlab("")
dev.off()

# Percentage of targeted persister genes by top TFs
ChEA3_TF_enrichment_persister = list_TF_enrichment$Genes_of_interest
ChEA3_TF_enrichment_persister$n_targets_in_persisters = sapply(ChEA3_TF_enrichment_persister$Overlapping_Genes,
                                                               function(x) length(unlist(str_split(x, pattern = ","))))
ChEA3_TF_enrichment_persister$Score = as.numeric(ChEA3_TF_enrichment_persister$Score)

scRNA = readxl::read_xlsx(file.path(maindir, "output","scRNAseq","MM468","UNC","Supervised","Tables","Differential_analysis_Limma_logFC_1.58.xlsx"))

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
                          axis.text.y = element_text(size = 14)) + xlab("") + ylim(c(0,60)) + ylab("% Targeted Persister Genes")
dev.off()


pdf(file.path(output,"Score_of_enriched_TFs.pdf"))
ChEA3_TF_enrichment_persister %>% head(20) %>%
    mutate(TF = forcats::fct_reorder(TF, 1/Score, .desc = T)) %>%
    ggplot(aes(x=TF, y= 1/Score)) + geom_bar(fill = "#009688", stat="identity") +
    theme_classic() + theme(axis.text.x = element_text(size = 15, angle=90),
                            axis.text.y = element_text(size = 14)) + xlab("") + ylab("Inverse of ChEA3 Mean Rank Score")
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

ChEA3_TF_enrichment_persister$CorrectedScore = 
    ChEA3_TF_enrichment_persister$Score / log2(scores_per_TF_df$mean_score[match(ChEA3_TF_enrichment_persister$TF,
                                                                                 scores_per_TF_df$TF)])
pdf(file.path(output,"CorrectedScore_of_enriched_TFs.pdf"))
ChEA3_TF_enrichment_persister %>% head(20) %>%
    mutate(TF = forcats::fct_reorder(TF, (1/CorrectedScore), .desc = T)) %>%
    ggplot(aes(x=TF, y= (1/CorrectedScore))) + geom_bar(fill = "#009688", stat="identity") +
    theme_classic() + theme(axis.text.x = element_text(size = 15, angle=90),
                            axis.text.y = element_text(size = 14)) + xlab("") + ylab("Inverse of ChEA3 Mean Rank Score Corrected by Random MeanRank")
dev.off()

bivalence_in_HBCx95 = read.csv(file.path(maindir,"output","scChIPseq","PDX", "Bivalence","K27_K4_df.csv"))
bivalence_in_HBCx95 = bivalence_in_HBCx95 %>% filter(K27 > 2 & K4 > 2 )
bivalence_in_HBCx95 = bivalence_in_HBCx95 %>% group_by(Gene) %>% slice_max(K27 + K4)

persister_DA$Is_bivalent = FALSE
persister_DA$Is_bivalent[which(persister_DA$Symbol %in% bivalence_in_HBCx95$Gene)] = TRUE

# Heatmap of TF targets in persister genes for top 10 to 100 TFs
for(top in 20){
    top_TF_targets = ChEA3_TF_enrichment_persister[1:top,c("TF","Overlapping_Genes","n_targets_in_persisters","Rank","Score")]
    bin_mat = matrix(0, nrow = length(persister_genes), ncol = top,
                     dimnames = list(persister_genes, top_TF_targets$TF))
    
    targets = sapply(top_TF_targets$Overlapping_Genes, function(x) unlist(str_split(x, pattern = ",")))
    for(i in 1:top){
        bin_mat[match(targets[[i]],rownames(bin_mat)),i] = 1
        # TF_genes_SCENIC. = TF_genes_SCENIC %>% filter(TF==colnames(bin_mat)[i])
        # bin_mat[which(rownames(bin_mat) %in% TF_genes_SCENIC.$gene),i] = 1
    }
    
    summary(rowSums(bin_mat))
    bin_mat = bin_mat[which(rowSums(bin_mat) > 0),]
    
    set.seed(47)
    hc_TF = hclust(as.dist(1-jaccard(t(bin_mat))), method = methHC)
    hc_persisters = hclust(as.dist(1-jaccard(bin_mat)),method = methHC)
    
    mat.ordered <- bin_mat[hc_persisters$order,hc_TF$order]
    
    top_TF_targets$Rank = as.numeric(top_TF_targets$Rank)
    anocol_TF = geco.annotToCol(top_TF_targets[match(colnames(mat.ordered),top_TF_targets$TF),c(3,4)])
    rownames(anocol_TF) = colnames(mat.ordered)
    
    annot_persister = persister_DA[match(rownames(mat.ordered),persister_DA$Symbol),c("Symbol","percent_expressed")]
    annot_persister$gene_modules = paste0("C",cutree(hc_persisters, 4)[match(annot_persister$Symbol, names(cutree(hc_persisters, 4)))])
    anocol_persister = geco.annotToCol4(as.data.frame(annot_persister))
    rownames(anocol_persister) = rownames(mat.ordered)
    
    png(file.path(output,paste0("Heatmap_top_",top,"_TFs.png")), height=2500,width=2500,res=300)
    print(geco.hclustAnnotHeatmapPlot.withColumn(
        x=mat.ordered,
        hc=hc_TF,
        hmColors=hmColors,
        anocol=anocol_TF,
        hc_row = hc_persisters,
        anorow = anocol_persister[,2:3],
        dendro.cex=0.01,
        xlab.cex=0.75,
        hmRowNames=TRUE,
        hmRowNames.cex=0.25,
        hmColNames=TRUE,
        hmColNames.cex=0.25,
        hmCategNamesRows = TRUE,
        hmCategNamesRows.cex = 0.5))
    dev.off()
}

n_unique_targets = c()
unique_targets = c()
for(i in 1:5){
    
    n_unique_targets = c(n_unique_targets,
                         length(setdiff(rownames(mat.ordered[which(mat.ordered[,i]==1),]),
                                        unique_targets))
    )
    names(n_unique_targets)[i] = colnames(mat.ordered)[i]
    unique_targets = c(unique_targets,
                       setdiff(rownames(mat.ordered[which(mat.ordered[,i]==1),]),
                               unique_targets)
    )
}


png(file.path(output,"Venn_diagramm_FOSL1_TP63_persisters.png"), height=2500,width=2500,res=300)
gplots::venn(list("FOSL1" = rownames(mat.ordered[which(mat.ordered[,"FOSL1"]==1),]),
                  "TP63" = rownames(mat.ordered[which(mat.ordered[,"TP63"]==1),])),
             universe = persister_genes)
dev.off()

png(file.path(output,"Pie_FOSL1.png"), height=2500,width=2500,res=300)
pie(c(length(rownames(mat.ordered[which(mat.ordered[,"FOSL1"]==1),])),
      length(persister_genes) -length(rownames(mat.ordered[which(mat.ordered[,"FOSL1"]==1),]))),
    labels = c("FOSL1", "All"), init.angle = 296, col = c("#009688","#afafafff"))
dev.off()

png(file.path(output,"Pie_TP63.png"), height=2500,width=2500,res=300)
pie(c(length(rownames(mat.ordered[which(mat.ordered[,"TP63"]==1),])),
      length(persister_genes) -length(rownames(mat.ordered[which(mat.ordered[,"TP63"]==1),]))),
    labels = c("TP63", "All"), init.angle = 276, col = c("#009688","#afafafff"))
dev.off()

png(file.path(output,"Percent_PDX_persisters_targeted_top5.png"), height=2500,width=2500,res=300)
df = data.frame(TF = names(sort(n_unique_targets,decreasing = T)),
                additional_targets = cumsum(sort(n_unique_targets,decreasing = T)))
df$total_targets =  colSums(mat.ordered)[match(df$TF,colnames(mat.ordered))]
df %>% ggplot(aes(x= forcats::fct_reorder(TF, additional_targets), y = 100*additional_targets/length(persister_genes))) +
    geom_point() + theme_classic() + ylim(0,110) +
    xlab("") + ylab("% PDX persisters targeted") +
    geom_hline(yintercept = 100, lty = 2, color="red") + 
    theme(text = element_text(size = 15)) + 
    geom_text(aes(label = total_targets), nudge_y = 5)
dev.off()



n = ChEA3_TF_enrichment_persister$n_targets_in_persisters[ChEA3_TF_enrichment_persister$TF=="FOSL1"]
tot = length(persister_genes)

pdf(file.path(output,"pie_FOSL1.pdf"))
pie(c(n, tot-n), labels =c(n, tot-n))
dev.off()
