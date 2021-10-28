# Intersection K27 -> K4 / K4 -> K27 MM468
library(here)

maindir= here()
resdir <- file.path(maindir,"output","bulk_ChIPseq","MM468","ChIPreChIP")
if(!dir.exists(resdir)) dir.create(resdir)
source(file.path(maindir,"Scripts","global_var.R"))

outputDir = file.path(maindir, "output","bulk_ChIPseq","MM468","ChIPreChIP")
classes = c("c2_curated", "c5_GO","hallmark")
pattern = "BREAST|MAMMARY|HALLMARK"

# K27 K4 files
qval_K27_K4 = 0.001
odd_ratio_threshold_K27_K4 = 0.15

summits_K27_K4_peak_significant = rtracklayer::import(file.path(outputDir,paste0("summits_K27_K4_peak_significant_",qval_K27_K4,"_",odd_ratio_threshold_K27_K4,".bed")))
fisher_pvalue_df_K27_K4 = read.csv(file.path(outputDir,paste0("fisher_pvalue_K27_K4_peak_",qval_K27_K4,"_",odd_ratio_threshold_K27_K4,".csv")))
bivalent_genes_K27_K4 = unique(unlist(strsplit(fisher_pvalue_df_K27_K4$gene[fisher_pvalue_df_K27_K4$Significant],split = ",")))
length(summits_K27_K4_peak_significant)
length(bivalent_genes_K27_K4)

Bivalent_pathways_K27_K4_all = readxl::read_xlsx(file.path(outputDir,"Enrichment_test_Bivalent_pathways_K27_K4.xlsx"))
Bivalent_pathways_K27_K4_all = Bivalent_pathways_K27_K4_all[which(Bivalent_pathways_K27_K4_all$Class %in% classes),]
Bivalent_pathways_K27_K4_all$`q-value` = as.numeric(Bivalent_pathways_K27_K4_all$`q-value`)
Bivalent_pathways_K27_K4_all$log10.qvalue = -log10(Bivalent_pathways_K27_K4_all$`q-value`)
Bivalent_pathways_K27_K4_all = Bivalent_pathways_K27_K4_all[which(Bivalent_pathways_K27_K4_all$`q-value` < 0.1),]
Bivalent_pathways_K27_K4_filtered = Bivalent_pathways_K27_K4_all[grep(pattern, Bivalent_pathways_K27_K4_all$Gene.Set),]


# K4 K27 files
qval_K4_K27 = 0.001
odd_ratio_threshold_K4_K27 = 4

summits_K4_K27_peak_significant = rtracklayer::import(file.path(outputDir,paste0("summits_K4_K27_peak_significant_",qval_K4_K27,"_",odd_ratio_threshold_K4_K27,".bed")))
fisher_pvalue_df_K4_K27 = read.csv(file.path(outputDir,paste0("fisher_pvalue_K4_K27_peak_",qval_K4_K27,"_",odd_ratio_threshold_K4_K27,".csv")))
# bivalent_genes_K4_K27 = unique(unlist(strsplit(fisher_pvalue_df_K4_K27$gene[fisher_pvalue_df_K4_K27$Significant],split = ",")))
length(summits_K4_K27_peak_significant)
length(bivalent_genes_K4_K27)

Bivalent_pathways_K4_K27_all = readxl::read_xlsx(file.path(outputDir,"Enrichment_test_Bivalent_pathways_K4_K27.xlsx"))
Bivalent_pathways_K4_K27_all = Bivalent_pathways_K4_K27_all[which(Bivalent_pathways_K4_K27_all$Class %in% classes),]
Bivalent_pathways_K4_K27_all$`q-value` = as.numeric(Bivalent_pathways_K4_K27_all$`q-value`)
Bivalent_pathways_K4_K27_all$log10.qvalue = -log10(Bivalent_pathways_K4_K27_all$`q-value`)
Bivalent_pathways_K4_K27_all = Bivalent_pathways_K4_K27_all[which(Bivalent_pathways_K4_K27_all$`q-value` < 0.1),]
Bivalent_pathways_K4_K27_filtered = Bivalent_pathways_K4_K27_all[grep(pattern, Bivalent_pathways_K4_K27_all$Gene.Set),]

annot_10k = read.table("annotation/gencode.v34.annotation.transcriptTSS_10k.bed")
all_genes = unique(annot_10k$V5)

# Fisher test
A_B = length(intersect(bivalent_genes_K27_K4,bivalent_genes_K4_K27))
A_nonB = length(bivalent_genes_K27_K4) - A_B
B_nonA =  length(bivalent_genes_K4_K27) - A_B
Omega = length(all_genes) - length(union(bivalent_genes_K27_K4,bivalent_genes_K4_K27))

contingency_table_A_U_B = data.frame("Inside B" = c(A_B,
                                     B_nonA),
                               "Outside B" = c(A_nonB,
                                     Omega), row.names = c("A","All Genes"))

fisher.test(contingency_table_A_U_B,alternative = "greater")

# Venn 
png(file.path(outputDir,"Intersection_bivalent_genes.png"),width = 1200, height = 1200, res = 200)
gplots::venn(
    list("MM468_K27_K4" = bivalent_genes_K27_K4,
         "MM468_K4_K27" = bivalent_genes_K4_K27),
    show.plot = T)
# text(200,50,"(K27 K4) & (K4 K27) vs All genes : Fisher test, p.value <  2.2e-16", cex = 0.5)
dev.off()

# fisher_pvalue_df_K27_K4$log10.qvalue = -log10(fisher_pvalue_df_K27_K4$qvalue)
# fisher_pvalue_df_K4_K27$log10.qvalue = -log10(fisher_pvalue_df_K4_K27$qvalue)
# 
# common_genes = intersect(fisher_pvalue_df_K27_K4$gene,
#                          fisher_pvalue_df_K4_K27$gene)
# 
# genes_combined = cbind(
#     fisher_pvalue_df_K27_K4[match(common_genes,fisher_pvalue_df_K27_K4$gene),c(7,3)],
#     fisher_pvalue_df_K4_K27[match(common_genes,fisher_pvalue_df_K4_K27$gene),c(3)]
# )
# colnames(genes_combined)[2:3] = c("K27_K4","K4_K27")
# 
# genes_combined = genes_combined[-which(is.infinite(genes_combined$K27_K4) | is.infinite(genes_combined$K4_K27)),]
# genes_combined = genes_combined[-which(is.na(genes_combined$K27_K4) | is.na(genes_combined$K4_K27)),]
# 
# png(file.path(outputDir,"DotPlot_K27_K4_vs_K4_K27_log10_qval_AllGenes.png"),width = 3000, height = 3000, res = 300)
# genes_combined %>% 
#     ggplot(aes(x = K27_K4, y =K4_K27 )) + geom_point() +
#     theme_classic() + theme(legend.position = "none",
#                             axis.text.x = element_text(angle=90))  +
#     ggrepel::geom_text_repel(aes(label=gene), box.padding = 1) +
#     annotate(x = 6, y = 60, geom = "text",
#              label = paste("Rho = ", round(cor(genes_combined$K27_K4, genes_combined$K4_K27),2) ))
# dev.off()

################################################################################
################################################################################
# Pathways filtered
A_B = length(intersect(Bivalent_pathways_K27_K4_filtered$Gene.Set,Bivalent_pathways_K4_K27_filtered$Gene.Set))
A_nonB = length(Bivalent_pathways_K27_K4_filtered$Gene.Set) - A_B
B_nonA =  length(Bivalent_pathways_K4_K27_filtered$Gene.Set) - A_B

library(msigdbr)

MSIG.gs_filtered = MSIG.gs[which(MSIG.gs$Class %in% classes ),]
MSIG.gs_filtered = MSIG.gs_filtered[grep(pattern, MSIG.gs_filtered$Gene.Set),]
length(unique(MSIG.gs_filtered$Gene.Set))
length(unique(MSIG.gs$Gene.Set))
Omega = length(unique(MSIG.gs_filtered$Gene.Set)) - length(union(Bivalent_pathways_K27_K4_filtered$Gene.Set,Bivalent_pathways_K4_K27_filtered$Gene.Set))

contingency_table_A_U_B = data.frame("Inside B" = c(A_B,
                                                    B_nonA),
                                     "Outside B" = c(A_nonB,
                                                     Omega), row.names = c("A","All Genes"))

fisher.test(contingency_table_A_U_B,alternative = "greater")


png(file.path(outputDir,"Intersection_bivalent_pathways.png"),width = 1200, height = 1200, res = 200)
gplots::venn(
    list("MM468_K27_K4" = Bivalent_pathways_K27_K4_filtered$Gene.Set,
         "MM468_K4_K27" = Bivalent_pathways_K4_K27_filtered$Gene.Set,
         "All (Hallmark|Breast|Mammary) Gene Sets" = MSIG.gs_filtered$Gene.Set),
    show.plot = T)
text(200,50,"(K27 K4) & (K4 K27) vs All Gene Sets: Fisher test, p.value <  2.2e-16", cex = 0.8)
dev.off()

common_GS = intersect(Bivalent_pathways_K27_K4_filtered$Gene.Set, Bivalent_pathways_K4_K27_filtered$Gene.Set)

GSA_combined = cbind(
    Bivalent_pathways_K27_K4_filtered[match(common_GS,Bivalent_pathways_K27_K4_filtered$Gene.Set),c(1,8)],
    Bivalent_pathways_K4_K27_filtered[match(common_GS,Bivalent_pathways_K4_K27_filtered$Gene.Set),c(8)]
)
colnames(GSA_combined)[2:3] = c("K27_K4","K4_K27")

png(file.path(outputDir,"DotPlot_K27_K4_vs_K4_K27_Bivalent_Pathways.png"),width = 3000, height = 3000, res = 300)
GSA_combined %>% 
    ggplot(aes(x = K27_K4, y =K4_K27 )) + geom_point() +
    theme_classic() + theme(legend.position = "none",
                            axis.text.x = element_text(angle=90))  +
    ggrepel::geom_text_repel(aes(label=Gene.Set), box.padding = 1) +
    annotate(x = 20, y = 15, geom = "text",
             label = paste("Rho = ", round(cor(GSA_combined$K27_K4, GSA_combined$K4_K27),2) ))
dev.off()

# Pathways non filtered
A_B = length(intersect(Bivalent_pathways_K27_K4_all$Gene.Set,Bivalent_pathways_K4_K27_all$Gene.Set))
A_nonB = length(Bivalent_pathways_K27_K4_all$Gene.Set) - A_B
B_nonA =  length(Bivalent_pathways_K4_K27_all$Gene.Set) - A_B

MSIG.gs_all = MSIG.gs[which(MSIG.gs$Class %in% classes ),]
length(unique(MSIG.gs_all$Gene.Set))

Omega = length(unique(MSIG.gs_all$Gene.Set)) - length(union(Bivalent_pathways_K27_K4_all$Gene.Set,Bivalent_pathways_K4_K27_all$Gene.Set))

contingency_table_A_U_B = data.frame("Inside B" = c(A_B,
                                                    B_nonA),
                                     "Outside B" = c(A_nonB,
                                                     Omega), row.names = c("A","All Genes"))
contingency_table_A_U_B
fisher.test(contingency_table_A_U_B,alternative = "greater")


png(file.path(outputDir,"Intersection_bivalent_pathways_all.png"),width = 1200, height = 1200, res = 200)
gplots::venn(
    list("MM468_K27_K4" = Bivalent_pathways_K27_K4_all$Gene.Set,
         "MM468_K4_K27" = Bivalent_pathways_K4_K27_all$Gene.Set,
         "All (C2|C5|Hallmark) Gene Sets" = MSIG.gs_all$Gene.Set),
    show.plot = T)
text(200,50,"(K27 K4) & (K4 K27) vs All Gene Sets: Fisher test, p.value <  2.2e-16", cex = 0.8)
dev.off()

common_GS = intersect(Bivalent_pathways_K27_K4_all$Gene.Set, Bivalent_pathways_K4_K27_all$Gene.Set)

GSA_combined = cbind(
    Bivalent_pathways_K27_K4_all[match(common_GS,Bivalent_pathways_K27_K4_all$Gene.Set),c(1,8)],
    Bivalent_pathways_K4_K27_all[match(common_GS,Bivalent_pathways_K4_K27_all$Gene.Set),c(8)]
)
colnames(GSA_combined)[2:3] = c("K27_K4","K4_K27")

png(file.path(outputDir,"DotPlot_K27_K4_vs_K4_K27_Bivalent_Pathways_all.png"),width = 3000, height = 3000, res = 300)
GSA_combined %>% 
    ggplot(aes(x = K27_K4, y =K4_K27 )) + geom_point() +
    theme_classic() + theme(legend.position = "none",
                            axis.text.x = element_text(angle=90))  +
    ggrepel::geom_text_repel(aes(label=Gene.Set), box.padding = 1) +
    annotate(x = 100, y = 15, geom = "text",
             label = paste("Rho = ", round(cor(GSA_combined$K27_K4, GSA_combined$K4_K27),2) ))
dev.off()

lower_than_500 = Bivalent_pathways_K4_K27_all$Gene.Set[which(Bivalent_pathways_K4_K27_all$Nb_of_genes < 500)]
GSA_combined = GSA_combined[which(GSA_combined$Gene.Set %in% lower_than_500),]

png(file.path(outputDir,"DotPlot_K27_K4_vs_K4_K27_Bivalent_Pathways_all_lower_than_500.png"),width = 2000, height = 2000, res = 300)
GSA_combined %>% 
    ggplot(aes(x = K27_K4, y =K4_K27 )) + geom_point() +
    theme_classic() + theme(legend.position = "none",
                            axis.text.x = element_text(angle=90))  +
    ggrepel::geom_text_repel(aes(label=Gene.Set), box.padding = 1) +
    annotate(x = 20, y = 5, geom = "text",
             label = paste("Rho = ", round(cor(GSA_combined$K27_K4, GSA_combined$K4_K27),2) ))
dev.off()


GSA_combined = rbind(Bivalent_pathways_K27_K4_filtered, Bivalent_pathways_K4_K27_filtered)

GSA_combined_save =  GSA_combined
GSA_combined_save = GSA_combined_save %>% group_by(Gene.Set) %>% mutate(occurence = n())

GSA_combined_save = GSA_combined_save %>% filter(occurence > 1) %>%
    group_by(Gene.Set) %>% summarise(
        Nb_of_genes = mean(Nb_of_genes),
        Nb_of_deregulated_genes = round(mean(Nb_of_deregulated_genes)),
        Class = unique(Class),
        mean_qval = mean(`q-value`),
        mean_log10qval = mean(log10.qvalue)
    )
GSA_combined_save = GSA_combined_save %>% arrange(desc(mean_log10qval))
WriteXLS::WriteXLS(GSA_combined_save, ExcelFileName = file.path(outputDir,"combined_GSA_bivalent_MM468_K4_K27_K27_K4.xlsx"))
