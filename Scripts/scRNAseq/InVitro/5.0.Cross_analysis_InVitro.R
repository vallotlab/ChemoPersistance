library(here)
# Directories -------------------------------------------------------------
maindir= here()
resdir <- file.path(maindir,"output","scRNAseq","MM468")
source(file.path(maindir,"Scripts","global_var.R"))
source(file.path(maindir,"Scripts","functions.R"))

output <- file.path(resdir, "CrossAnalysis_inVitro2") ; if(!file.exists(output)){dir.create(output)}

################################################################################
# scRNA persister genes
################################################################################

# Overexpressed genes in MM468
load(file.path(maindir,"output","scRNAseq","MM468","Persister","Unsupervised","RData","LogCounts.RData"))
LogCounts = LogCounts[,grep("persister",colnames(LogCounts))]
dim(LogCounts)
BinCounts = LogCounts
BinCounts[BinCounts>0] = 1
rowSums(BinCounts)

persister_DA_MM468 = readxl::read_xlsx(file.path(maindir,"output","scRNAseq","MM468","Persister","Supervised",
                                                 "Tables","Differential_analysis_Limma_logFC_1.58.xlsx"),3)
persister_DA_MM468$percent_expressed_MM468 = rowSums(BinCounts[match(persister_DA_MM468$Symbol,rownames(BinCounts)),])/ncol(BinCounts)
persister_DA_MM468$qval.C2_pers = as.numeric(persister_DA_MM468$qval.C2_pers)
persister_DA_MM468 = persister_DA_MM468 %>% filter(log2FC.C2_pers > log2(3) & qval.C2_pers < 0.01 )
persister_genes_MM468 = persister_DA_MM468$Symbol

# Overexpressed genes in BT20
load(file.path(maindir,"output","scRNAseq","BT20","Persister","Unsupervised","RData","LogCounts.RData"))
LogCounts = LogCounts[,grep("persister",colnames(LogCounts))]
dim(LogCounts)
BinCounts = LogCounts
BinCounts[BinCounts>0] = 1

persister_DA_BT20 = readxl::read_xlsx(file.path(maindir,"output","scRNAseq","BT20","Persister","Supervised",
                                                 "Tables","Differential_analysis_Limma_logFC_1.58.xlsx"),3)
persister_DA_BT20$qval.BT20 = as.numeric(persister_DA_BT20$qval.persister)
persister_DA_BT20$percent_expressed_BT20 = rowSums(BinCounts[match(persister_DA_BT20$Symbol,rownames(BinCounts)),])/ncol(BinCounts)
persister_DA_BT20 = persister_DA_BT20 %>% filter(log2FC.persister > log2(3) & qval.persister   < 0.01 )
persister_genes_BT20 = persister_DA_BT20$Symbol

# Overexpressed genes in HCC38
load(file.path(maindir,"output","scRNAseq","HCC38","Persister","Unsupervised","RData","LogCounts.RData"))
LogCounts = LogCounts[,grep("persister",colnames(LogCounts))]
dim(LogCounts)
BinCounts = LogCounts
BinCounts[BinCounts>0] = 1

persister_DA_HCC38 = readxl::read_xlsx(file.path(maindir,"output","scRNAseq","HCC38","Persister","Supervised",
                                                  "Tables","Differential_analysis_Limma_logFC_1.58.xlsx"),3)
persister_DA_HCC38$qval.HCC38 = as.numeric(persister_DA_HCC38$qval.persister)
persister_DA_HCC38$percent_expressed_HCC38 = rowSums(BinCounts[match(persister_DA_HCC38$Symbol,rownames(BinCounts)),])/ncol(BinCounts)
persister_DA_HCC38 = persister_DA_HCC38 %>% filter(log2FC.persister > log2(3) & qval.persister   < 0.01)
persister_genes_HCC38 = persister_DA_HCC38$Symbol

persister_genes = intersect(intersect(persister_genes_MM468, persister_genes_BT20),persister_genes_HCC38)

common_genes = intersect(intersect(persister_DA_HCC38$Symbol, persister_DA_BT20$Symbol), persister_DA_MM468$Symbol)

colnames(persister_DA_MM468) = gsub("C2_pers","MM468", colnames(persister_DA_MM468))
colnames(persister_DA_BT20) = gsub("persister","BT20", colnames(persister_DA_BT20))
colnames(persister_DA_HCC38) = gsub("persister","HCC38", colnames(persister_DA_HCC38))

# Venn diagrams of genes
png(file.path(output,"Intersection_InVitro_persister_genes.png"),width = 1200, height = 1200, res = 200)
gplots::venn(
    list("MM468" = persister_genes_MM468,
         "BT20" = persister_genes_BT20,
         "HCC38" = persister_genes_HCC38),
    show.plot = T)
dev.off()

################################################################################
# scRNA persister pathways
################################################################################

# Load each Gene Set Analysis # Filter
classes = c("c2_curated", "c5_GO","hallmark")
pattern = "BREAST|MAMMARY|HALLMARK"

# MM468
GSA_MM468 = readxl::read_xlsx(file.path(maindir,"output","scRNAseq","MM468","Persister","Supervised",
                                        "Tables","Gene_set_analysis", "Enrichment_test_C2_pers_vs_DMSO_MSigDB_logFC1.58.xlsx"))
GSA_MM468$log10.qvalue = -log10(GSA_MM468$`q-value`)
colnames(GSA_MM468)[c(4,6,7,8)] = paste0("InVitro_MM468.",colnames(GSA_MM468)[c(4,6,7,8)])

# BT20
GSA_BT20 = readxl::read_xlsx(file.path(maindir,"output","scRNAseq","BT20","Persister","Supervised",
                                        "Tables","Gene_set_analysis", "Enrichment_test_persister_vs_chemonaive_MSigDB_logFC1.58.xlsx"))
GSA_BT20$log10.qvalue = -log10(GSA_BT20$`q-value`)
colnames(GSA_BT20)[c(4,6,7,8)] = paste0("InVitro_BT20.",colnames(GSA_BT20)[c(4,6,7,8)])

# HCC38
GSA_HCC38 = readxl::read_xlsx(file.path(maindir,"output","scRNAseq","HCC38","Persister","Supervised",
                                         "Tables","Gene_set_analysis", "Enrichment_test_persister_vs_chemonaive_MSigDB_logFC1.58.xlsx"))
GSA_HCC38$log10.qvalue = -log10(GSA_HCC38$`q-value`)
colnames(GSA_HCC38)[c(4,6,7,8)] = paste0("InVitro_HCC38.",colnames(GSA_HCC38)[c(4,6,7,8)])

# Combining GSA together
GSA_BT20 = GSA_BT20[order(GSA_BT20$Gene.Set),]
GSA_combined = cbind(GSA_MM468[match(GSA_BT20$Gene.Set, GSA_MM468$Gene.Set),1:3],
                     GSA_MM468[match(GSA_BT20$Gene.Set, GSA_MM468$Gene.Set),c(4,6,7,8)],
                     GSA_BT20[,c(4,6,7,8)],
                     GSA_HCC38[match(GSA_BT20$Gene.Set, GSA_HCC38$Gene.Set),c(4,6,7,8)])
GSA_combined =  mutate(GSA_combined,
                       mean_qvalue = rowMeans(dplyr::select(GSA_combined, ends_with("q-value")),na.rm = T),
                       mean_log10.qvalue = rowMeans(dplyr::select(GSA_combined, ends_with("log10.qvalue")),na.rm = T),
                       mean_Nb_Deregulated_Genes = rowMeans(dplyr::select(GSA_combined, ends_with("Nb_of_deregulated_genes")),na.rm = T)
)

GSA_combined = GSA_combined[,c(1:3,16:18,4:15)]
GSA_combined = GSA_combined %>% arrange(mean_qvalue)

# Formating GSA for saving
GSA_combined_save =  GSA_combined
GSA_combined_save = GSA_combined_save %>% group_by(Gene.Set) %>% mutate(occurence = n())

GSA_combined_curated = GSA_combined_save[grep(pattern, GSA_combined_save$Gene.Set),]
GSA_combined_curated = GSA_combined_curated[grep("AMPLICON", invert = TRUE, GSA_combined_curated$Gene.Set),]
GSA_combined_curated = GSA_combined_curated[which(GSA_combined_curated$Class %in% classes),]

png(file.path(output,"Intersection_InVitro_persister_GSA.png"),width = 1200, height = 1200, res = 200)
gplots::venn(
    list("MM468" = GSA_combined_curated$Gene.Set[which(GSA_combined_curated$`InVitro_MM468.q-value` < 0.1)],
         "BT20" = GSA_combined_curated$Gene.Set[which(GSA_combined_curated$`InVitro_BT20.q-value` < 0.1)],
         "HCC38" = GSA_combined_curated$Gene.Set[which(GSA_combined_curated$`InVitro_HCC38.q-value` < 0.1)]),
    show.plot = T)
dev.off()

library(SuperExactTest)
SuperExactTest::MSET(list("MM468" = GSA_combined_curated$Gene.Set[which(GSA_combined_curated$`InVitro_MM468.q-value` < 0.1)],
                          "BT20" = GSA_combined_curated$Gene.Set[which(GSA_combined_curated$`InVitro_BT20.q-value` < 0.1)],
                          "HCC38" = "EMPTY"),
                     n = length(unique(GSA_combined_curated$Gene.Set)),
                     lower.tail = FALSE)

GSA_combined_save = GSA_combined_save %>% arrange(mean_qvalue)
WriteXLS::WriteXLS(GSA_combined_save, ExcelFileName = file.path(output,"combined_GSA_scRNA_InVitros_unfiltered.xlsx"))

################################################################################
# scRNA TF enrichment (ChEA3)
################################################################################

# Load each CheA3 analysis separately
load(file.path(maindir,"output","scRNAseq","MM468","Persister","Unsupervised","TFmotif","list_TF_enrichment.RData") )
TF_MM468_all = list_TF_enrichment$Genes_of_interest
TF_MM468 = head(list_TF_enrichment$Genes_of_interest,100)

load(file.path(maindir,"output","scRNAseq","BT20","Persister","TFmotif","list_TF_enrichment.RData") )
TF_BT20_all = list_TF_enrichment$Genes_of_interest
TF_BT20 = head(list_TF_enrichment$Genes_of_interest,100)

load(file.path(maindir,"output","scRNAseq","HCC38","Persister","TFmotif","list_TF_enrichment.RData") )
TF_HCC38_all = list_TF_enrichment$Genes_of_interest
TF_HCC38 = head(list_TF_enrichment$Genes_of_interest,100)

png(file.path(output,"Intersection_InVitro_persister_TF_motif.png"),width = 1200, height = 1200, res = 200)
gplots::venn(
    list("MM468" = TF_MM468$TF,
         "BT20" = TF_BT20$TF,
         "HCC38" = TF_HCC38$TF),
    show.plot = T)
dev.off()

SuperExactTest::MSET(list("MM468" = TF_MM468$TF,
                          "BT20" = TF_BT20$TF,
                          "HCC38" = TF_HCC38$TF),
                     n = length(unique(TF_HCC38_all$TF)),
                     lower.tail = FALSE
)

################################################################################
# Sequential ChIP-seq - bivalent genes
################################################################################

# MM468
Bivalent_MM468 = read.csv(file.path("output","bulk_ChIPseq","MM468","ChIPreChIP",paste0("fisher_pvalue_K4_K27_peak_0.001_4.csv")))
Bivalent_genes_MM468  = unique(unlist(strsplit(Bivalent_MM468 $gene[which(Bivalent_MM468$Significant)], split = ",")))

# BT20
Bivalent_BT20 = read.csv(file.path("output","bulk_ChIPseq","BT20","ChIPreChIP",paste0("fisher_pvalue_K4_K27_peak_0.001_4.csv")))
Bivalent_genes_BT20  = unique(unlist(strsplit(Bivalent_BT20 $gene[which(Bivalent_BT20$Significant)], split = ",")))


# HCC38
Bivalent_HCC38 = read.csv(file.path("output","bulk_ChIPseq","HCC38","ChIPreChIP",paste0("fisher_pvalue_K4_K27_peak_0.001_4.csv")))
Bivalent_genes_HCC38 = unique(unlist(strsplit(Bivalent_HCC38$gene[which(Bivalent_HCC38$Significant)], split = ",")))

png(file.path(output,"Intersection_InVitro_bivalent_genes.png"),width = 1200, height = 1200, res = 200)
gplots::venn(
    list("MM468" = Bivalent_genes_MM468,
         "BT20" = Bivalent_genes_BT20,
         "HCC38" = Bivalent_genes_HCC38),
    show.plot = T)
dev.off()

annot10k = read.table("annotation/gencode.v34.annotation.transcriptTSS_10k.bed")
SuperExactTest::MSET(list("MM468" = Bivalent_genes_MM468,
                          "BT20" = Bivalent_genes_BT20,
                          "HCC38" = Bivalent_genes_HCC38),
                     n = length(unique(annot10k$V5)),
                     lower.tail = FALSE
)

################################################################################
# Sequential ChIP-seq - bivalent pathways
################################################################################

# Load each Bivalent analysis separately
classes = c("c2_curated", "c5_GO","hallmark")
pattern = "BREAST|MAMMARY|HALLMARK"

# MM468
Bivalent_GSA_MM468 = readxl::read_xlsx(file.path("output","bulk_ChIPseq","MM468","ChIPreChIP","Enrichment_test_Bivalent_pathways_K4_K27.xlsx"))
Bivalent_GSA_MM468$log10.qvalue = -log10(Bivalent_GSA_MM468$`q-value`)
colnames(Bivalent_GSA_MM468)[c(4,6,7,8)] = paste0("InVitro_MM468.",colnames(Bivalent_GSA_MM468)[c(4,6,7,8)])

# BT20
Bivalent_GSA_BT20 = readxl::read_xlsx(file.path("output","bulk_ChIPseq","BT20","ChIPreChIP","Enrichment_test_Bivalent_pathways_K4_K27.xlsx"))
Bivalent_GSA_BT20$log10.qvalue = -log10(Bivalent_GSA_BT20$`q-value`)
colnames(Bivalent_GSA_BT20)[c(4,6,7,8)] = paste0("InVitro_BT20.",colnames(Bivalent_GSA_BT20)[c(4,6,7,8)])

# HCC38
Bivalent_GSA_HCC38 = readxl::read_xlsx(file.path("output","bulk_ChIPseq","HCC38","ChIPreChIP","Enrichment_test_Bivalent_pathways_K4_K27.xlsx"))
Bivalent_GSA_HCC38$log10.qvalue = -log10(Bivalent_GSA_HCC38$`q-value`)
colnames(Bivalent_GSA_HCC38)[c(4,6,7,8)] = paste0("InVitro_HCC38.",colnames(Bivalent_GSA_HCC38)[c(4,6,7,8)])


Bivalent_GSA_BT20 = Bivalent_GSA_BT20[order(Bivalent_GSA_BT20$Gene.Set),]
GSA_combined = cbind(Bivalent_GSA_MM468[match(Bivalent_GSA_BT20$Gene.Set, Bivalent_GSA_MM468$Gene.Set),1:3],
                     Bivalent_GSA_MM468[match(Bivalent_GSA_BT20$Gene.Set, Bivalent_GSA_MM468$Gene.Set),c(4,6,7,8)],
                     Bivalent_GSA_BT20[,c(4,6,7,8)],
                     Bivalent_GSA_HCC38[match(Bivalent_GSA_BT20$Gene.Set, Bivalent_GSA_HCC38$Gene.Set),c(4,6,7,8)])

GSA_combined =  mutate(GSA_combined,
                       mean_qvalue = rowMeans(dplyr::select(GSA_combined, ends_with("q-value")),na.rm = T),
                       mean_log10.qvalue = rowMeans(dplyr::select(GSA_combined, ends_with("log10.qvalue")),na.rm = T),
                       mean_Nb_Deregulated_Genes = rowMeans(dplyr::select(GSA_combined, ends_with("Nb_of_deregulated_genes")),na.rm = T)
)

GSA_combined = GSA_combined[,c(1:3,16:18,4:15)]
GSA_combined = GSA_combined %>% arrange(mean_qvalue)

GSA_combined_save =  GSA_combined
GSA_combined_save = GSA_combined_save %>% group_by(Gene.Set) %>% mutate(occurence = n())
GSA_combined_save = GSA_combined_save %>% arrange(mean_qvalue)

GSA_combined_curated = GSA_combined[grep(pattern, GSA_combined$Gene.Set),]
GSA_combined_curated = GSA_combined_curated[grep("AMPLICON", invert = TRUE, GSA_combined_curated$Gene.Set),]
GSA_combined_curated = GSA_combined_curated[which(GSA_combined_curated$Class %in% classes),]

png(file.path(output,"Intersection_InVitro_bivalent_GSA.png"),width = 1200, height = 1200, res = 200)
gplots::venn(
    list("MM468" = GSA_combined_curated$Gene.Set[which(GSA_combined_curated$`InVitro_MM468.q-value` < 0.1)],
         "BT20" = GSA_combined_curated$Gene.Set[which(GSA_combined_curated$`InVitro_BT20.q-value` < 0.1)],
         "HCC38" = GSA_combined_curated$Gene.Set[which(GSA_combined_curated$`InVitro_HCC38.q-value` < 0.1)]),
    show.plot = T)
dev.off()

SuperExactTest::MSET(list("MM468" = GSA_combined_curated$Gene.Set[which(GSA_combined_curated$`InVitro_MM468.q-value` < 0.1)],
                          "BT20" = GSA_combined_curated$Gene.Set[which(GSA_combined_curated$`InVitro_BT20.q-value` < 0.1)],
                          "HCC38" = GSA_combined_curated$Gene.Set[which(GSA_combined_curated$`InVitro_HCC38.q-value` < 0.1)]),
                     n = length(unique(GSA_combined_curated$Gene.Set)),
                     lower.tail = FALSE)


WriteXLS::WriteXLS(GSA_combined_save, ExcelFileName = file.path(output,"combined_GSA_Bivalence_InVitros_unfiltered.xlsx"))






