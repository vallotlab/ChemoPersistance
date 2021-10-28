library(here)
# Directories -------------------------------------------------------------
maindir= here()
resdir <- file.path(maindir,"output","scRNAseq","PDX_BC976")
source(file.path(maindir,"Scripts","global_var.R"))
source(file.path(maindir,"Scripts","functions.R"))


output <- file.path(resdir, "CrossAnalysis_PDX") ; if(!file.exists(output)){dir.create(output)}

################################################################################
# scRNA persister genes
################################################################################

# Overexpressed genes in BC976
load(file.path(maindir,"output","scRNAseq","PDX_BC976","Unsupervised","RData","LogCounts.RData"))
LogCounts = LogCounts[,grep("persister",colnames(LogCounts))]
dim(LogCounts)
BinCounts = LogCounts
BinCounts[BinCounts>0] = 1
rowSums(BinCounts)

persister_DA_BC976 = readxl::read_xlsx(file.path(maindir,"output","scRNAseq","PDX_BC976","Supervised",
                                           "Tables","Differential_analysis_Limma_logFC_1.58.xlsx"),3)
persister_DA_BC976$percent_expressed_BC976 = rowSums(BinCounts[match(persister_DA_BC976$Symbol,rownames(BinCounts)),])/ncol(BinCounts)
persister_DA_BC976 = persister_DA_BC976 %>% filter(log2FC.persister_6_vs_UNT > log2(3) & qval.persister_6_vs_UNT < 0.01 )
persister_genes_BC976 = persister_DA_BC976$Symbol

# Overexpressed genes in BC408
load(file.path(maindir,"output","scRNAseq","PDX_BC408","Unsupervised","RData","LogCounts.RData"))
LogCounts = LogCounts[,grep("persister",colnames(LogCounts))]
dim(LogCounts)
BinCounts = LogCounts
BinCounts[BinCounts>0] = 1

persister_DA_BC408 = readxl::read_xlsx(file.path(maindir,"output","scRNAseq","PDX_BC408","Supervised",
                                                 "Tables","Differential_analysis_Limma_logFC_Pers_vs_chemonaive_1.58.xlsx"),3)
persister_DA_BC408$percent_expressed_BC408 = rowSums(BinCounts[match(persister_DA_BC408$Symbol,rownames(BinCounts)),])/ncol(BinCounts)
persister_DA_BC408 = persister_DA_BC408 %>% filter(log2FC.persister > log2(3) & qval.persister   < 0.01 )
persister_genes_BC408 = persister_DA_BC408$Symbol

# Overexpressed genes in BC1224
load(file.path(maindir,"output","scRNAseq","PDX_BC1224","Unsupervised","RData","LogCounts.RData"))
LogCounts = LogCounts[,grep("persister",colnames(LogCounts))]
dim(LogCounts)
BinCounts = LogCounts
BinCounts[BinCounts>0] = 1

persister_DA_BC1224 = readxl::read_xlsx(file.path(maindir,"output","scRNAseq","PDX_BC1224","Supervised",
                                                 "Tables","Differential_analysis_Limma_logFC_Pers_vs_chemonaive_1.58.xlsx"),3)
persister_DA_BC1224$percent_expressed_BC1224 = rowSums(BinCounts[match(persister_DA_BC1224$Symbol,rownames(BinCounts)),])/ncol(BinCounts)
persister_DA_BC1224 = persister_DA_BC1224 %>% filter(log2FC.persister > log2(3) & qval.persister   < 0.01)
persister_genes_BC1224 = persister_DA_BC1224$Symbol

persister_genes = intersect(intersect(persister_genes_BC976, persister_genes_BC408),persister_genes_BC1224)

common_genes = intersect(intersect(persister_DA_BC1224$Symbol, persister_DA_BC408$Symbol), persister_DA_BC976$Symbol)

colnames(persister_DA_BC976) = gsub("persister_6_vs_UNT","BC976", colnames(persister_DA_BC976))
colnames(persister_DA_BC408) = gsub("persister","BC408", colnames(persister_DA_BC408))
colnames(persister_DA_BC1224) = gsub("persister","BC1224", colnames(persister_DA_BC1224))

# Venn diagrams of genes
png(file.path(output,"Intersection_PDX_persister_genes.png"),width = 1200, height = 1200, res = 200)
gplots::venn(
    list("BC976" = persister_genes_BC976,
         "BC408" = persister_genes_BC408,
         "BC1224" = persister_genes_BC1224),
    show.plot = T)
dev.off()

################################################################################
# scRNA persister pathways
################################################################################

# Load each Gene Set Analysis # Filter
classes = c("c2_curated", "c5_GO","hallmark")
pattern = "BREAST|MAMMARY|HALLMARK"

# BC976
GSA_BC976 = readxl::read_xlsx(file.path(maindir,"output","scRNAseq","PDX_BC976","Supervised",
                                        "Tables","Gene_set_analysis", "Enrichment_test_persister_vs_UNT_vs_UNT_MSigDB_logFC1.58_unfiltered.xlsx"))
GSA_BC976$log10.qvalue = -log10(GSA_BC976$`q-value`)
colnames(GSA_BC976)[c(4,6,7,8)] = paste0("PDX_BC976.",colnames(GSA_BC976)[c(4,6,7,8)])


# BC408
GSA_BC408 = readxl::read_xlsx(file.path(maindir,"output","scRNAseq","PDX_BC408","Supervised",
                                        "Tables","Gene_set_analysis", "Enrichment_test_persister_vs_chemonaive_MSigDB_logFC1.58_unfiltered.xlsx"))
GSA_BC408$log10.qvalue = -log10(GSA_BC408$`q-value`)
colnames(GSA_BC408)[c(4,6,7,8)] = paste0("PDX_BC408.",colnames(GSA_BC408)[c(4,6,7,8)])

# BC1224
GSA_BC1224 = readxl::read_xlsx(file.path(maindir,"output","scRNAseq","PDX_BC1224","Supervised",
                                         "Tables","Gene_set_analysis", "Enrichment_test_persister_vs_chemonaive_MSigDB_logFC1.58_unfiltered.xlsx"))
GSA_BC1224$log10.qvalue = -log10(GSA_BC1224$`q-value`)
colnames(GSA_BC1224)[c(4,6,7,8)] = paste0("PDX_BC1224.",colnames(GSA_BC1224)[c(4,6,7,8)])


# Combining GSA together
GSA_BC408 = GSA_BC408[order(GSA_BC408$Gene.Set),]
GSA_combined = cbind(GSA_BC976[match(GSA_BC408$Gene.Set, GSA_BC976$Gene.Set),1:3],
                     GSA_BC976[match(GSA_BC408$Gene.Set, GSA_BC976$Gene.Set),c(4,6,7,8)],
                     GSA_BC408[,c(4,6,7,8)],
                     GSA_BC1224[match(GSA_BC408$Gene.Set, GSA_BC1224$Gene.Set),c(4,6,7,8)])
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

png(file.path(output,"Intersection_PDX_persister_GSA.png"),width = 1200, height = 1200, res = 200)
gplots::venn(
    list("BC976" = GSA_combined_curated$Gene.Set[which(GSA_combined_curated$`PDX_BC976.q-value` < 0.1)],
         "BC408" = GSA_combined_curated$Gene.Set[which(GSA_combined_curated$`PDX_BC408.q-value` < 0.1)],
         "BC1224" = GSA_combined_curated$Gene.Set[which(GSA_combined_curated$`PDX_BC1224.q-value` < 0.1)]),
    show.plot = T)
dev.off()

library(SuperExactTest)
SuperExactTest::MSET(list("BC976" = GSA_combined_curated$Gene.Set[which(GSA_combined_curated$`PDX_BC976.q-value` < 0.1)],
                          "BC408" = GSA_combined_curated$Gene.Set[which(GSA_combined_curated$`PDX_BC408.q-value` < 0.1)],
                          "BC1224" = GSA_combined_curated$Gene.Set[which(GSA_combined_curated$`PDX_BC1224.q-value` < 0.1)]),
                     n = length(unique(GSA_combined_curated$Gene.Set)),
                     lower.tail = FALSE)


GSA_combined_save = GSA_combined_save %>% arrange(desc(mean_log10qval))
WriteXLS::WriteXLS(GSA_combined_save, ExcelFileName = file.path(output,"combined_GSA_scRNA_PDXs_unfiltered.xlsx"))

################################################################################
# scRNA TF enrichment (ChEA3)
################################################################################

# Load each CheA3 analysis separately
load(file.path(maindir,"output","scRNAseq","PDX_BC976","Unsupervised","TFmotif","list_TF_enrichment.RData") )
TF_BC976_all = list_TF_enrichment$Genes_of_interest
TF_BC976 = head(list_TF_enrichment$Genes_of_interest,100)

load(file.path(maindir,"output","scRNAseq","PDX_BC408","Unsupervised","TFmotif","list_TF_enrichment.RData") )
TF_BC408_all = list_TF_enrichment$Genes_of_interest
TF_BC408 = head(list_TF_enrichment$Genes_of_interest,100)

load(file.path(maindir,"output","scRNAseq","PDX_BC1224","Unsupervised","TFmotif","list_TF_enrichment.RData") )
TF_BC1224_all = list_TF_enrichment$Genes_of_interest
TF_BC1224 = head(list_TF_enrichment$Genes_of_interest,100)


png(file.path(output,"Intersection_PDX_persister_TF_motif.png"),width = 1200, height = 1200, res = 200)
gplots::venn(
    list("BC976" = TF_BC976$TF,
         "BC408" = TF_BC408$TF,
         "BC1224" = TF_BC1224$TF),
    show.plot = T)
dev.off()

SuperExactTest::MSET(list("BC976" = TF_BC976$TF,
                          "BC408" = TF_BC408$TF,
                          "BC1224" = TF_BC1224$TF),
                     n = length(unique(TF_BC1224_all$TF)),
                     lower.tail = FALSE
                     )

################################################################################
# Sequential ChIP-seq - bivalent genes
################################################################################

# BC976
Bivalent_BC976 = read.csv(file.path("output","bulk_ChIPseq","BC976","ChIPreChIP",paste0("fisher_pvalue_K4_K27_peak_1e-15_4.csv")))
Bivalent_genes_BC976  = unique(unlist(strsplit(Bivalent_BC976 $gene[which(Bivalent_BC976$Significant)], split = ",")))

# BC408
Bivalent_BC408 = read.csv(file.path("output","bulk_ChIPseq","BC408","ChIPreChIP",paste0("fisher_pvalue_K4_K27_peak_0.001_4.csv")))
Bivalent_genes_BC408  = unique(unlist(strsplit(Bivalent_BC408 $gene[which(Bivalent_BC408$Significant)], split = ",")))


# BC1224
Bivalent_BC1224 = read.csv(file.path("output","bulk_ChIPseq","BC1224","ChIPreChIP",paste0("fisher_pvalue_K4_K27_peak_0.001_4.csv")))
Bivalent_genes_BC1224 = unique(unlist(strsplit(Bivalent_BC1224$gene[which(Bivalent_BC1224$Significant)], split = ",")))

png(file.path(output,"Intersection_PDX_bivalent_genes.png"),width = 1200, height = 1200, res = 200)
gplots::venn(
    list("BC976" = Bivalent_genes_BC976,
         "BC408" = Bivalent_genes_BC408,
         "BC1224" = Bivalent_genes_BC1224),
    show.plot = T)
dev.off()

annot10k = read.table("annotation/gencode.v34.annotation.transcriptTSS_10k.bed")
SuperExactTest::MSET(list("BC976" = Bivalent_genes_BC976,
                          "BC408" = Bivalent_genes_BC408,
                          "BC1224" = Bivalent_genes_BC1224),
                     n = length(unique(annot10k$V5)),
                     lower.tail = FALSE
)

################################################################################
# Sequential ChIP-seq - bivalent pathways
################################################################################

# Load each Bivalent analysis separately
classes = c("c2_curated", "c5_GO","hallmark")
pattern = "BREAST|MAMMARY|HALLMARK"

# BC976
Bivalent_GSA_BC976 = readxl::read_xlsx(file.path("output","bulk_ChIPseq","BC976","ChIPreChIP","Enrichment_test_Bivalent_pathways_K4_K27.xlsx"))
Bivalent_GSA_BC976$log10.qvalue = -log10(Bivalent_GSA_BC976$`q-value`)
colnames(Bivalent_GSA_BC976)[c(4,6,7,8)] = paste0("PDX_BC976.",colnames(Bivalent_GSA_BC976)[c(4,6,7,8)])

# BC408
Bivalent_GSA_BC408 = readxl::read_xlsx(file.path("output","bulk_ChIPseq","BC408","ChIPreChIP","Enrichment_test_Bivalent_pathways_K4_K27.xlsx"))
Bivalent_GSA_BC408$log10.qvalue = -log10(Bivalent_GSA_BC408$`q-value`)
colnames(Bivalent_GSA_BC408)[c(4,6,7,8)] = paste0("PDX_BC408.",colnames(Bivalent_GSA_BC408)[c(4,6,7,8)])

# BC1224
Bivalent_GSA_BC1224 = readxl::read_xlsx(file.path("output","bulk_ChIPseq","BC1224","ChIPreChIP","Enrichment_test_Bivalent_pathways_K4_K27.xlsx"))
Bivalent_GSA_BC1224$log10.qvalue = -log10(Bivalent_GSA_BC1224$`q-value`)
colnames(Bivalent_GSA_BC1224)[c(4,6,7,8)] = paste0("PDX_BC1224.",colnames(Bivalent_GSA_BC1224)[c(4,6,7,8)])


Bivalent_GSA_BC408 = Bivalent_GSA_BC408[order(Bivalent_GSA_BC408$Gene.Set),]
GSA_combined = cbind(Bivalent_GSA_BC976[match(Bivalent_GSA_BC408$Gene.Set, Bivalent_GSA_BC976$Gene.Set),1:3],
                     Bivalent_GSA_BC976[match(Bivalent_GSA_BC408$Gene.Set, Bivalent_GSA_BC976$Gene.Set),c(4,6,7,8)],
                     Bivalent_GSA_BC408[,c(4,6,7,8)],
                     Bivalent_GSA_BC1224[match(Bivalent_GSA_BC408$Gene.Set, Bivalent_GSA_BC1224$Gene.Set),c(4,6,7,8)])

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

png(file.path(output,"Intersection_PDX_bivalent_GSA.png"),width = 1200, height = 1200, res = 200)
gplots::venn(
    list("BC976" = GSA_combined_curated$Gene.Set[which(GSA_combined_curated$`PDX_BC976.q-value` < 0.1)],
         "BC408" = GSA_combined_curated$Gene.Set[which(GSA_combined_curated$`PDX_BC408.q-value` < 0.1)],
         "BC1224" = GSA_combined_curated$Gene.Set[which(GSA_combined_curated$`PDX_BC1224.q-value` < 0.1)]),
    show.plot = T)
dev.off()

SuperExactTest::MSET(list("BC976" = GSA_combined_curated$Gene.Set[which(GSA_combined_curated$`PDX_BC976.q-value` < 0.1)],
                          "BC408" = GSA_combined_curated$Gene.Set[which(GSA_combined_curated$`PDX_BC408.q-value` < 0.1)],
                          "BC1224" = GSA_combined_curated$Gene.Set[which(GSA_combined_curated$`PDX_BC1224.q-value` < 0.1)]),
                     n = length(unique(GSA_combined_curated$Gene.Set)),
                     lower.tail = FALSE)


WriteXLS::WriteXLS(GSA_combined_save, ExcelFileName = file.path(output,"combined_GSA_Bivalence_PDXs_unfiltered2.xlsx"))






