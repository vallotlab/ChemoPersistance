# BC408 Gene Based Intersection
library(here)
# Directories -------------------------------------------------------------
maindir= here()
resdir <- file.path(maindir,"output","scRNAseq","PDX_BC408")
source(file.path(maindir,"Scripts","global_var.R"))
source(file.path(maindir,"Scripts","functions.R"))

output <- file.path(resdir)

load(file.path(maindir,"output","scRNAseq","PDX_BC408",  "Unsupervised","RData","LogCounts.RData"))
LogCounts = LogCounts[,grep("persister",colnames(LogCounts))]
dim(LogCounts)
BinCounts = LogCounts
BinCounts[BinCounts>0] = 1

scRNA = readxl::read_xlsx(file.path(maindir,"output","scRNAseq","PDX_BC408", "Supervised",
                                    "Tables","Differential_analysis_Limma_logFC_Pers_vs_chemonaive_1.58.xlsx"),3)
scRNA$qval.persister = as.numeric(scRNA$qval.persister)
scRNA$log2FC.persister = as.numeric(scRNA$log2FC.persister)
scRNA$percent_expressed_BC408 = rowSums(BinCounts[match(scRNA$Symbol,rownames(BinCounts)),])/ncol(BinCounts)

persister_genes_BC408 = scRNA$Symbol[which(scRNA$log2FC.persister > log2(3) & 
                                               scRNA$qval.persister < 0.01)]
length(persister_genes_BC408)

# TF enrichment
load(file.path(maindir,"output","scRNAseq","PDX_BC408", "Unsupervised","TFmotif","list_TF_enrichment.RData") )
TF_BC408 = list_TF_enrichment$Genes_of_interest

intersect(persister_genes_BC408, TF_BC408$TF[1:100])

# Bivalency H3K27me3 - H3K4me3 / H3K4me3 - H3K27me3
K4_K27 = read.csv(file.path(maindir,"output","bulk_ChIPseq","BC408", "ChIPreChIP", "fisher_pvalue_K4_K27_peak_0.001_4.csv"))

# Group by gene, take most significant peak
K4_K27 = K4_K27 %>% tidyr::separate_rows(gene, sep =",") %>%
    group_by(gene) %>% summarise(
        odd_ratio = max(odd_ratio),,
        qvalue = min(qvalue),
        Significant = any(Significant)
    )
length(intersect(persister_genes_BC408, K4_K27$gene[which(K4_K27$Significant)]))
length(K4_K27$gene[which(K4_K27$Significant)])

### Build Table 

## TF
scRNA$TF_ChEA3_MeanRank = as.numeric(TF_BC408$Rank[match(scRNA$Symbol, TF_BC408$TF)])
scRNA$TF_ChEA3_Score = as.numeric(TF_BC408$Score[match(scRNA$Symbol, TF_BC408$TF)])

# Verif
View(scRNA[grep("THAP1", scRNA$Symbol),])

## Bivalency
scRNA$is_bivalent_K4_K27 = FALSE
scRNA$is_bivalent_K4_K27 = K4_K27$Significant[match(scRNA$Symbol, K4_K27$gene)]

# Verif
View(scRNA[grep("FOSL1|FOXQ1|CNKSR2", scRNA$Symbol),])

WriteXLS(c("scRNA"),
         ExcelFileName = file.path(output,"Multiomic_gene_based_table_PDX_BC408.xlsx"),
         perl = "perl", verbose = FALSE, row.names = FALSE,
         col.names = TRUE, AdjWidth = T, AutoFilter = TRUE,
         BoldHeaderRow = TRUE, na = "", FreezeRow = 1, FreezeCol = 1)

scRNA = readxl::read_xlsx(file.path(output,"Multiomic_gene_based_table_PDX_BC408.xlsx"))

scRNA_persister = scRNA[match(persister_genes_BC408, scRNA$Symbol),]

list = list(
    "Bivalent_K4_K27" = scRNA_persister$Symbol[which(scRNA_persister$is_bivalent_K4_K27 == "TRUE")],
    "Top_100_TF_ChEA3" = scRNA_persister$Symbol[which(scRNA_persister$TF_ChEA3_MeanRank < 100)]
)

sapply(list, length)
round(100*sapply(list, length) / length(persister_genes_BC408) ,2)

# K27 K4
length(intersect(list$Bivalent_K4_K27, list$Top_100_TF_ChEA3))
