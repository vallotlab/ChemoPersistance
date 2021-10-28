# MM468 Gene Based Intersection
library(here)
# Directories -------------------------------------------------------------
maindir= here()
resdir <- file.path(maindir,"output","scRNAseq","MM468")
source(file.path(maindir,"Scripts","global_var.R"))
source(file.path(maindir,"Scripts","functions.R"))

output <- file.path(resdir)

load(file.path(maindir,"output","scRNAseq","MM468", "Persister", "Unsupervised","RData","LogCounts.RData"))
LogCounts = LogCounts[,grep("5FU1_day33",colnames(LogCounts))]
dim(LogCounts)
BinCounts = LogCounts
BinCounts[BinCounts>0] = 1

scRNA = readxl::read_xlsx(file.path(maindir,"output","scRNAseq","MM468", "Persister","Supervised",
                                                 "Tables","Differential_analysis_Limma_logFC_1.58.xlsx"),3)
scRNA$qval.C2_pers = as.numeric(scRNA$qval.C2_pers)
scRNA$log2FC.C2_pers = as.numeric(scRNA$log2FC.C2_pers)
scRNA$percent_expressed_MM468 = rowSums(BinCounts[match(scRNA$Symbol,rownames(BinCounts)),])/ncol(BinCounts)
persister_genes_MM468 = scRNA$Symbol[which(scRNA$log2FC.C2_pers > log2(3) & 
                                                            scRNA$qval.C2_pers < 0.01)]
length(persister_genes_MM468)

# TF enrichment
load(file.path(maindir,"output","scRNAseq","MM468", "Persister","Unsupervised","TFmotif","list_TF_enrichment.RData") )
TF_MM468 = list_TF_enrichment$Genes_of_interest

intersect(persister_genes_MM468, TF_MM468$TF[1:100])

# TF enrichment
library(AUCell)
library(SCENIC)
regulon_SCENIC = readRDS("/media/pacome/Depic_bioinfo_1/InstitutCurie/Documents/Data/results/SCENIC_ChemoPersistence/MM468/int/3.4_regulonAUC.Rds")
cellInfo = readRDS("/media/pacome/Depic_bioinfo_1/InstitutCurie/Documents/Data/results/SCENIC_ChemoPersistence/MM468/int/cellInfo.Rds")
rss_SCENIC <- calcRSS(AUC=getAUC(regulon_SCENIC), cellAnnotation=cellInfo[, "sample_id"])

rss_SCENIC_extended = rss_SCENIC[grep("_extended", rownames(rss_SCENIC)),]
rownames(rss_SCENIC_extended) = gsub("_extended","", rownames(rss_SCENIC_extended))
rss_SCENIC = rss_SCENIC[grep("_extended", rownames(rss_SCENIC), invert = TRUE),]

rss_SCENIC = rbind(rss_SCENIC, rss_SCENIC_extended[which(!rownames(rss_SCENIC_extended) %in% rownames(rss_SCENIC)),])

rssPlot <- plotRSS(rss_SCENIC)
TF_SCENIC = rssPlot$df[which(rssPlot$df$cellType == "MM468_5FU6_day33"),]


# Gene-based H3K27me3 differential analysis (persister vs chemonaive)
ChIPseq_gene = readxl::read_xlsx(file.path(maindir,"output","scChIPseq", "ChromSCape_analyses","MM468_H3K27me3_peaks",
                                                 "Diff_Analysis_Gene_Sets","DA_grouped_Persister_vs_DSMO_with_bulk_TSS.xlsx"),3)
ChIPseq_gene = ChIPseq_gene %>% group_by(Gene) %>% slice_max(abs(log2FC.Persister))


# Peak-based H3K27me3 differential analysis (persister vs chemonaive)
ChIPseq_peak = readxl::read_xlsx(file.path(maindir,"output","scChIPseq", "ChromSCape_analyses","MM468_H3K27me3_peaks",
                                      "Diff_Analysis_Gene_Sets","DA_grouped_Persister_vs_DSMO_with_bulk_peaks.xlsx"),3)
ChIPseq_peak = ChIPseq_peak %>% tidyr::separate_rows(Gene, sep =",") %>%
    group_by(Gene) %>% slice_max(abs(log2FC.Persister))

# Bivalency H3K27me3 - H3K4me3 / H3K4me3 - H3K27me3
K4_K27 = read.csv(file.path(maindir,"output","bulk_ChIPseq","MM468", "ChIPreChIP", "fisher_pvalue_K4_K27_peak_0.001_4.csv"))
K27_K4 = read.csv(file.path(maindir,"output","bulk_ChIPseq","MM468", "ChIPreChIP", "fisher_pvalue_K27_K4_peak_0.001_0.15.csv"))

# Group by gene, take most significant peak
K4_K27 = K4_K27 %>% tidyr::separate_rows(gene, sep =",") %>%
    group_by(gene) %>% summarise(
        odd_ratio = max(odd_ratio),,
        qvalue = min(qvalue),
        Significant = any(Significant)
    )
length(intersect(persister_genes_MM468, K4_K27$gene[which(K4_K27$Significant)]))

# Group by gene, take most significant peak
K27_K4 = K27_K4 %>% tidyr::separate_rows(gene, sep =",") %>%
    group_by(gene) %>% summarise(
        odd_ratio = max(odd_ratio),,
        qvalue = min(qvalue),
        Significant = any(Significant)
    )
length(intersect(persister_genes_MM468, K27_K4$gene[which(K27_K4$Significant)]))

# Differential Expression upon K27 perturbation
scRNA_UNC = readxl::read_xlsx(file.path(maindir,"output","scRNAseq","MM468", "UNC","Supervised",
                                    "Tables","Differential_analysis_Limma_logFC_1.58.xlsx"),3)

### Build Table 

## TF
scRNA$TF_ChEA3_MeanRank = as.numeric(TF_MM468$Rank[match(scRNA$Symbol, TF_MM468$TF)])
scRNA$TF_ChEA3_Score = as.numeric(TF_MM468$Score[match(scRNA$Symbol, TF_MM468$TF)])
scRNA$TF_SCENIC_RSS = TF_SCENIC$RSS[match(scRNA$Symbol, TF_SCENIC$Topic)]
scRNA$TF_SCENIC_Z = TF_SCENIC$Z[match(scRNA$Symbol, TF_SCENIC$Topic )]

# Verif
View(scRNA[grep("FOSL1", scRNA$Symbol),])

## ChIPseq differential gene
scRNA$logFC_K27_gene = ChIPseq_gene$log2FC.Persister[match(scRNA$Symbol, ChIPseq_gene$Gene)]
scRNA$qval_K27_gene = ChIPseq_gene$qval.Persister[match(scRNA$Symbol, ChIPseq_gene$Gene)]

scRNA$logFC_K27_peak = ChIPseq_peak$log2FC.Persister[match(scRNA$Symbol, ChIPseq_peak$Gene)]
scRNA$qval_K27_peak = ChIPseq_peak$qval.Persister[match(scRNA$Symbol, ChIPseq_peak$Gene)]

# Verif
View(scRNA %>% filter(Symbol %in% persister_genes_MM468, scRNA$logFC_K27_peak < -1 | scRNA$logFC_K27_gene < -1))
scRNA %>% filter(Symbol %in% persister_genes_MM468, (scRNA$logFC_K27_peak < -1  & scRNA$qval_K27_peak < 0.1) | 
                     (scRNA$logFC_K27_gene    & scRNA$qval_K27_gene < 0.1))

scRNA$Depleted_in_H3K27me3 = FALSE
scRNA$Depleted_in_H3K27me3[which((scRNA$logFC_K27_peak < -1  & scRNA$qval_K27_peak < 0.1) | 
                (scRNA$logFC_K27_gene    & scRNA$qval_K27_gene < 0.1))] = TRUE

## Bivalency
scRNA$is_bivalent_K4_K27 = FALSE
scRNA$is_bivalent_K4_K27 = K4_K27$Significant[match(scRNA$Symbol, K4_K27$gene)]
scRNA$odd_ratio_K4_K27 = K4_K27$odd_ratio[match(scRNA$Symbol, K4_K27$gene)]
scRNA$qvalue_K4_K27 = K4_K27$qvalue[match(scRNA$Symbol, K4_K27$gene)]

scRNA$is_bivalent_K27_K4 = FALSE
scRNA$is_bivalent_K27_K4 = K27_K4$Significant[match(scRNA$Symbol, K27_K4$gene)]
scRNA$odd_ratio_K27_K4 = K27_K4$odd_ratio[match(scRNA$Symbol, K27_K4$gene)]
scRNA$qvalue_K27_K4 = K27_K4$qvalue[match(scRNA$Symbol, K27_K4$gene)]

# Verif
View(scRNA[grep("FOSL1|FOXQ1|CNKSR2", scRNA$Symbol),])

# Expression upon perturbation
scRNA$UNC_logFC = scRNA_UNC$log2FC.UNC[match(scRNA$Symbol, scRNA_UNC$Symbol)]
scRNA$UNC_qvalue = scRNA_UNC$qval.UNC[match(scRNA$Symbol, scRNA_UNC$Symbol)]

# Verif
length(intersect(persister_genes_MM468, scRNA$Symbol[which(scRNA$UNC_logFC > log2(3) & scRNA$UNC_qvalue < 0.01)]))

WriteXLS(c("scRNA"),
         ExcelFileName = file.path(output,"Multiomic_gene_based_table_MM468.xlsx"),
         perl = "perl", verbose = FALSE, row.names = FALSE,
         col.names = TRUE, AdjWidth = T, AutoFilter = TRUE,
         BoldHeaderRow = TRUE, na = "", FreezeRow = 1, FreezeCol = 1)

scRNA = readxl::read_xlsx(file.path(output,"Multiomic_gene_based_table_MM468.xlsx"))
scRNA_persister = scRNA[match(persister_genes_MM468, scRNA$Symbol),]

list = list(
    "Depleted_K27" = scRNA_persister$Symbol[which( (scRNA_persister$logFC_K27_gene < -1 & scRNA_persister$qval_K27_gene < 0.1) |  
                                                       (scRNA_persister$logFC_K27_peak < -1 & scRNA_persister$qval_K27_peak < 0.1) )],
    "Bivalent_K27_K4" = scRNA_persister$Symbol[which(scRNA_persister$is_bivalent_K27_K4 == "TRUE")],
    "Bivalent_K4_K27" = scRNA_persister$Symbol[which(scRNA_persister$is_bivalent_K4_K27 == "TRUE")],
    "Union_Bivalent" = scRNA_persister$Symbol[which(scRNA_persister$is_bivalent_K4_K27 == "TRUE" | 
                                                        scRNA_persister$is_bivalent_K27_K4 == "TRUE")],
    "Top_100_TF_ChEA3" = scRNA_persister$Symbol[which(scRNA_persister$TF_ChEA3_MeanRank < 100)],
    "Top_SCENIC" = scRNA_persister$Symbol[which(scRNA_persister$TF_SCENIC_RSS > 0.5 & scRNA_persister$TF_SCENIC_Z > 1.5)],
    "ReExpressed_UNC" = scRNA_persister$Symbol[which(scRNA_persister$UNC_logFC > log2(3) & scRNA_persister$UNC_qvalue < 0.01)]
)

gplots::venn(
    list[c(1:4,6)]
)
sapply(list, length)
round(100*sapply(list, length) / 168 ,2)


# Depletion
length(intersect(list$Depleted_K27, list$Bivalent_K27_K4))
length(intersect(list$Depleted_K27, list$Bivalent_K4_K27))
length(intersect(list$Depleted_K27, list$Union_Bivalent))
length(intersect(list$Depleted_K27, list$Top_100_TF_ChEA3))
length(intersect(list$Depleted_K27, list$Top_SCENIC))
length(intersect(list$Depleted_K27, list$ReExpressed_UNC))

# K27 K4
length(intersect(list$Bivalent_K27_K4, list$Bivalent_K4_K27))
length(intersect(list$Bivalent_K27_K4, list$Top_100_TF_ChEA3))
length(intersect(list$Bivalent_K27_K4, list$Top_SCENIC))
length(intersect(list$Bivalent_K27_K4, list$ReExpressed_UNC))


# Union 
length(intersect(list$Union_Bivalent, list$Bivalent_K4_K27))
length(intersect(list$Union_Bivalent, list$Top_100_TF_ChEA3))
length(intersect(list$Union_Bivalent, list$Top_SCENIC))
length(intersect(list$Union_Bivalent, list$ReExpressed_UNC))


# K4 K27
length(intersect(list$Bivalent_K4_K27, list$Bivalent_))
length(intersect(list$Bivalent_K27_K4, list$Top_100_TF_ChEA3))
length(intersect(list$Bivalent_K27_K4, list$Top_SCENIC))
length(intersect(list$Bivalent_K27_K4, list$ReExpressed_UNC))

length(intersect(list$Top_100_TF_ChEA3, list$ReExpressed_UNC))


