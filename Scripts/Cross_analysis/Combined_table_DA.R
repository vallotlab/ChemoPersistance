
library(here)
source(file.path(here(),"Scripts","functions.R"))
source(file.path(here(),"Scripts","global_var.R"))

maindir = here()

## MM468
scRNA_pers = readxl::read_xlsx(file.path("output","scRNAseq","MM468","Persister","Supervised","Tables","Differential_analysis_Limma_logFC_1.58.xlsx"),sheet = 3)
scRNA_UNC = readxl::read_xlsx(file.path("output","scRNAseq","MM468","UNC","Supervised","Tables","Differential_analysis_Limma_logFC_1.58.xlsx"),sheet = 3)
scChIP_K27 = readxl::read_xlsx(file.path("output","scChIPseq","ChromSCape_analyses","MM468_H3K27me3_peaks","Diff_Analysis_Gene_Sets","DA_grouped_Persister_vs_DSMO_with_bulk_10k.xlsx"),sheet = 3)
bivalent_df = read.csv(file.path("output","bulk_ChIPseq","MM468","ChIPreChIP","fisher_pvalue_K27_K4_peak_0.001_0.15.csv"))
colnames(bivalent_df)[c(3,5,6)] = c("OddRatio_bivalent","q.value_bivalent","Significantly_bivalent")
bivalent_df = bivalent_df %>% tidyr::separate_rows(gene)
bivalent_df = bivalent_df %>% group_by(gene) %>%
    slice_min(q.value_bivalent)

scChIP_K27_byGene = scChIP_K27 %>% group_by(Gene) %>%
    slice_max(abs(log2FC.Persister))

all_genes = (union(scRNA_UNC$Symbol,union(scChIP_K27_byGene$Gene,scRNA_pers$Symbol)))
length(all_genes)

combined_table = tibble(Gene = all_genes)
combined_table = left_join(combined_table, scRNA_pers[,c(1,4,7)], by = c("Gene"="Symbol"))
combined_table = left_join(combined_table, scRNA_UNC[,c(1,4,7)], by = c("Gene"="Symbol"))
combined_table = left_join(combined_table, scChIP_K27_byGene[,c(5,7,8)], by = c("Gene"="Gene")) %>% unique
combined_table = left_join(combined_table, bivalent_df[,c(3,5,6,7)], by = c("Gene"="gene")) %>% unique
colnames(combined_table) = c("Gene", "log2FC.scRNA.Persister","qval.scRNA.Persister",
                             "log2FC.scRNA.UNC","qval.scRNA.UNC",
                             "log2FC.scChIP.H3K27me3","qval.scChIP.H3K27me3",
                             "OddRatio_bivalent", "q.value_bivalent", "Significantly_bivalent")
WriteXLS(combined_table, ExcelFileName = file.path("output","combined_DA_MM468.xlsx"))

## MM468
scRNA_PDX = readxl::read_xlsx(file.path("output","scRNAseq","PDX","Supervised","Tables","Differential_analysis_Limma_logFC_1.58.xlsx"),sheet = 3)
all_genes = union(scRNA_PDX$Symbol,bivalent_df$gene)
length(all_genes)

combined_table = tibble(Gene = all_genes)
combined_table = left_join(combined_table, scRNA_PDX[,c(1,4,7)], by = c("Gene"="Symbol"))
combined_table = left_join(combined_table, bivalent_df[,c(3,5,6,7)], by = c("Gene"="gene")) %>% unique
colnames(combined_table) = c("Gene", "log2FC.scRNA.PDX","qval.scRNA.PDX",
                             "OddRatio_bivalent", "q.value_bivalent", "Significantly_bivalent")
WriteXLS(combined_table, ExcelFileName = file.path("output","combined_DA_PDX.xlsx"))

over_PDX = combined_table %>% filter(log2FC.scRNA.PDX > log2(3) & qval.scRNA.PDX < 0.01)
table(over_PDX$Significantly_bivalent)
over_PDX$Gene[which(over_PDX$Significantly_bivalent)]

scRNA_pers$qval.C2_pers = as.numeric(scRNA_pers$qval.C2_pers)
over_MM468  = scRNA_pers %>% filter(log2FC.C2_pers > log2(3) & qval.C2_pers < 0.01)
