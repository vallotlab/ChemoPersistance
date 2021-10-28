#### Transcriptomic Analysis TNBC from TCGA
# Analysis of TCGA Metastatic Breast Cancer Cohort (https://mbcproject.org/)
# containing 5 TNBC patients with paired bulk RNA-seq samples of primary and
# metastatic tumors.
library(here)
source(file.path(here(),"Scripts","global_var.R"))
maindir = here()

log2FC_threshold = log2(3)
outdir <- file.path(maindir,"output","scRNAseq","TCGA", "Metastatic_cohort")
if(!dir.exists(outdir)) dir.create(outdir)

metastatic_cohort_dir = file.path(maindir,"input","TCGA", "brca_mbcproject_wagle_2017")

mycolramp <- c("white",viridis(n=4))

# Retrieve MM468 "persister" genes
DA_MM468 = readxl::read_xlsx(file.path("output","bulk_ChIPseq","MM468","ChIPreChIP","Persister_K27_bivalence.xlsx"))
persister_genes_MM468 = DA_MM468$Symbol

# Retrieve PDX "persister" genes
DA_PDX = readxl::read_xlsx(file.path("output","scRNAseq","MM468_PDX","Combined_persister_vs_chemonaive_DA.xlsx"))
persister_genes_PDX = DA_PDX$Genes[which(DA_PDX$log2FC.HBCx95 > log2(3) & DA_PDX$log2FC.BC408 > log2(3) &
                                                 DA_PDX$qval.HBCx95 < 0.01 & DA_PDX$qval.MM468 < 0.01 )]


# Find TNBC patients
# Load patient data
patient_data = readxl::read_xlsx(file.path(metastatic_cohort_dir, "data_clinical_patient.xlsx"))
TNBC_patient_data = patient_data %>% dplyr::filter(TRIPLE_NEG_STATUS == "YES")

# Load sample data
sample_data = readxl::read_xlsx(file.path(metastatic_cohort_dir, "data_clinical_sample.xlsx"))
TNBC_sample_data = sample_data %>% dplyr::filter(PATIENT_ID %in% TNBC_patient_data$PATIENT_ID)
TNBC_sample_data$SAMPLE_ID = gsub("-",".",TNBC_sample_data$SAMPLE_ID)

# Load gene expression datamatrix
gene_mat = read.table(file.path(metastatic_cohort_dir, "data_RNA_Seq_v2_expression_median.txt"),
                      sep = "\t",header = T)
# Deduplicate gene names / remove empty gene names
gene_mat = gene_mat[,-1]
gene_mat = unique(gene_mat)
gene_mat = gene_mat[-which(gene_mat$Hugo_Symbol == "" | is.na(gene_mat$Hugo_Symbol)),]
gene_mat = gene_mat[-which(duplicated(gene_mat$Hugo_Symbol)),]
rownames(gene_mat) = gene_mat$Hugo_Symbol
gene_mat = as.matrix(gene_mat[,-1])
dim(gene_mat)

patient_with_pair = names(table(TNBC_sample_data$PATIENT_ID)[
    which((table(TNBC_sample_data$PATIENT_ID)>1))])  # 5 TNBC patients with pair gene data
TNBC_gene_samples = intersect(colnames(gene_mat),
                              TNBC_sample_data$SAMPLE_ID[which(TNBC_sample_data$PATIENT_ID %in% patient_with_pair)]) # 23 samples
TNBC_sample_data = TNBC_sample_data[match(TNBC_gene_samples, TNBC_sample_data$SAMPLE_ID),]
gene_mat_TNBC = gene_mat[, match(TNBC_gene_samples, colnames(gene_mat))] # 10 samples

colnames(gene_mat_TNBC) = TNBC_sample_data$SAMPLE_TIMEPOINT
colnames(gene_mat_TNBC) = gsub("T2B|T3","T2", colnames(gene_mat_TNBC))
gene_mat_TNBC

## MM468
gene_mat_TNBC_persister_MM468 = gene_mat_TNBC[match(persister_genes_MM468, rownames(gene_mat_TNBC)),]

plot_df = as.data.frame(gene_mat_TNBC_persister_MM468)
plot_df$Gene = rownames(plot_df)
plot_df = plot_df %>% tidyr::gather(sample, mean_expression, - Gene) 
plot_df = plot_df %>% mutate("type" = ifelse(grepl("T1",sample),"residual", "metastatic"))
plot_df$type = factor(plot_df$type, levels = c("residual", "metastatic"))
plot_df %>% ggplot() + geom_violin(aes(y = mean_expression, x = type)) +
    geom_jitter(aes(y = mean_expression, x = type)) + scale_y_log10() 


pdf(file.path(outdir, "persister_genes_MM468_metastatic_pair_all.pdf"))
p = plot_df %>%
    ggplot(aes(y = mean_expression, x = type)) +
    geom_jitter(height = 0, width = 0.15, cex=0.3,alpha=0.3) +
    ggtitle("All Persister Genes") +
    ggpubr::stat_compare_means(paired = F, method = "t.test", label.x.npc = 0.4) +
    scale_y_log10() + theme_classic() + ylab("Gene Expression")
print(ggpubr::add_summary(p, fun = "mean_sd"))
dev.off()

pdf(file.path(outdir, "persister_genes_MM468_metastatic_pair_bivalent_genes.pdf"))
p = plot_df %>% filter(Gene %in% DA_MM468$Symbol[which(DA_MM468$bivalent=="TRUE")]) %>%
    ggplot(aes(y = mean_expression, x = type)) +
    geom_jitter(height = 0, width = 0.15, cex=0.3,alpha=0.3) +
    ggtitle("Bivalent Persister Genes") +
    ggpubr::stat_compare_means(paired = F, method = "t.test", label.x.npc = 0.4) +
    scale_y_log10() + theme_classic() + ylab("Gene Expression")
print(ggpubr::add_summary(p, fun = "mean_sd"))
dev.off()

pdf(file.path(outdir, "persister_genes_MM468_metastatic_pair_non_bivalent_genes.pdf"))
p = plot_df %>% filter(!Gene %in% DA_MM468$Symbol[which(DA_MM468$bivalent=="TRUE")]) %>%
    ggplot(aes(y = mean_expression, x = type)) +
    geom_jitter(height = 0, width = 0.15, cex=0.3,alpha=0.3) +
    ggtitle("Non Bivalent Persister Genes") +
    ggpubr::stat_compare_means(paired = F, method = "t.test", label.x.npc = 0.4) +
    scale_y_log10() + theme_classic() + ylab("Gene Expression")
print(ggpubr::add_summary(p, fun = "mean_sd"))
dev.off()


pdf(file.path(outdir, "persister_genes_MM468_metastatic_pair_all_genes.pdf"))
for(i in persister_genes){
    tryCatch({
        p = plot_df %>% filter(Gene == i)  %>% 
            ggplot(aes(y = mean_expression, x = type)) +
        geom_jitter(height = 0, width = 0.15) +
        ggtitle(paste0(i," - ", ifelse(DA_MM468$bivalent[which(DA_MM468$Symbol==i)], "bivalent", "non bivalent"))) +
        ggpubr::stat_compare_means(method = "t.test", label.x.npc = 0.4) +
        scale_y_log10() + theme_classic() + ylab("Gene Expression")
    print(ggpubr::add_summary(p, fun = "mean_sd"))},
    error = function(e){print(NULL)})
}
dev.off()

## PDX
gene_mat_TNBC_persister_PDX = gene_mat_TNBC[match(persister_genes_PDX, rownames(gene_mat_TNBC)),]

plot_df = as.data.frame(gene_mat_TNBC_persister_PDX)
plot_df$Gene = rownames(plot_df)
plot_df = plot_df %>% tidyr::gather(sample, mean_expression, - Gene) 
plot_df = plot_df %>% mutate("type" = ifelse(grepl("T1",sample),"residual", "metastatic"))
plot_df$type = factor(plot_df$type, levels = c("residual", "metastatic"))
plot_df %>% ggplot() + geom_violin(aes(y = mean_expression, x = type)) +
    geom_jitter(aes(y = mean_expression, x = type)) + scale_y_log10()  + ylab("Gene Expression")


pdf(file.path(outdir, "persister_genes_PDX_metastatic_pair_all.pdf"))
p = plot_df %>%
    ggplot(aes(y = mean_expression, x = type)) +
    geom_jitter(height = 0, width = 0.15, cex=0.3,alpha=0.3) +
    ggtitle("All PDX Persister Genes") +
    ggpubr::stat_compare_means(paired = F, method = "t.test", label.x.npc = 0.4) +
    scale_y_log10() + theme_classic() + ylab("Gene Expression")
print(ggpubr::add_summary(p, fun = "mean_sd"))
dev.off()


pdf(file.path(outdir, "persister_genes_PDX_metastatic_pair_all_genes.pdf"))
for(i in persister_genes_PDX){
    tryCatch({
        p = plot_df %>% filter(Gene == i)  %>% 
            ggplot(aes(y = mean_expression, x = type)) +
            geom_jitter(height = 0, width = 0.15) +
            ggtitle(i) +
            ggpubr::stat_compare_means(method = "t.test", label.x.npc = 0.4) +
            scale_y_log10() + theme_classic() + ylab("Gene Expression")
        print(ggpubr::add_summary(p, fun = "mean_sd"))},
        error = function(e){print(NULL)})
}
dev.off()

## persisting chemonaive vs non persisting chemonaives
gene_mat_TNBC_persister_p_vs_non_p_chemonaive = gene_mat_TNBC[match(c("S100A2","LDHB"), rownames(gene_mat_TNBC)),]

plot_df = as.data.frame(gene_mat_TNBC_persister_p_vs_non_p_chemonaive)
plot_df$Gene = rownames(plot_df)
plot_df = plot_df %>% tidyr::gather(sample, mean_expression, - Gene) 
plot_df = plot_df %>% mutate("type" = ifelse(grepl("T1",sample),"residual", "metastatic"))
plot_df$type = factor(plot_df$type, levels = c("residual", "metastatic"))
plot_df %>% ggplot() + geom_violin(aes(y = mean_expression, x = type)) +
    geom_jitter(aes(y = mean_expression, x = type)) + scale_y_log10()  + ylab("Gene Expression")


pdf(file.path(outdir, "persister_genes_persister_p_vs_non_p_chemonaive_metastatic_pair_all_genes.pdf"))
for(i in c("S100A2","LDHB")){
    tryCatch({
        p = plot_df %>% filter(Gene == i)  %>% 
            ggplot(aes(y = mean_expression, x = type)) +
            geom_jitter(height = 0, width = 0.15) +
            ggtitle(i) +
            ggpubr::stat_compare_means(method = "t.test", label.x.npc = 0.4) +
            scale_y_log10() + theme_classic() + ylab("Gene Expression")
        print(ggpubr::add_summary(p, fun = "mean_sd"))},
        error = function(e){print(NULL)})
}
dev.off()

