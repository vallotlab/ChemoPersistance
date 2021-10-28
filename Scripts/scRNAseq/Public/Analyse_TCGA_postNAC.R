#### AnalysisTranscriptomic TNBC from Magbanua et al.
# I-SPY 1 TRIAL (Magbuana et al., Breast Cancer Research 2015) containing 18
# TNBC patients with paired samples before and post NAC. 

library(here)
source(file.path(here(),"Scripts","global_var.R"))
maindir = here()

log2FC_threshold = log2(3)
outdir <- file.path(maindir,"output","scRNAseq","TCGA", "PostNAC_cohort")
if(!dir.exists(outdir)) dir.create(outdir)

postNAC_cohort_dir = file.path("/media/pacome/LaCie/InstitutCurie/Documents/Data/results/Magbanua_el_al/GSE32603_RAW/")

mycolramp <- c("white",viridis(n=4))

# Retrieve "persister" genes
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
name_to_id = readxl::read_xlsx(file.path(postNAC_cohort_dir, "Name_to_ID.xlsx"), col_names = F)
name_to_id = name_to_id[,c(1,3,4)]
colnames(name_to_id) = c("ID","Gene", "Detail")

sample_data = read.table(file.path(postNAC_cohort_dir, "sample_metadata.tsv"),
                         sep = "\t", header = FALSE)
files = list.files(postNAC_cohort_dir, pattern = '.gpr')
list_gene_data = sapply(files, function(file){
    gene_data = read.table(file.path(postNAC_cohort_dir, file), sep = "\t", skip = 29)
    Gene = name_to_id$Gene[match(gene_data$V5, name_to_id$ID)]
    gene_data = gene_data[which(!is.na(Gene)),9, drop = F]
    gene_data$Gene = Gene[which(!is.na(Gene))]
    gene_data = as.data.frame(
        gene_data %>% dplyr::group_by(Gene) %>% summarise(V9 = mean(as.numeric(V9)))
    )
    rownames(gene_data) = gene_data$Gene
    gene_data = gene_data[,-1,drop=F]
    gene_data = as.matrix((gene_data))
    colnames(gene_data) = gsub(".gpr.gz","",basename(file))
    gene_data
})

sapply(list_gene_data, nrow)
genes = rownames(list_gene_data[["GSM808307.gpr.gz"]])
gene_mat = sapply(list_gene_data, function(x){
    x = x[genes,,drop=F]
    x
} )
rownames(gene_mat) = genes
colnames(gene_mat) = gsub(".gpr.gz","",colnames(gene_mat))
colnames(gene_mat) = sample_data$V2[match(colnames(gene_mat), sample_data$V1)]

write.csv(gene_mat, file.path(outdir, "gene_mat.csv"), row.names = TRUE)
gene_mat = read.csv(file.path(outdir, "gene_mat.csv"), row.names = 1, header = T)

# Number of sample before treatment
length(grep("T1", colnames(gene_mat)))
# Number of sample 24-96hours after inititaion of chemotherapy
length(grep("T2", colnames(gene_mat)))
# Number of sample at surgery
length(grep("TS", colnames(gene_mat)))

# Find triple negative breast cancer patient:
df = as.data.frame(t(rbind(gene_mat["ESR1",,drop=F], gene_mat["PGR",,drop=F], gene_mat["ERBB2",,drop=F])))
df = df %>% tidyr::gather("Gene", "median_expression")

# Plot TNBC status of the three receptors
pdf(file.path(outdir,"TNBC_status_genes.pdf"))
df %>% ggplot(aes(y = median_expression, x = Gene, fill = Gene)) + geom_violin() +
    scale_y_log10() + theme_classic()
dev.off()

### Find TNBC with PAM 50 
log2_gene_mat = log2(gene_mat)
write.table(log2_gene_mat, file.path(outdir, "log2_gene_mat.csv"), row.names = TRUE, sep = "\t", col.names = T)
log2_gene_mat = read.table(file.path(outdir, "log2_gene_mat.csv"),
                           row.names = 1, sep = "\t", header = T)

###
# input variables for the subtype prediction script
###
library(ctc)
library(heatmap.plus)

paramDir<- "output/scRNAseq/TCGA/PostNAC_cohort/bioclassifier_R" # the location of unchanging files such as the function library and main program
inputDir<- "output/scRNAseq/TCGA/PostNAC_cohort/"  # the location of the data matrix, and where output will be located

inputFile<- "log2_gene_mat.csv" # the input data matrix as a tab delimited text file
short <- "PostNAC_cohort" # short name that will be used for output files

calibrationParameters<- NA 	#the column of the "mediansPerDataset.txt" file to use for calibration; 
#NA will force centering within the test set & -1 will not do any 
#adjustment (when adjustment performed by used)

hasClinical<-FALSE 	#may include tumor size as second row, with 'T' as the gene name, 
#and encoded as binary (0 for size <= 2cm or 1 for size > 2cm)
#set this variable to FALSE if tumor size is not available

collapseMethod<-"mean" # can be mean or iqr (probe with max iqr is selected)
# typically, mean is preferred for long oligo and
# iqr is preferred for short oligo platforms


####
# run the assignment algorithm
####
source(paste(paramDir,"subtypePrediction_functions.R",sep="/"))
source(paste(paramDir,"subtypePrediction_distributed.R",sep="/"))

BreastCancer_subtype = read.table(file.path(outdir, "PostNAC_cohort_pam50scores.txt"), header = T)
table(BreastCancer_subtype$Call)

pdf(file.path(outdir, "repartition_subtypes.pdf"))
pie(table(BreastCancer_subtype$Call))
dev.off()

TNBC = which(BreastCancer_subtype$Call == "Basal")
gene_mat_TNBC = gene_mat[,TNBC]

# Number of sample before treatment
length(grep("T1", colnames(gene_mat_TNBC)))
# Number of sample 24-96hours after inititaion of chemotherapy
length(grep("T2", colnames(gene_mat_TNBC)))
# Number of sample at surgery
length(grep("TS", colnames(gene_mat_TNBC)))

TNBC_sample_data = data.frame("sample" = colnames(gene_mat_TNBC))
TNBC_sample_data$timepoint = gsub(".*_timePt|_patient_.*","", TNBC_sample_data$sample)
TNBC_sample_data$patient = gsub(".*_patient","patient", TNBC_sample_data$sample)
table(TNBC_sample_data$patient)

#####
patient_with_pair = names(table(TNBC_sample_data$patient)[
    which((table(TNBC_sample_data$patient)>1))])  # 5 TNBC patients with pair gene data

TNBC_gene_samples = intersect(colnames(gene_mat_TNBC),
                              TNBC_sample_data$sample[which(TNBC_sample_data$patient %in% patient_with_pair)]) # 23 samples
TNBC_sample_data = TNBC_sample_data[match(TNBC_gene_samples, TNBC_sample_data$sample),]
gene_mat_TNBC = gene_mat_TNBC[, match(TNBC_gene_samples, colnames(gene_mat_TNBC))] # 36 samples

# MM468
persister_genes_MM468 = intersect(persister_genes_MM468,rownames(gene_mat_TNBC))
gene_mat_TNBC_persister_MM468 = gene_mat_TNBC[match(persister_genes_MM468, rownames(gene_mat_TNBC)),]

plot_df = as.data.frame(gene_mat_TNBC_persister_MM468)
plot_df$Gene = rownames(plot_df)
plot_df = plot_df %>% tidyr::gather(sample, mean_expression, - Gene) 
plot_df = plot_df %>% mutate("type" = ifelse(grepl("T1",sample),"primary", "postNAC/surgery"))
plot_df$type = factor(plot_df$type, levels = c("primary", "postNAC/surgery"))

pdf(file.path(outdir, "persister_genes_all_MM468.pdf"))
p = plot_df %>%
    ggplot(aes(y = mean_expression, x = type)) +
    geom_jitter(height = 0, width = 0.15,cex=0.3,alpha=0.3) +
    ggtitle("All MM468 Persister Genes") +
    ggpubr::stat_compare_means(paired = F, method = "t.test", label.x.npc = 0.4) +
    scale_y_log10() + theme_classic() + ylab("Gene Expression")
print(ggpubr::add_summary(p, fun = "mean_sd")) 
dev.off()

pdf(file.path(outdir, "persister_genes_pair_bivalent_genes_MM468.pdf"))
p = plot_df %>% filter(Gene %in% DA_MM468$Symbol[which(DA_MM468$bivalent=="TRUE")]) %>%
    ggplot(aes(y = mean_expression, x = type)) +
    geom_jitter(height = 0, width = 0.15,cex=0.3,alpha=0.3) +
    ggtitle("Bivalent MM468 Persister Genes") +
    ggpubr::stat_compare_means(paired = F, method = "t.test", label.x.npc = 0.4) +
    scale_y_log10() + theme_classic() + ylab("Gene Expression")
print(ggpubr::add_summary(p, fun = "mean_sd"))
dev.off()

pdf(file.path(outdir, "persister_genes_pair_non_bivalent_genes_MM468.pdf"))
p = plot_df %>% filter(!Gene %in% DA_MM468$Symbol[which(DA_MM468$bivalent=="TRUE")]) %>%
    ggplot(aes(y = mean_expression, x = type)) +
    geom_jitter(height = 0, width = 0.15, cex=0.3,alpha=0.3) +
    ggtitle("Non Bivalent MM468 Persister Genes") +
    ggpubr::stat_compare_means(paired = F, method = "t.test", label.x.npc = 0.4) +
    scale_y_log10() + theme_classic()  + ylab("Gene Expression")
print(ggpubr::add_summary(p, fun = "mean_sd"))
dev.off()

pdf(file.path(outdir, "persister_genes_gene_by_gene_MM468.pdf"))
for(i in persister_genes_MM468){
    tryCatch({p = plot_df %>% filter(Gene == i) %>% ggplot(aes(y = mean_expression, x = type)) +
        geom_jitter(height = 0, width = 0.15) +
        ggtitle(paste0(i," - ", ifelse(DA_MM468$bivalent[which(DA_MM468$Symbol==i)], "bivalent", "non bivalent"))) +
        ggpubr::stat_compare_means(paired = F, method = "t.test", label.x.npc = 0.4) +
        scale_y_log10() + theme_classic()  + ylab("Gene Expression")
    print(ggpubr::add_summary(p, fun = "mean_sd"))},
    error = function(e){print(NULL)})
}
dev.off()

# PDX
persister_genes_PDX = intersect(persister_genes_PDX,rownames(gene_mat_TNBC))
gene_mat_TNBC_persister_PDX = gene_mat_TNBC[match(persister_genes_PDX, rownames(gene_mat_TNBC)),]

plot_df = as.data.frame(gene_mat_TNBC_persister_PDX)
plot_df$Gene = rownames(plot_df)
plot_df = plot_df %>% tidyr::gather(sample, mean_expression, - Gene) 
plot_df = plot_df %>% mutate("type" = ifelse(grepl("T1",sample),"primary", "postNAC/surgery"))
plot_df$type = factor(plot_df$type, levels = c("primary", "postNAC/surgery"))

pdf(file.path(outdir, "persister_genes_all_PDX.pdf"))
p = plot_df %>%
    ggplot(aes(y = mean_expression, x = type)) +
    geom_jitter(height = 0, width = 0.15,cex=0.3,alpha=0.3) +
    ggtitle("All PDX Persister Genes") +
    ggpubr::stat_compare_means(paired = F, method = "t.test", label.x.npc = 0.4) +
    scale_y_log10() + theme_classic() + ylab("Gene Expression")
print(ggpubr::add_summary(p, fun = "mean_sd")) 
dev.off()

pdf(file.path(outdir, "persister_genes_gene_by_gene_PDX.pdf"))
for(i in persister_genes_PDX){
    tryCatch({p = plot_df %>% filter(Gene == i) %>% ggplot(aes(y = mean_expression, x = type)) +
        geom_jitter(height = 0, width = 0.15) +
        ggtitle(paste0(i)) +
        ggpubr::stat_compare_means(paired = F, method = "t.test", label.x.npc = 0.4) +
        scale_y_log10() + theme_classic()  + ylab("Gene Expression")
    print(ggpubr::add_summary(p, fun = "mean_sd"))},
    error = function(e){print(NULL)})
}
dev.off()

# Genes overexpressed in 'persisting' chemonaives
genes_ov_persisting_chemonaive = c("S100A2", "LDHB")
gene_mat_TNBC_ov_persisting_chemonaive= gene_mat_TNBC[match(genes_ov_persisting_chemonaive, rownames(gene_mat_TNBC)),]

plot_df = as.data.frame(gene_mat_TNBC_ov_persisting_chemonaive)
plot_df$Gene = rownames(plot_df)
plot_df = plot_df %>% tidyr::gather(sample, mean_expression, - Gene) 
plot_df = plot_df %>% mutate("type" = ifelse(grepl("T1",sample),"primary", "postNAC/surgery"))
plot_df$type = factor(plot_df$type, levels = c("primary", "postNAC/surgery"))

# Genes overexpressed in 'resisting' persister
genes_ov_persisting_chemonaive = c("S100A2", "LDHB", "MAGED1", "TPM2", "AARD")
gene_mat_TNBC_ov_persisting_chemonaive= gene_mat_TNBC[match(genes_ov_persisting_chemonaive, rownames(gene_mat_TNBC)),]

plot_df = as.data.frame(gene_mat_TNBC_ov_persisting_chemonaive)
plot_df$Gene = rownames(plot_df)
plot_df = plot_df %>% tidyr::gather(sample, mean_expression, - Gene) 
plot_df = plot_df %>% mutate("type" = ifelse(grepl("T1",sample),"primary", "postNAC/surgery"))
plot_df$type = factor(plot_df$type, levels = c("primary", "postNAC/surgery"))

pdf(file.path(outdir, "ov_persisting_chemonaive_gene_by_gene_PDX.pdf"))
for(i in c("S100A2", "LDHB")){
    tryCatch({p = plot_df %>% filter(Gene == i) %>% ggplot(aes(y = mean_expression, x = type)) +
        geom_jitter(height = 0, width = 0.15) +
        ggtitle(paste0(i)) +
        ggpubr::stat_compare_means(paired = F, method = "t.test", label.x.npc = 0.4) +
        scale_y_log10() + theme_classic()  + ylab("Gene Expression")
    print(ggpubr::add_summary(p, fun = "mean_sd"))},
    error = function(e){print(NULL)})
}
dev.off()
