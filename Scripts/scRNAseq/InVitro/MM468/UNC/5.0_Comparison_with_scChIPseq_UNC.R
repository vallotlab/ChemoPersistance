library(here)
library(ggpubr)
library(dplyr)

maindir = here()
source(file.path(maindir,"Scripts","global_var.R"))
QCdir <- file.path(maindir,"output","scRNAseq","MM468","QC")
resdir <- file.path(maindir,"output","scRNAseq","MM468","UNC")
resSUBdir <- file.path(resdir, paste0("ComparisonChIPseq")) ;if(!file.exists(resSUBdir)){dir.create(resSUBdir)}
resSUBdir_K27 <- file.path(resSUBdir, paste0("K27")) ;if(!file.exists(resSUBdir_K27)){dir.create(resSUBdir_K27)}
resSUBdir_K4 <- file.path(resSUBdir, paste0("K4")) ;if(!file.exists(resSUBdir_K4)){dir.create(resSUBdir_K4)}
rdatadir <- file.path(resdir, "Supervised","RData")

#Comparing scRNAseq datasets and ChIPseq datasets
load(file.path(rdatadir,"Supervised_res_object_edgeR.Rdata"))
res.scRNA = my.res

annot10k_K27 <- read.table(
  unz(file.path(
    maindir,"annotation","gencode.v34.transcripts10k_K27.zip"),
    "gencode.v34.transcripts10k_K27.bed"), sep="\t")[,-c(6,9)]
colnames(annot10k_K27) <- c("chr","start","end","transcripts","gene","log2FC_K27","qvalue_K27","K27_status")

# Keep only transcipts overlapping K27 peaks
annot10k_K27_ovlp = annot10k_K27[which(annot10k_K27$K27_status !="Not Overlapping"),]
annot10k_K27_byGene = annot10k_K27_ovlp %>% group_by(gene) %>%  slice_max(order_by = abs(log2FC_K27), n = 1) # Select only the transcript associated with top differential K27 log2FC 

sc_bulk_ChIP_DA = readxl::read_xlsx(file.path(maindir,"output","scChIPseq","ChromSCape_analyses",
                                              "MM468_H3K27me3_peaks","Diff_Analysis_Gene_Sets",
                                              "DA_grouped_Persister_vs_DSMO_with_bulk_10k.xlsx"), sheet = 3)

sc_bulk_ChIP_DA = sc_bulk_ChIP_DA %>% group_by(Gene) %>% slice_max(abs(log2FC.Persister))

# before grouping
annot10k_K27[grep("DAPK1",annot10k_K27$gene),]
# after grouping
annot10k_K27_byGene[grep("DAPK1",annot10k_K27_byGene$gene),]

log2FC_thresholds <- log2(c(2,3,4))

res.scRNA$log2FC_ChIP <- sc_bulk_ChIP_DA$log2FC.Persister[match(res.scRNA$Symbol,sc_bulk_ChIP_DA$Gene)]
res.scRNA$qval_K27 <- sc_bulk_ChIP_DA$qval.Persister[match(res.scRNA$Symbol,sc_bulk_ChIP_DA$Gene)]
res.scRNA$K27_status <- annot10k_K27_ovlp$K27_status[match(res.scRNA$Symbol,annot10k_K27_ovlp$gene)]


# Pie Chart K27 in overexpressed
subset_int = res.scRNA[res.scRNA$log2FC.UNC>log2(3) & res.scRNA$qval.UNC<0.01,]
subset_int$K27_status[is.na(subset_int$K27_status)] = "No K27 Peak"
subset_int$K27_status[which(subset_int$K27_status != "No K27 Peak")] = "Not differential"
subset_int$K27_status[which(subset_int$K27_status != "No K27 Peak" & subset_int$log2FC_ChIP < -log2(1.5) & subset_int$qval_K27 < 0.1)] = "Depleted FC < -1.5"
subset_int$K27_status[which(subset_int$K27_status != "No K27 Peak" & subset_int$log2FC_ChIP < -log2(2) & subset_int$qval_K27 < 0.1)] = "Depleted FC < -2"
subset_int$K27_status[which(subset_int$K27_status != "No K27 Peak" & subset_int$log2FC_ChIP < -log2(3) & subset_int$qval_K27 < 0.1)] = "Depleted FC < -3"
subset_int$K27_status[which(subset_int$K27_status != "No K27 Peak" & subset_int$log2FC_ChIP > log2(1.5) & subset_int$qval_K27 < 0.1)] = "Overexpressed - FC > 1.5"
subset_int$K27_status[which(subset_int$K27_status != "No K27 Peak" & subset_int$log2FC_ChIP > log2(2) & subset_int$qval_K27 < 0.1)] = "Overexpressed - FC > 2"
subset_int$K27_status[which(subset_int$K27_status != "No K27 Peak" & subset_int$log2FC_ChIP > log2(3) & subset_int$qval_K27 < 0.1)] = "Overexpressed - FC > 3"

tab = subset_int %>% dplyr::select(
  Symbol, log2FC.UNC, qval.UNC, log2FC_ChIP, qval_K27, K27_status)
WriteXLS(tab, file.path(resSUBdir_K27,"persister_K27_status.xls"), SheetNames = "Persister_K27")

vec = (table(subset_int$K27_status))[c(4,5,1,2,3)]

pdf(file.path(resSUBdir_K27,paste0("NumPeaksK27_scRNAseq_bulk_sc_UNC.pdf")),height=5,width=5)
pie(as.numeric(vec), labels = paste0(names(vec)," ; ",as.numeric(vec)),
    col = c("#999999ff","#D1CE68FD","#67CF85", "#33BD5A", "#00A64B"),
    cex=0.5)
dev.off()

# Pie Chart overexpressed in K27 
subset_int = annot10k_K27_byGene[annot10k_K27_byGene$log2FC_K27 < -log2(3) & annot10k_K27_byGene$qvalue_K27 < 0.1,]
dim(my.res)
subset_int = subset_int[which(subset_int$gene %in% my.res$Symbol),]
subset_int$scRNA_status = "Not Overexpressed"
subset_int$log2FC_scRNA = res.scRNA$log2FC.UNC[match(subset_int$gene, res.scRNA$Symbol)]
subset_int$qvalue_scRNA = res.scRNA$qval.UNC[match(subset_int$gene, res.scRNA$Symbol)]

subset_int$scRNA_status[which(subset_int$log2FC_scRNA > log2(1.5))] = "Overexpressed - FC > 1.5"
subset_int$scRNA_status[which(subset_int$log2FC_scRNA > log2(2))] = "Overexpressed - FC > 2"
subset_int$scRNA_status[which(subset_int$log2FC_scRNA > log2(3))] = "Overexpressed - FC > 3"

subset_int$scRNA_status[which(subset_int$log2FC_scRNA < -log2(2) )] = "Underexpressed - FC < -2"

table(subset_int$scRNA_status)

vec = (table(subset_int$scRNA_status))

pdf(file.path(resSUBdir_K27,paste0("NumGenes_scRNAseq_inK27_10k_",sum(vec),"genes.pdf")),height=5,width=5)
pie(as.numeric(vec), labels = paste0(names(vec)," ; ",as.numeric(vec)),
    col = c("#999999ff","#78b7c5ff","#56BCD6","#2E767A", "#BD4444"),
    cex=0.5)
text(0,0.975,paste0("Total = ",sum(vec)," genes"))
dev.off()
