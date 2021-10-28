library(here)
source(file.path(here(),"Scripts","functions.R"))
source(file.path(here(),"Scripts","global_var.R"))

# Directories -------------------------------------------------------------
maindir = here()
resdir <- file.path(maindir,"output","scChIPseq","MM468","scH3K4me3"); if(!file.exists(resdir)){dir.create(resdir)}
unsupervised_dir  <- file.path(resdir,"Unsupervised","RData")
input_dir  <- file.path(maindir, "input","scChIPseq","MM468")

# Unsupervised directories
resdir <- file.path(resdir, "Supervised") ; if(!file.exists(resdir)){dir.create(resdir)}
resdir_tables = file.path(resdir,"Tables"); if(!dir.exists(resdir_tables)) dir.create(resdir_tables)
resdir_enrich <- file.path(resdir_tables,"Gene_set_analysis");if(!file.exists(resdir_enrich)){dir.create(resdir_enrich)}

dataset_name = "MM468_H3K4me3_10k_TSS"
ref_genome ="hg38"

metadata = data.frame(sample_id = c(
    "MM468_DMSO1_day0", "MM468_5FU1_day60"),
    sample_id_color = c("#afafafff", "#118675ff"))


# Differential analysis
qval.th = 0.1
logFC.th = 1
enrichment_qval = 0.1

# Supervised
scExp_Pers_vs_DMSO = scExp_cf
scExp_Pers_vs_DMSO$cell_cluster[which(scExp_Pers_vs_DMSO$sample_id == "MM468_5FU1_day60")] = "C2"
scExp_Pers_vs_DMSO$cell_cluster[which(scExp_Pers_vs_DMSO$sample_id != "MM468_5FU1_day60")] = "C1"
de_type  = "one_vs_rest"
scExp_Pers_vs_DMSO = differential_analysis_scExp(scExp_Pers_vs_DMSO,de_type = de_type)

### Gene Set Analysis
scExp_Pers_vs_DMSO = gene_set_enrichment_analysis_scExp(scExp_Pers_vs_DMSO,enrichment_qval = 0.1,
                                                        qval.th = qval.th, cdiff.th = logFC.th)

diff_Pers_vs_DMSO = scExp_Pers_vs_DMSO@metadata$diff$res
diff_Pers_vs_DMSO = diff_Pers_vs_DMSO[,c(1:4,10:14)]
colnames(diff_Pers_vs_DMSO) = gsub("C2","Pers_vs_DMSO",colnames(diff_Pers_vs_DMSO))

GSA_Pers_vs_DMSO = scExp_Pers_vs_DMSO@metadata$enr

qs::qsave(diff_Pers_vs_DMSO, file = file.path(resdir_tables,
                                                             paste0("diff_Pers_vs_DMSO_",qval.th,"_",logFC.th,"_",de_type,".qs"))) #save the data 
qs::qsave(GSA_Pers_vs_DMSO, file = file.path(resdir_tables,
                                                             paste0("GSA_Pers_vs_DMSO_",qval.th,"_",logFC.th,"_",de_type,".qs"))) #save the data 

annot_10k = read.table(gzfile(file.path(maindir,"annotation",
                                        "gencode.v34.annotation.transcriptTSS_10k.bed.gz")),
                       sep="\t", header = F)
colnames(annot_10k) = c("chr","start","end","transcripts","gene","strand")

diff_Pers_vs_DMSO$gene = 
    annot_10k$gene[match(rownames(diff_Pers_vs_DMSO),paste0(annot_10k$chr,"_",annot_10k$start,"_",annot_10k$end))]

under_res = diff_Pers_vs_DMSO %>% dplyr::filter(cdiff.Pers_vs_DMSO < -logFC.th & qval.Pers_vs_DMSO < qval.th) %>% 
    dplyr::arrange(qval.Pers_vs_DMSO) %>% dplyr::select(chr,start,end,,cdiff.Pers_vs_DMSO,qval.Pers_vs_DMSO,gene)
over_res = diff_Pers_vs_DMSO %>% dplyr::filter(cdiff.Pers_vs_DMSO > logFC.th & qval.Pers_vs_DMSO < qval.th) %>% 
    dplyr::arrange(qval.Pers_vs_DMSO) %>% dplyr::select(chr,start,end,cdiff.Pers_vs_DMSO,qval.Pers_vs_DMSO,gene)
under_res = data.frame
WriteXLS(c("over_res","under_res","diff_Pers_vs_DMSO"),
         ExcelFileName = file.path(resdir_tables,paste0("DA_Pers_vs_DMSO_logFC_",
                                                                   round(logFC.th,2),".xlsx")),
         SheetNames = c(
             paste0("Over_",round(logFC.th,2),"_",qval.th,"_n",
                    nrow(over_res)),paste0(
                        "Under_-",round(logFC.th,2),"_",qval.th,
                        "_n",nrow(under_res)),"All"),
         perl = "perl", verbose = FALSE, row.names = FALSE,
         col.names = TRUE, AdjWidth = T, AutoFilter = TRUE,
         BoldHeaderRow = TRUE, na = "", FreezeRow = 1, FreezeCol = 1)

Overexpressed = GSA_Pers_vs_DMSO$Overexpressed[[2]]
Underexpressed = GSA_Pers_vs_DMSO$Underexpressed[[2]]

WriteXLS(
    c("Overexpressed", "Underexpressed"),
    ExcelFileName = file.path(resdir_tables,
                              paste0("Enrichment_test_Pers_vs_DMSO_logFC",
                                     round(logFC.th,2),"_TSS.xlsx")), 
    SheetNames = c( paste0("Overexp_in_Pers"), paste0("Underexp_in_Pers") ),
    perl = "perl", verbose = FALSE, row.names = FALSE, col.names = TRUE,
    AdjWidth = TRUE, AutoFilter = TRUE, BoldHeaderRow = TRUE, na = "",
    FreezeRow = 1, FreezeCol = 1)
