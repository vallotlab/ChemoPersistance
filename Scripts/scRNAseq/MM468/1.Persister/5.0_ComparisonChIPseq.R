library(here)
library(ggpubr)
library(dplyr)

maindir = here()
source(file.path(maindir,"Scripts","global_var.R"))
QCdir <- file.path(maindir,"output","scRNAseq","MM468","QC")
resdir <- file.path(maindir,"output","scRNAseq","MM468","Persister")
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

annot10k_K27_ovlp = annot10k_K27[which(annot10k_K27$K27_status != "Not Overlapping"),]

sc_bulk_ChIP_DA = readxl::read_xlsx(file.path(maindir,"output","scChIPseq","ChromSCape_analyses",
                                      "MM468_H3K27me3_peaks","Diff_Analysis_Gene_Sets",
                                      "DA_grouped_Persister_vs_DSMO_with_bulk_TSS.xlsx"), sheet = 3)

sc_bulk_ChIP_DA = sc_bulk_ChIP_DA %>% group_by(Gene) %>% slice_max(abs(log2FC.Persister))

log2FC_thresholds <- log2(3)

res.scRNA$log2FC_ChIP <- sc_bulk_ChIP_DA$log2FC.Persister[match(res.scRNA$Symbol,sc_bulk_ChIP_DA$Gene)]
res.scRNA$qval_K27 <- sc_bulk_ChIP_DA$qval.Persister[match(res.scRNA$Symbol,sc_bulk_ChIP_DA$Gene)]
res.scRNA$K27_status <- annot10k_K27_ovlp$K27_status[match(res.scRNA$Symbol,annot10k_K27_ovlp$gene)]

# Quantile of expression
res.scRNA$decile <- ntile(res.scRNA$log2FC.C2_pers, 10)  
figure <- res.scRNA %>% group_by(decile) %>% 
  summarise(ChIP=length(which(log2FC_ChIP<(-1) & qval_K27<0.1)),
            ChIP_res=length(which(log2FC_ChIP<(-1) & qval_K27<0.1)),
            expression=max(log2FC.C2_pers))

pdf(file.path(resSUBdir_K27,paste0("ExpressionDecileforSignificantDepletedPeaks_logFC1_10k.pdf")),height=5,width=5)
print(
  ggplot(data=figure, aes(x=decile , y=ChIP, fill=decile)) +
    geom_rect(data=figure,aes(xmin=decile-0.5,xmax=decile+0.5,ymin=-Inf,ymax=Inf),fill=alpha(colorRampPalette(c("royalblue","white","indianred1"))(10),0.8)) + 
    geom_point(data=figure, aes(x=decile, y=ChIP)) + theme_classic() + theme(legend.position = "none") + xlab("log2FC persister vs DMSO decile") +
    ylab("number of significantly depleted peaks")  + geom_text(aes(x=decile,y=40,label=round(expression,2)))
)
dev.off()

# Pie Chart K27 in overexpressed
subset_int = res.scRNA[res.scRNA$log2FC.C2_pers>log2(3) & res.scRNA$qval.C2_pers<0.01,]
subset_int$K27_status[is.na(subset_int$K27_status)] = "No K27 Peak"
subset_int$K27_status[which(subset_int$K27_status != "No K27 Peak")] = "Not differential"
subset_int$K27_status[which(subset_int$K27_status != "No K27 Peak" & subset_int$log2FC_ChIP < -log2(1.5) & subset_int$qval_K27 < 0.1)] = "Depleted FC < -1.5"
subset_int$K27_status[which(subset_int$K27_status != "No K27 Peak" & subset_int$log2FC_ChIP < -log2(2) & subset_int$qval_K27 < 0.1)] = "Depleted FC < -2"
subset_int$K27_status[which(subset_int$K27_status != "No K27 Peak" & subset_int$log2FC_ChIP < -log2(3) & subset_int$qval_K27 < 0.1)] = "Depleted FC < -3"
subset_int$K27_status[which(subset_int$K27_status != "No K27 Peak" & subset_int$log2FC_ChIP > log2(1.5) & subset_int$qval_K27 < 0.1)] = "Overexpressed - FC > 1.5"
subset_int$K27_status[which(subset_int$K27_status != "No K27 Peak" & subset_int$log2FC_ChIP > log2(2) & subset_int$qval_K27 < 0.1)] = "Overexpressed - FC > 2"
subset_int$K27_status[which(subset_int$K27_status != "No K27 Peak" & subset_int$log2FC_ChIP > log2(3) & subset_int$qval_K27 < 0.1)] = "Overexpressed - FC > 3"

tab = subset_int %>% dplyr::select(
  Symbol, log2FC.C2_pers, qval.C2_pers, log2FC_ChIP, qval_K27, K27_status)
WriteXLS(tab, file.path(resSUBdir_K27,"persister_K27_status.xls"), SheetNames = "Persister_K27")

vec = (table(subset_int$K27_status))[c(3,4,1,2)]

pdf(file.path(resSUBdir_K27,paste0("NumPeaksK27_scRNAseq_bulk_sc.pdf")),height=5,width=5)
pie(as.numeric(vec), labels = paste0(names(vec)," ; ",as.numeric(vec)),
    col = c("#999999ff","#D1CE68FD","#67CF85", "#33BD5A", "#00A64B"),
    cex=0.5)
dev.off()



