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
                                      "DA_grouped_Persister_vs_DSMO_with_bulk_10k.xlsx"), sheet = 3)

sc_bulk_ChIP_DA = sc_bulk_ChIP_DA %>% group_by(Gene) %>% slice_max(abs(log2FC.Persister))

log2FC_thresholds <- log2(c(2,3,4))

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

vec = (table(subset_int$K27_status))[c(4,5,1,2,3)]

pdf(file.path(resSUBdir_K27,paste0("NumPeaksK27_scRNAseq_bulk_sc.pdf")),height=5,width=5)
pie(as.numeric(vec), labels = paste0(names(vec)," ; ",as.numeric(vec)),
    col = c("#999999ff","#D1CE68FD","#67CF85", "#33BD5A", "#00A64B"),
    cex=0.5)
dev.off()

# Pie Chart overexpressed in K27 
subset_int = annot10k_K27_byGene[annot10k_K27_byGene$log2FC_K27 < -log2(3) & annot10k_K27_byGene$qvalue_K27 < 0.1,]
dim(my.res)
subset_int = subset_int[which(subset_int$gene %in% my.res$Symbol),]
subset_int$scRNA_status = "Not Overexpressed"
subset_int$log2FC_scRNA = res.scRNA$log2FC.C2_pers[match(subset_int$gene, res.scRNA$Symbol)]
subset_int$qvalue_scRNA = res.scRNA$qval.C2_pers[match(subset_int$gene, res.scRNA$Symbol)]

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

# K27, K4, scRNA 
differential_results_K4_all <- readxl::read_xlsx(file.path(
  maindir,"output","bulk_ChIPseq","MM468","K4_transcripts_10k","Supervised","Tables",
  "Differential_analysis_Limma.xlsx"), sheet = 3)
differential_results_K4_all$ID = sub("-",":",gsub("_","-",differential_results_K4_all$ID))
differential_results_K4_byGene = differential_results_K4_all %>% group_by(Gene) %>%  slice_max(order_by = abs(log2FC.X5FU6), n = 1)

res.scRNA_over = res.scRNA %>% filter(log2FC.C2_pers > log2(3) & qval.C2_pers < 0.01)

differential_results_K4_byGene = differential_results_K4_byGene[match(res.scRNA_over$Symbol,differential_results_K4_byGene$Gene),]
differential_results_K27_byGene = annot10k_K27_byGene[match(res.scRNA_over$Symbol,annot10k_K27_byGene$gene),]

res.scRNA_over$log2FC_K27 = -differential_results_K27_byGene$log2FC_K27
res.scRNA_over$log2FC_K4 = differential_results_K4_byGene$log2FC.X5FU6

res.scRNA_over %>% tidyr::gather(log2FC,Mark,log2FC_K27,log2FC_K4)
ggplot(res.scRNA_over %>% tidyr::gather(Mark,log2FC,log2FC_K27,log2FC_K4)) +
  geom_col(aes(x=Symbol, y = log2FC, fill = Mark))

# K27 peak based intersection
res.scRNA <- my.res
annot_K27_gene = read.table(file.path(maindir,"annotation","MM468_peaks_K27_gene.bed"), sep="\t")
load(file.path(maindir,"output","bulk_ChIPseq","MM468","K27_peaks_K27","Supervised","RData","Supervised_res_Limma_object.Rdata"))
res.ChIP <- res[,c(1,2,3,4,19,21)]
res.ChIP$gene = annot_K27_gene$V4

res.ChIP_splitGene = res.ChIP %>% tidyr::separate_rows(gene)
res.ChIP_splitGene = res.ChIP_splitGene[!is.na(res.ChIP_splitGene$gene),]

#Within expressed genes in the model system, n=14924 genes, 2083 have a H3K27me3 peak over their tss or gene_body
res.scRNA$log2FC_ChIP <- res.ChIP_splitGene$log2FC.X5FU2_3_5[match(res.scRNA$Symbol,res.ChIP_splitGene$gene)]
res.scRNA$qvalue_ChIP <- res.ChIP_splitGene$qval.X5FU2_3_5[match(res.scRNA$Symbol,res.ChIP_splitGene$gene)]

# Pie Chart K4 in overexpressed
subset_int = res.scRNA[res.scRNA$log2FC.C2_pers>log2(3) & res.scRNA$qval.C2_pers<0.01,]

# Read in K4 peaks
annot_K4 <- read.table(
  gzfile(file.path(
    maindir,"annotation","MM468_peaks_K4_29714_peaks.bed.gz"),
    "gencode.v34.transcripts10k_K27.bed"), sep="\t")
colnames(annot_K4) <- c("chr","start","end","length")
annot_K4 = as(annot_K4,"GRanges")
gencodeTSS = as(annot10k_K27[,c(1:3,5)],"GRanges")
hits = findOverlaps(gencodeTSS, annot_K4)
gencodeTSS$K4_status = FALSE
gencodeTSS$K4_status[queryHits(hits)] = TRUE

table(gencodeTSS$K4_status)
gencodeTSS_byGene = as.data.frame(gencodeTSS) %>% group_by(gene) %>%
  summarise(K4_status = any(K4_status)) # Select only the transcript associated with top differential K27 log2FC 
table(gencodeTSS_byGene$K4_status)
subset_int$K4_status = gencodeTSS_byGene$K4_status[match(subset_int$Symbol,gencodeTSS_byGene$gene)]
subset_int$K4_status[!subset_int$K4_status] = "No K4 peak"
subset_int$K4_status[which(subset_int$K4_status=="TRUE")] = "Not differential"
subset_int$log2FC_K4 = differential_results_K4_byGene$log2FC.X5FU6[match(subset_int$Symbol,differential_results_K4_byGene$Gene)]
subset_int$log2FC_K4[which(subset_int$K4_status == "No K4 peak")] = NA
subset_int$K4_status[which(subset_int$log2FC_K4 > log2(1.5))] = "Enriched FC > 1.5"
subset_int$K4_status[which(subset_int$log2FC_K4 > log2(2))] = "Enriched FC > 2"
subset_int$K4_status[which(subset_int$log2FC_K4 > log2(3))] = "Enriched FC > 3"
subset_int$K4_status[which(subset_int$log2FC_K4 < -log2(1.5))] = "Depleted FC > 1.5"
subset_int$K4_status[which(subset_int$log2FC_K4 < -log2(2))] = "Depleted FC > 2"
subset_int$K4_status[which(subset_int$log2FC_K4 < -log2(3))] = "Depleted FC > 3"
table(subset_int$K4_status)
vec = (table(subset_int$K4_status))[c(7,8,1:6)]

pdf(file.path(resSUBdir_K4,paste0("NumPeaksK4_scRNAseq_10k.pdf")),height=5,width=5)
pie(as.numeric(vec), labels = paste0(names(vec)," ; ",as.numeric(vec)),
    col = c("#999999ff","#D1CE68FD","#67CF85", "#33BD5A", "#00A64B","#FC8A8A", "#FA3F3F", "#FF0000"),
    cex=0.5)
dev.off()

for(log2FC_threshold in log2FC_thresholds){
  
  sig.scRNA <- res.scRNA[res.scRNA$log2FC.C2_pers>log2FC_threshold & res.scRNA$qval.C2_pers<0.01,]
  sig.peaks <- res.scRNA$Symbol[which(abs(res.scRNA$log2FC_ChIP)>2 & res.scRNA$qvalue_ChIP<0.1)]

  #Study the correlation between expression and changes in H3K27me3 enrichment on peaks that are diff enriched
  subset_int <- res.scRNA[res.scRNA$Symbol %in% sig.peaks,]
  
  n_tot = nrow(sig.scRNA)
  n_peaks = length(which(!is.na(sig.scRNA$log2FC_ChIP)))
  n_differential_peaks = length(which(sig.scRNA$log2FC_ChIP < -log2FC_threshold & sig.scRNA$qvalue_ChIP <0.1))
  
  pdf(file.path(resSUBdir_K27,paste0("NumPeaksK27_",round(log2FC_threshold,2),"_scRNAseq.pdf")),height=5,width=5)
  pie(c(n_tot - n_peaks, n_peaks - n_differential_peaks, n_differential_peaks),
      labels = c(paste0(n_tot - n_peaks," - No peaks"),
                 paste0(n_peaks - n_differential_peaks, " - K27 Peaks"),
                 paste0(n_differential_peaks, " - K27 Depleted Peaks")),
      col = c("#999999ff", "#3b9ab2ff", "#78b7c5ff"))
  dev.off()
  
  top_demethylated <- head(subset_int$Symbol[order(subset_int$log2FC.C2_pers,decreasing=T)],n=15)
  
  pdf(file.path(resSUBdir_K27,paste0("Log2FC_",round(log2FC_threshold,2),"_scRNAseq_for_signigicantPeaks.pdf")),height=5,width=5)
  sp <- ggscatter(subset_int, y = "log2FC.C2_pers", x = "log2FC_ChIP",
                  add = "reg.line",  # Add regressin line
                  add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
                  conf.int = TRUE,
                  label="Symbol",repel=TRUE,
                  label.select= top_demethylated) #add size of dot function of initial expression in DMSO
  # Add correlation coefficient
  print(sp + stat_cor(method = "pearson",label.x=1,label.sep = "\n"))
  dev.off()
  
  res.scRNA$decile <- ntile(res.scRNA$log2FC.C2_pers, 10)  
  figure <- res.scRNA %>% group_by(decile) %>% 
    summarise(ChIP=length(which(log2FC_ChIP<(-log2FC_threshold) & qvalue_ChIP<0.1)),
              ChIP_res=length(which(log2FC_ChIP<(-log2FC_threshold) & qvalue_ChIP<0.1)),
              expression=max(log2FC.C2_pers))
  
  pdf(file.path(resSUBdir_K27,paste0("ExpressionDecileforSignificantDepletedPeaks_logFC",round(log2FC_threshold,2),".pdf")),height=5,width=5)
  
  print(
    ggplot(data=figure, aes(x=decile , y=ChIP, fill=decile)) +
      geom_rect(data=figure,aes(xmin=decile-0.5,xmax=decile+0.5,ymin=-Inf,ymax=Inf),fill=alpha(colorRampPalette(c("royalblue","white","indianred1"))(10),0.8))+ 
      geom_point(data=figure, aes(x=decile, y=ChIP)) + theme_classic() + theme(legend.position = "none") + xlab("log2FC persister vs DMSO decile") + ylab("number of significantly depleted peaks")  + geom_text(aes(x=decile,y=40,label=round(expression,2)))
  )
  dev.off()
}

##### Comparison with H3K4me3 #####

load(file.path(maindir,"output","bulk_ChIPseq","MM468","K4_transcripts_10k","Supervised","RData","Supervised_res_Limma_object.Rdata"))


res.scRNA <- my.res
res.ChIP <- res

trial <- function(x) {
  indice <- grep(pattern = x, res.ChIP$transcripts)
  res.ChIP$ID[indice]
}
test <- sapply(res.scRNA$Symbol,FUN=trial)

res.scRNA$ChIP <- 0
for(i in 1:length(res.scRNA$Symbol)){
  if(length(test[[i]])!=0) res.scRNA$ChIP[i] <- c(test[[i]])
}

#Within expressed genes in the model system, n=14924 genes, 2083 have a H3K27me3 peak over their tss or gene_body
res.scRNA$peak_affectation <- res.ChIP$peak_affectation[match(res.scRNA$ChIP,res.ChIP$ID)]
res.scRNA$log2FC_ChIP <- res.ChIP$log2FC.X5FU6[match(res.scRNA$ChIP,res.ChIP$ID)]
res.scRNA$qvalue_ChIP <- res.ChIP$qval.X5FU6[match(res.scRNA$ChIP,res.ChIP$ID)]
res.scRNA$bivalent <- res.ChIP$bivalent[match(res.scRNA$ChIP,res.ChIP$ID)]

for(log2FC_threshold in log2FC_thresholds){
  #sig.scRNA <- res.scRNA[abs(res.scRNA$log2FC.persisterall)>1 & res.scRNA$qval.persisterall<0.01,]
  sig.peaks <- res.ChIP$ID[abs(res.ChIP$log2FC.X5FU6)>log2FC_threshold & res.ChIP$qval.X5FU6<0.1]
  
  #Study the correlation between expression and changes in H3K27me3 enrichment on peaks that are diff enriched
  subset_int <- res.scRNA[res.scRNA$ChIP %in% sig.peaks,]
  
  top_demethylated <- head(subset_int$Symbol[order(subset_int$log2FC.persisterall,decreasing=T)],n=15)
  
  pdf(file.path(resSUBdir_K4,paste0("Log2FC_",round(log2FC_threshold,2),"_scRNAseq_for_signigicantPeaks.pdf")),height=5,width=5)
  sp <- ggscatter(subset_int, y = "log2FC.persisterall", x = "log2FC_ChIP",
                  add = "reg.line",  # Add regressin line
                  add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
                  conf.int = TRUE,
                  label="Symbol",repel=TRUE,
                  label.select= top_demethylated) #add size of dot function of initial expression in DMSO
  # Add correlation coefficient
  print(sp + stat_cor(method = "pearson",label.x=1,label.sep = "\n"))
  
  sp <- ggscatter(subset_int, y = "log2FC.persisterall", x = "log2FC_ChIP",
                  add = "reg.line",  # Add regressin line
                  add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                  conf.int = TRUE,
                  label="Symbol",repel=TRUE,
                  label.select= top_demethylated) #add size of dot function of initial expression in DMSO
  # Add correlation coefficient
  print(sp + stat_cor(method = "pearson",label.x=1,label.sep = "\n"))
  
  sp <- ggscatter(subset_int, y = "log2FC.persisterall", x = "log2FC_ChIP",
                  add = "reg.line",  # Add regressin line
                  add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
                  conf.int = TRUE) #add size of dot function of initial expression in DMSO
  # Add correlation coefficient
  print(sp + stat_cor(method = "pearson",label.x=1,label.sep = "\n"))
  dev.off()
  
  # sp + geom_density_2d()
  # 
  # sp + stat_density_2d(aes(fill = ..level..), geom = "polygon")
  # 
  
  #Evaluate how much of diff expression is linked to H3K27me3
  
  res.scRNA$decile <- ntile(res.scRNA$log2FC.persisterall, 10)  
  figure <- res.scRNA %>% group_by(decile) %>% 
    summarise(ChIP=length(which(log2FC_ChIP>(log2FC_threshold) & qvalue_ChIP<0.1)),
              ChIP_res=length(which(log2FC_ChIP>(log2FC_threshold) & qvalue_ChIP<0.1)),
              expression=max(log2FC.persisterall))
  
  pdf(file.path(resSUBdir_K4,paste0("ExpressionDecileforSignificantDepletedPeaks_logFC",round(log2FC_threshold,2),".pdf")),height=5,width=5)
  
  print(
    ggplot(data=figure, aes(x=decile , y=ChIP, fill=decile)) +
      geom_rect(data=figure,aes(xmin=decile-0.5,xmax=decile+0.5,ymin=-Inf,ymax=Inf),fill=alpha(colorRampPalette(c("royalblue","white","indianred1"))(10),0.8))+ 
      geom_point(data=figure, aes(x=decile, y=ChIP)) + theme_classic() + theme(legend.position = "none") + xlab("log2FC persister vs DMSO decile") + ylab("number of significantly enriched peaks")  + geom_text(aes(x=decile,y=40,label=round(expression,2)))
  )
  dev.off()
}
