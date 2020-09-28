library(scran)
library(devtools)
library(dplyr)
library(DropletUtils)
library(irlba)
# library(corrplot)
# library(geco.utils)
# library(geco.visu)
# library(ConsensusClusterPlus)
# library(geco.unsupervised)
# library(seriation)
# library(heatmap.plus)
# library(scatterplot3d)
library(scater)
library(rscrublet)
# library(colorRamps)
# library(monocle3)
# library(viridis)
# library(colorRamps)
# library(RColorBrewer)

options(stringsAsFactors = F)

# Directories -------------------------------------------------------------

QCdir <- "/media/pacome/LaCie/InstitutCurie/Z_server/Manuscripts/2020_ChemoTolerance/Raw_analysis/scRNAseq/QC/"; if(!dir.exists(QCdir)){dir.create(QCdir)}
resdir <- "/media/pacome/LaCie/InstitutCurie/Z_server/Manuscripts/2020_ChemoTolerance/Raw_analysis/scRNAseq/Results_PDX/"
resSUBdir <- file.path(resdir, paste0("Unsupervised")) ; if(!file.exists(resSUBdir)){dir.create(resSUBdir)}
RDatadir <- file.path(resSUBdir,"RData") ; if(!file.exists(RDatadir)){dir.create(RDatadir)}
RDatadirSamples <- file.path(RDatadir,"RData_perSample") ; if(!file.exists(RDatadirSamples)){dir.create(RDatadirSamples)}


# Parameters --------------------------------------------------------------



# Load Data ---------------------------------------------------------------


BCmodel <- "HBCx95"
#aal have been analyzed with hg19 and mm10 alignments

HBCx95_m00_UNT_ini <- read10xCounts("/media/pacome/LaCie/InstitutCurie/Z_server/Manuscripts/2020_ChemoTolerance/Raw_analysis/scRNAseq/Results_PDX/Data/HBCx95_m00_UNT/hg19_mm10/outs/filtered_gene_bc_matrices_HBCx-95/hg19/")
HBCx95_m30_CAPAR_ini <- read10xCounts("/media/pacome/LaCie/InstitutCurie/Z_server/Manuscripts/2020_ChemoTolerance/Raw_analysis/scRNAseq/Results_PDX/Data/HBCx95_m30_CAPAR/filtered_feature_bc_matrix/")
HBCx95_m40_CAPAR_ini <- read10xCounts("/media/pacome/LaCie/InstitutCurie/Z_server/Manuscripts/2020_ChemoTolerance/Raw_analysis/scRNAseq/Results_PDX/Data/HBCx95_m40_CAPAR/hg19_mm10/outs/filtered_gene_bc_matrices_HBCx-95_CAPAR/hg19/")
HBCx95_m6_CAPAS_ini <- read10xCounts("/media/pacome/LaCie/InstitutCurie/Z_server/Manuscripts/2020_ChemoTolerance/Raw_analysis/scRNAseq/Results_PDX/Data/HBCx95_m6_CAPAS/filtered_gene_bc_matrices/hg19/")
HBCx95_m9_50_REC_ini <- read10xCounts("/media/pacome/LaCie/InstitutCurie/Z_server/Manuscripts/2020_ChemoTolerance/Raw_analysis/scRNAseq/Results_PDX/Data/HBCx95_m9_50_REC/filtered_feature_bc_matrix/")
HBCx95_persister_6souris_ini <- read10xCounts("/media/pacome/LaCie/InstitutCurie/Z_server/Manuscripts/2020_ChemoTolerance/Raw_analysis/scRNAseq/Results_PDX/Data/HBCx95_persister_6souris_hg19/filtered_feature_bc_matrix_merged/")
HBCx95_persister_2souris_ini <- read10xCounts("/media/pacome/LaCie/InstitutCurie/Z_server/Manuscripts/2020_ChemoTolerance/Raw_analysis/scRNAseq/Results_PDX/Data/HBCx95_persister_2souris/filtered_feature_bc_matrix/")
HBCx95_v3_UNT_ini <- read10xCounts("/media/pacome/LaCie/InstitutCurie/Z_server/Manuscripts/2020_ChemoTolerance/Raw_analysis/scRNAseq/Results_PDX/Data/HBCx95_p26/filtered_feature_bc_matrix/")
gc()

# QC, normalization  ------------------------------------------------------
for(i in 1:8){
  nameS <- c("HBCx95_m00_UNT","HBCx95_m30_CAPAR","HBCx95_m40_CAPAR","HBCx95_m6_CAPAS","HBCx95_m9_50_REC","HBCx95_persister_6souris","HBCx95_persister_2souris","HBCx95_v3_UNT")[i]
  kit_version <- c("v2","v3","v2","v3","v3","v3","v3","v3")[i]
  umi <- list(HBCx95_m00_UNT_ini,HBCx95_m30_CAPAR_ini,HBCx95_m40_CAPAR_ini,HBCx95_m6_CAPAS_ini,HBCx95_m9_50_REC_ini,HBCx95_persister_6souris_ini,HBCx95_persister_2souris_ini,HBCx95_v3_UNT_ini)[[i]]
  keep_feat <- rowSums(counts(umi)>0)>0 #remove genes that are not expressed in any cells
  umi <- umi[keep_feat,]
  
  genome <- "hg19"
  hg19 <- grepl(paste0("^",genome),rowData(umi)[[2]])
  umi <- umi[hg19,]
  
  mt <- grepl("^hg19_MT-", rowData(umi)[[2]])
  ercc <- grepl("ERCC", rowData(umi)[[2]])
  isSpike(umi,"ERCC") <- ercc
  rb.genes <- rowData(umi)[[1]][grep("^RP[SL]",rowData(umi)[[2]])]
  
  umi <- scater::calculateQCMetrics(umi,feature_controls=list(ERCC=ercc,MT=mt))
  
  ##Manual quality filters:
  #Distribution of number of total detected features should be normal, get rid of the heavy right tail (<6000)
  umi$filter_by_expr_features <- (umi$total_features_by_counts<8000 & umi$total_features_by_counts>2500) 
  #Distribution of number of total number of counts
  umi$filter_by_lib_size <- (umi$total_counts<100000) 
  #MT contamination
  umi$filter_by_MT_conta <- (umi$pct_counts_MT<15) 
  #Spike-in 
  umi$filter_by_spikein <- (umi$pct_counts_ERCC<0.05) 
  
  umi$use <- (umi$filter_by_lib_size & umi$filter_by_MT_conta & umi$filter_by_spikein & umi$filter_by_expr_features)
  
  
  umi.qc <- umi[,umi$use]
  # rowData(umi.qc)$use2 <- apply(counts(umi.qc),1,function(x) length(x[x>1])>=2) #add a field to select genes with more than 1 transcript in at least 2 non-outlier ce
  # umi.qc <- umi.qc[rowData(umi.qc)$use2,] 
  
  #remove MT & ERCC genes from dataset for further analysis 
  umi.qc <- umi.qc[!rowData(umi.qc)$is_feature_control,]
  dim(umi.qc)
  
  #Normalization by library size only
  sizeFactors(umi.qc) <- librarySizeFactors(umi.qc)
  umi.qc <- scater::normalize(umi.qc) #normalized log2 expression values
  umi.qc <- scater::normalize(umi.qc,return_log=FALSE) #normalized expression values
  
  #Detect the most variant genes
  # fit <- trendVar(umi.qc, parametric=TRUE, use.spikes=FALSE) 
  # dec <- decomposeVar(umi.qc, fit)
  
  pdf(file.path(QCdir,paste0(nameS,"_QC_scRNAseq.pdf")))
  scater::plotColData(umi,x="total_features_by_counts",y="pct_counts_MT")
  scater::plotColData(umi,x="total_features_by_counts",y="pct_counts_ERCC")
  hist(umi$total_features_by_counts,breaks=100,xlim=c(0,10000),main=" Number of detected genes",xlab="NUMBER of GENES")
  abline(v=8000,col="red")
  abline(v=2500,col="red")
  hist(umi$total_counts,breaks=100, xlim=c(0,300000),main="Library size",xlab="TOTAL COUNTS")
  M <- median(umi$total_counts)
  abline(v=100000,col="red")
  
  dev.off()
  
  
  # dec.umi.qc <- dec
  # dec.umi.qc$Symbol <- rowData(umi.qc)$Symbol
  # dec.umi.qc <- dec.umi.qc[order(dec.umi.qc$bio, decreasing=TRUE),]
  # head(dec.umi.qc)
  
  assign(paste0("filt.",nameS),umi.qc)
  #assign(paste0("dec.",nameS),dec.umi.qc)
  
  unc <- counts(umi.qc)
  rownames(unc) <- rowData(umi.qc)$Symbol
  assign(paste0("counts.",nameS),unc)
  
  # doubletsCell <- scrubDoublets(as.matrix(unc))
  metadata <- data.frame(barcode= substr(colData(umi.qc)$Barcode,0,16),cell_id=sprintf(paste0(nameS,"_c%d"), 1:dim(umi.qc)[2]), sample_id=rep(nameS,dim(umi.qc)[2]), total_features=colData(umi.qc)$total_features_by_counts,
                         total_counts=colData(umi.qc)$total_counts,kit=rep(kit_version,dim(umi.qc)[2]) )
  
  
  assign(paste0("metadata.",nameS),metadata)
  save(list=c(paste0("counts.",nameS),paste0("metadata.",nameS)),file=file.path(RDatadirSamples,paste0(nameS,".RData")))
  
}





