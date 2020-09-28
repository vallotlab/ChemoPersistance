library(scran)
library(devtools)
library(dplyr)
library(DropletUtils)
library(irlba)
library(scater)
library(rscrublet)

# library(corrplot)
# library(geco.utils)
# library(geco.visu)
# library(ConsensusClusterPlus)
# library(geco.unsupervised)
# library(seriation)
# library(heatmap.plus)
# library(scatterplot3d)
# library(colorRamps)
# library(monocle3)
# library(viridis)
# library(colorRamps)
# library(RColorBrewer)

options(stringsAsFactors = F)

# Directories -------------------------------------------------------------

QCdir <- "~/Desktop/scRNAseq_data_local/QC/"
resdir <- "~/Desktop/scRNAseq_data_local/Results_MM468/"
resSUBdir <- file.path(resdir, paste0("Unsupervised_persister")) ; if(!file.exists(resSUBdir)){dir.create(resSUBdir)}
RDatadir <- file.path(resSUBdir,"RData") ; if(!file.exists(RDatadir)){dir.create(RDatadir)}
RDatadirSamples <- file.path(RDatadir,"RData_perSample") ; if(!file.exists(RDatadirSamples)){dir.create(RDatadirSamples)}


# Parameters --------------------------------------------------------------


## packages Seurat or scater attach an enormous amount of bystander packages, they cannot be loaded in the same session without crashing the DLLs.
##Hierarchical clustering
distHC <- c("distPearson","distCosine","euclidean","maximum","manhattan","canberra","binary","minkowski")[1]
methHC <- c("ward","ward.D","ward.D2","single","complete","average")[2]

##ConsClust
repsCC <- 1000
pItemCC <- 0.8
pFeatureCC <- 1
clusterAlgCC <- c("hc","pam","km","kmdist")[1]
distCC <- c("pearson","distCosine","euclidean","manhattan")[1]
innerLinkageCC <- c("ward","ward.D","ward.D2","single","complete","average")[2]
finalLinkageCC <- c("ward","ward.D","ward.D2","single","complete","average")[2]

## Heatmap
chRangeHM <-TRUE # Should be set to TRUE for expression data, FALSE for methylation data
hmColors <- colorRampPalette(c("royalblue","white","indianred1"))(256)
corColors <- colorRampPalette(c("royalblue1","white","indianred1"))(256)





# Load Data ---------------------------------------------------------------


BCmodel <- "MM468"


MM468_initial_ini <- read10xCounts("~/Desktop/scRNAseq_data_local/MM468_initial/filtered_feature_bc_matrix/")


MM468_DMSO5_day67_ini <- read10xCounts("~/Desktop/scRNAseq_data_local/MM468_DMSO5_day67/filtered_feature_bc_matrix/")
MM468_DMSO3_day50_ini <- read10xCounts("~/Desktop/scRNAseq_data_local/MM468_DMSO3_day50/filtered_feature_bc_matrix/")
MM468_5FU3_day50_ini <- read10xCounts("~/Desktop/scRNAseq_data_local/MM468_5FU3_day50/filtered_feature_bc_matrix/")
MM468_5FU3_day77_ini <- read10xCounts("~/Desktop/scRNAseq_data_local/MM468_5FU3_day77/filtered_feature_bc_matrix/")
MM468_5FU3_day202_ini <- read10xCounts("~/Desktop/scRNAseq_data_local/MM468_5FU3_day_202/filtered_feature_bc_matrix/")
MM468_5FU5_day67_ini <- read10xCounts("~/Desktop/scRNAseq_data_local/MM468_5FU5_day67/filtered_feature_bc_matrix/")
MM468_5FU5_day171_ini <- read10xCounts("~/Desktop/scRNAseq_data_local/MM468_5FU5_day171/filtered_feature_bc_matrix/")

MM468_5FU6_day33_ini <- read10xCounts("~/Desktop/scRNAseq_data_local/MM468_5FU6_day33/filtered_feature_bc_matrix/")
MM468_5FU6_day214_ini <- read10xCounts("~/Desktop/scRNAseq_data_local/MM468_5FU6_day214/filtered_feature_bc_matrix/")

MM468_5FU6_UNC_day33_ini <- read10xCounts("~/Desktop/scRNAseq_data_local/MM468_5FU6_UNC_day33/filtered_feature_bc_matrix/")
MM468_UNC_day33_ini <- read10xCounts("~/Desktop/scRNAseq_data_local/MM468_UNC_day33/filtered_feature_bc_matrix/")

#MM468_5FU6_dayXX_ini <- read10xCounts("~/Google Drive/scRNAseq_data/MM468_5FU6_dayXX/filtered_feature_bc_matrix/")

# QC, normalization  ------------------------------------------------------


for(i in 8:10){
  nameS <- c("MM468_DMSO3_day50","MM468_5FU3_day50","MM468_5FU3_day77","MM468_5FU3_day202","MM468_DMSO5_day67","MM468_5FU5_day67","MM468_5FU5_day171","MM468_5FU6_day33","MM468_5FU6_UNC_day33","MM468_UNC_day33","MM468_5FU6_day214","MM468_initial")[i]
  umi <- list(MM468_DMSO3_day50_ini,MM468_5FU3_day50_ini,MM468_5FU3_day77_ini,MM468_5FU3_day202_ini,MM468_DMSO5_day67_ini,MM468_5FU5_day67_ini,MM468_5FU5_day171_ini,MM468_5FU6_day33_ini,MM468_5FU6_UNC_day33_ini,MM468_UNC_day33_ini,MM468_5FU6_day214_ini,MM468_initial_ini)[[i]]
  #umi <- list(MM468_5FU6_day214_ini,MM468_initial_ini)[[i-10]]
  
  keep_feat <- rowSums(counts(umi)>0)>0 #remove genes that are not expressed in any cells
  umi <- umi[keep_feat,]

  mt <- grepl("^MT-", rowData(umi)[[2]])
  # mt <- grepl("MT-", rowData(umi)[[2]])
  # mt[which(rowData(umi)[[2]]=="INMT-MINDY4")] <- FALSE
  ercc <- grepl("ERCC", rowData(umi)[[2]])
  isSpike(umi,"ERCC") <- ercc
  rb.genes <- rowData(umi)[[1]][grep("^RP[SL]",rowData(umi)[[2]])]

  umi <- scater::calculateQCMetrics(umi,feature_controls=list(ERCC=ercc,MT=mt))

  ##Manual quality filters:
  #Distribution of number of total detected features should be normal, get rid of the heavy right tail (<6000)
  umi$filter_by_expr_features <- (umi$total_features_by_counts<8000 & umi$total_features_by_counts>3000)
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


  pdf(file.path(QCdir,paste0(nameS,"_QC_scRNAseq.pdf")))
  scater::plotColData(umi,x="total_features_by_counts",y="pct_counts_MT")
  scater::plotColData(umi,x="total_features_by_counts",y="pct_counts_ERCC")
  hist(umi$total_features_by_counts,breaks=100,xlim=c(0,10000),main=" Number of detected genes",xlab="NUMBER of GENES")
  abline(v=6000,col="red")
  abline(v=1500,col="red")
  hist(umi$total_counts,breaks=100, xlim=c(0,300000),main="Library size",xlab="TOTAL COUNTS")
  M <- median(umi$total_counts)
  abline(v=100000,col="red")

  dev.off()


  #assign(paste0("filt.",nameS),umi.qc)

  unc <- counts(umi.qc)
  rownames(unc) <- rowData(umi.qc)$Symbol
  assign(paste0("counts.",nameS),unc)
  set.seed(100)
  test <- doubletCells(umi.qc)
  metadata <- data.frame(barcode= substr(colData(umi.qc)$Barcode,0,16),cell_id=sprintf(paste0(nameS,"_c%d"), 1:dim(umi.qc)[2]), sample_id=rep(nameS,dim(umi.qc)[2]), total_features=colData(umi.qc)$total_features_by_counts,
                         total_counts=colData(umi.qc)$total_counts,
                         doublet_score=test)

  ###############################################
  #Add lineage info###############################
  ###############################################

 
  reference_lenti <- read.csv("~/Documents/Vallot_lab/R_custom_scripts/LineageTracing/LG22_filtered.fa",comment.char = ">",header=F)
  fl <- c("CTAGAACACTCGAGATCAG","TGATACATCATACCACAT","CTGATCTCGAGTGTTCTAG","ATGTGGTATGATGTATCA") # Top and Bot, and RTop, RBot flanking sequence
  I = 1
  D = 1
  S = 1

  drop_barcode <- read.table(file = paste0("~/Desktop/scRNAseq_data_local/LineageBarcodes/matched_sequences_",nameS),sep="\t")
  #drop_barcode$barcode <- substr(drop_barcode$V3,0,16)
  
  drop_barcode$barcode <- substr(drop_barcode$V2,0,16) #for files generated by Pacôme
  
  length(unique(drop_barcode$barcode))
  lineage_barcode <- read.table(file = paste0("~/Desktop/scRNAseq_data_local/LineageBarcodes/sequence_withTop_",nameS))

  pos_top = attr(adist(fl[1], lineage_barcode[,1], costs = c(I,D,S), fixed=F, partial=T, counts=T), "offsets")
  lineage_barcode$pos_FTop <- pos_top[,,"first"]
  pos_bot = attr(adist(fl[4], lineage_barcode[1,1], costs = c(I,D,S), fixed=F, partial=T, counts=T), "offsets")
  lineage_barcode$pos_RBot <- pos_bot[,,"first"]
  lineage_barcode$count <- (lineage_barcode$pos_RBot - lineage_barcode$pos_FTop)
  lineage_barcode$LBC <- substr(lineage_barcode$V1,lineage_barcode$pos_FTop+19,lineage_barcode$pos_FTop+38)
  lineage_barcode$length_barcode <-  nchar(lineage_barcode$LBC)
  
  lineage_barcode$LBC[lineage_barcode$length_barcode < 14] <- "short"
  

  unique_lineage <- unique(lineage_barcode$LBC[which(lineage_barcode$LBC!="short")]);

  list_agrep  = sapply(unique_lineage,function(x) agrep(x, reference_lenti$V1, max.distance=2, costs = c(I,D,S), value=F))
  indice_agrep <- unlist(list_agrep)

  if(length(grep("\\d",names(indice_agrep)))>0) {
    indice_agrep = indice_agrep[-grep("\\d",names(indice_agrep))] # remove sequence matching multiple ref
  }

  lineage_barcode$lineage_indice <- indice_agrep[match(substring(lineage_barcode$LBC,1,20),names(indice_agrep))]
  lineage_barcode$tag <- reference_lenti$V1[lineage_barcode$lineage_indice]
  lineage_barcode$lineage_indice[is.na(lineage_barcode$lineage_indice) & lineage_barcode$length_barcode==20] <- "unknown"
  lineage_barcode$lineage_indice[is.na(lineage_barcode$lineage_indice)] <- "N"
  lineage_barcode$tag[lineage_barcode$lineage_indice=="unknown"] <- lineage_barcode$LBC[lineage_barcode$lineage_indice=="unknown"]
  
  
  dim(lineage_barcode)

  metadata$lineage_barcode <- lineage_barcode$tag[match(metadata$barcode,drop_barcode$barcode)]
  metadata$lineage_barcode_ref <- lineage_barcode$lineage_indice[match(metadata$barcode,drop_barcode$barcode)]
  metadata$lineage_barcode_ref[metadata$lineage_barcode_ref=="N"] <- NA
  
  #metadata$lineage_barcode <- NA
  ################################################
  ################################################
  
  assign(paste0("metadata.",nameS),metadata)
  save(list=c(paste0("counts.",nameS),paste0("metadata.",nameS)),file=file.path(RDatadirSamples,paste0(nameS,".RData")))
  
  
}

################################################
############Add AMarie' lineage info############
################################################

for(i in 5:11){
  nameS <- c("MM468_DMSO3_day50","MM468_5FU3_day50","MM468_5FU3_day77","MM468_5FU3_day202","MM468_5FU5_day67","MM468_5FU5_day171","MM468_5FU6_day33","MM468_5FU6_UNC_day33","MM468_UNC_day33","MM468_5FU6_day214","MM468_initial")[i]
  load(file.path(RDatadirSamples,paste0(nameS,".RData")))
  lineage_AM <- read.csv(paste0("~/Desktop/scRNAseq_data_local/LineageBarcodes/","consensus_cell_10xbc_and_vbc_",nameS,"_20_4_annot.csv"))
  lineage_AM$BC_10x <- substr(lineage_AM$BC_10x,0,16)
  
  test <- get(paste0("metadata.",nameS))
  colnames(test)[1] <- "BC_10x"
  test <- left_join(test,lineage_AM)
  
  assign(paste0("metadata.",nameS),test)
  save(list=c(paste0("counts.",nameS),paste0("metadata.",nameS)),file=file.path(RDatadirSamples,paste0(nameS,".RData")))
  rm(list=(paste0("counts.",nameS)))
}



