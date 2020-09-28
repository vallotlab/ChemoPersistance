
##############################################################################################
###############			LIBRARIES AND FUNCTIONS					##############################
##############################################################################################
options(stringsAsFactors=FALSE, width=180)
organism <- "hg38"
GencodeVersion <- ""
chRangeHM=T
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

library(geco.supervised);
library(geco.utils)
library(geco.visu)
library(edgeR); library(ggplot2); library(rgl); library(RColorBrewer); library(genefilter); library(xtable); library(geco.RNAseq); library(WriteXLS); library(data.table); library(stringr); library(limma); library(edgeR);library(dplyr) ####################################################################
library(scTools)
library(Matrix)
library(dplyr)
library(corrplot)
library(geco.utils)
library(geco.visu)
library(geco.unsupervised)
library(scatterplot3d)
library(scater)
library(Rtsne)
library(ccRemover)
library(colorRamps)
library(geco.supervised)
library(viridis)
library(colorRamps)
library(geco.supervised);  library(edgeR); library(ggplot2); library(rgl); library(RColorBrewer); library(genefilter); library(xtable); library(geco.RNAseq); library(WriteXLS); library(data.table); library(stringr); library(limma); library(edgeR)

library(scTools)
library(monocle3)
library(dplyr)
library(WriteXLS)
library(clValid)
library(ape)
library(ConsensusClusterPlus)
devtools::load_all("/media/pprompsy/LaCie/InstitutCurie/Documents/GitLab/R_packages/GeCo.Rpackages/geco.visu/")

MSigDBFile1 <- "/media/pprompsy/LaCie/InstitutCurie/Documents/GitLab/ChromSCape_devel/data/hg38.MSIG.gs.rda" 
MSigDBFile2 <- "/media/pprompsy/LaCie/InstitutCurie/Documents/GitLab/ChromSCape_devel/data/hg38.MSIG.ls.rda"
load(MSigDBFile1)
load(MSigDBFile2)
MSIG.ls = hg38.MSIG.ls
MSIG.gs = hg38.MSIG.gs
### PARAMETERS	
####################################################################


useClusterInfoFromUnsupp <- TRUE ## TRUE , FALSE
####################################################################
### DIRECTORIES and FILES 
####################################################################

######VARIABLES
options(width=180)
maindir <- "/media/pprompsy/LaCie/InstitutCurie/Z_server/Manuscripts/2020_ChemoTolerance/Raw_analysis/scRNAseq/"
resdir <- "/media/pprompsy/LaCie/InstitutCurie/Z_server/Manuscripts/2020_ChemoTolerance/Raw_analysis/scRNAseq/RNA_ChIP";if(!file.exists(resdir)){dir.create(resdir)}
RDatadir_PDX <- "/media/pprompsy/LaCie/InstitutCurie/Z_server/Manuscripts/2020_ChemoTolerance/Raw_analysis/scRNAseq/Results_PDX/Supervised/RData/"
RDatadir_MM468 <- "/media/pprompsy/LaCie/InstitutCurie/Z_server/Manuscripts/2020_ChemoTolerance/Raw_analysis/scRNAseq/Results_MM468/Supervised_persister/RData/"
RDatadir_ChIP_K27 <- "/media/pprompsy/LaCie/InstitutCurie/Z_server/Manuscripts/2020_ChemoTolerance/Figures/ChIP_seq_Bulk_Persisters/K27_peaks_K27/results/Supervised//RData/"
RDatadir_ChIP_K4 <- "/media/pprompsy/LaCie/InstitutCurie/Z_server/Manuscripts/2020_ChemoTolerance/Figures/ChIP_seq_Bulk_Persisters/K4_transcripts_10k/results/Supervised/RData/"
RDatadir_unsup <- "/media/pprompsy/LaCie/InstitutCurie/Z_server/Manuscripts/2020_ChemoTolerance/Raw_analysis/scRNAseq/Results_PDX/Unsupervised//RData/"
RDataSupdir <-  file.path(resdir,"RData");if(!file.exists(RDataSupdir)){dir.create(RDataSupdir)}

## Import annotation file
load("common_over_genes_pers_vs_unt.RData")

K27 = new.env()
K4 = new.env()
PDX = new.env()
MM468 = new.env()
load(file.path(RDatadir_PDX,"Supervised_res_object_edgeR.Rdata"),PDX)
load(file.path(RDatadir_MM468,"Supervised_res_object_edgeR.Rdata"),MM468)
load(file.path(RDatadir_ChIP_K27,"Supervised_res_Limma_object.Rdata"),K27)
load(file.path(RDatadir_ChIP_K4,"Supervised_res_Limma_object.Rdata"),K4)

PDX$overexpressed = PDX$my.res$gene_short_name[which(PDX$my.res$log2FC.persister_6_vs_UNT > 2 & PDX$my.res$qval.persister_6_vs_UNT < 0.01)]
MM468$overexpressed = MM468$my.res$gene_short_name[which(MM468$my.res$log2FC.C2_C4_pers > 2 & MM468$my.res$qval.C2_pers < 0.01)]
K27$depleted = unique(as.character(unlist(str_split_fixed(K27$res$gene_affectation[which(K27$res$log2FC.X5FU2_3_5 < -2 & K27$res$log2FC.X5FU2_3_5 < 0.1)],",",100))))
K27$depleted = K27$depleted[-which(K27$depleted == "")]
K4$enriched = unique(as.character(unlist(str_split_fixed(K4$res$Gene[which(K4$res$log2FC.X5FU6 > 2 & K4$res$qval.X5FU6 < 0.1)],",",100))))
K4$enriched = K4$enriched[-which(K4$enriched == "")]

intersect(K27$depleted,MM468$overexpressed)
intersect(K27$depleted,PDX$overexpressed)
intersect(K4$enriched,MM468$overexpressed)
intersect(K4$enriched,PDX$overexpressed)
