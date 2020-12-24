###### LIBRARIES ######

# Packages
library(ChromSCape)
library(devtools)
library(DropletUtils)
library(irlba)
library(corrplot)
library(R.utils)
library(scater)
library(Rtsne)
library(ccRemover)
library(viridis)
library(colorRamps)
library(RColorBrewer)
library(edgeR)
library(gplots)
library(ggplot2)
library(RColorBrewer)
library(genefilter)
library(xtable)
library(WriteXLS)
library(data.table)
library(stringr)
library(limma)
#library(edgeR)
library(monocle3)
library(dplyr)
library(Seurat)
library(dendextend)
library(clValid)
library(ape)
library(ConsensusClusterPlus)
library(Matrix)
library(genefilter)
library(ggpubr)
library(eulerr)
library(GenomicRanges)
library(Sushi)
library(kableExtra)
library(colorspace)
library(forcats)

#Geco packages
library(geco.supervised)
library(geco.RNAseq)
library(geco.ChIPseq)
library(geco.utils)
library(geco.visu)
library(scTools)
library(geco.unsupervised)

###### GLOBAL VARIABLES ######
pcaText <- TRUE
annotText <- "Sample"
hcText <- "Name"  ## column used in hierarchical clustering # change names from Sample_x to actual sample name
centering <- c("none","mean","median")[1]

##Hierarchical clustering
distHC <- c("distPearson","distCosine","euclidean","maximum","manhattan","canberra","binary","minkowski")[1]
methHC <- c("ward","ward.D","ward.D2","single","complete","average")[2]
GencodeVersion <- ""

##ConsClust
repsCC <- 1000
pItemCC <- 0.8
pFeatureCC <- 1
clusterAlgCC <- c("hc","pam","km","kmdist")[1]
distCC <- c("pearson","distCosine","euclidean","manhattan")[1]
innerLinkageCC <- c("ward","ward.D","ward.D2","single","complete","average")[2]
finalLinkageCC <- c("ward","ward.D","ward.D2","single","complete","average")[2]
useClusterInfoFromUnsupp <- TRUE

## Heatmap
chRangeHM <-FALSE # Should be set to TRUE for expression data, FALSE for methylation data
hmColors <- colorRampPalette(c("royalblue","white","indianred1"))(256)
corColors <- colorRampPalette(c("royalblue","white","indianred1"))(256)
# corColors <- colorRampPalette(c("indianred","white","forestgreen"))(256)

#metadataFile <- file.path(ProjectDir, "Scripts", expType, paste0("metadata_",expType,"_", annoType, ".txt") )
maxKHC <- 10
maxKCC <- 10

###### PATHWAYS ######
hallmark_EMT <- c("FN1","TAGLN","NNMT","ELN","INHBA","MYL9","FBLN2","TGFBR3","CRLF1")

KEGG_TGFB_FC1 <- c("INHBB","INHBA","CREBBP","ROCK1","ID4","RBL2","SMAD4","SMAD1","SMAD2","MAPK1","SMURF1","EP300","TGFB2","CDKN2B","PPP2R1A","SMAD5","FST","TGFBR2","SMAD7","BMPR2","BMPR1B","ACVR2B","ACVR2A","SP1","BMP6","THBS1")
TGF_pathway <- c("TGFBR3","TGFBR2","TGFB2","INHBA","INHBB")

KEGG_TGFB_FC2 <-  c("CDKN2B","ID4","BMP6","INHBA","INHBB")

hallmark_TNFA <- c("CXCL2","ATF3","NFKBIA","TNFAIP3","IER3","CCL20","MAFF","BIRC3","PLAUR","ZFP36","PPP1R15A","NR4A2","MAP3K8","IL6","FOSL1","NR4A1","REL","NR4A3","GEM","GADD45A","PLK2","BHLHE40","NFIL3","TRIB1","TIPARP","INHBA","SGK1","FOSL2","CDKN1A","CYR61","IFIT2","IRF1","FOS","KLF4","LAMB3","CCNL1","CLCF1","IRS2","DUSP4","MXD1")

hallmark_IFN <- c("MX1","BST2","OASL","IFI44L","IFI27","PTGS2","HLA-B")
KRT_signature <- c("KRT6A","KRT6B","KRT79","KRT14","KRT80","KRT23","KRT16","KRT15")
hallmark_apoptosis <- c("BIRC3","BIK","CDKN1A","DDIT3","IL18","TNFRSF12A","HSPB1","ATF3","HMOX1","IER3","TGFB2","TXNIP","ANXA1")


NRG1_signalling <- c("NFIL3","ID1","TIPARP","MIR22HG","AVPI1","ATF3","CLCF1","PPP1R15A","GADD45A","CYP24A1","ZFP36","USP36","CYR61","TRIB1","MAP3K8","BHLHE40","ZNF331","GEM","TNFRSF11B","PLAUR","MAFF","NR4A2","DLX2","IER3","SGK1","FOS","REL","DUSP4","CTH","FOSL1","IL6R","SPRR1B","CCNL1","IRS2","NR4A3","NR4A1","DDIT3")

MYB_targets <- c("HSPB1","MAP1LC3B","S100A7","CDKN1A","HMOX1","SERPINB2","RRAD","ANXA1","KRT17","ANXA3","MYL12A","TUBA4A","KRT16","IL8","ARID5B","IGFBP7","KRT14","AHNAK2","KRT6A","PURA","NNMT")

SMID_LUMINAL_DN <- c("ALDH1A3","NNMT","ID4","KYNU","KRT6A","NOV","KRT23","DSC3","TNFRSF11B","CLIC3","IL8","ANXA1","ANXA3","BIRC3","TMEM45A","KLK7","DSC2","KRT17","ADM","KLK6","KLK10","ID1","KRT16","S100A7","KRT14","TAGLN","PERP")

TGFB_KEGG <- c('CHRD', 'NOG', 'NBL1', 'MICOS10-NBL1', 'GREM1', 'GREM2'
               , 'THBS1', 'DCN', 'FMOD', 'LEFTY1', 'LEFTY2', 'FST', 'BMP2', 
               'BMP4', 'BMP6', 'INHBB', 'BMP5', 'BMP7', 'BMP8B', 'BMP8A', 
               'GDF5', 'GDF6', 'GDF7', 'AMH', 'THSD4', 'FBN1', 'LTBP1', 'TGFB1',
               'TGFB2', 'TGFB3', 'INHBA', 'INHBC', 'INHBE', 'NODAL', 'NEO1', 'HJV',
               'BMPR1A', 'BMPR1B', 'ACVR1', 'BMPR2', 'ACVR2A', 'RGMA', 'RGMB',
               'AMHR2', 'TGFBR1', 'TGFBR2', 'ACVR1B', 'ACVR2B', 'ACVR1C', 'BAMBI',
               'SMAD1', 'SMAD5', 'SMAD9', 'SMAD2', 'SMAD3', 'SMAD4', 'SMAD6',
               'SMAD7', 'SMURF1', 'SMURF2', 'ZFYVE9', 'ZFYVE16', 'HAMP', 'ID1',
               'ID2', 'ID3', 'ID4', 'RBL1', 'E2F4', 'E2F5', 'TFDP1', 'CREBBP',
               'EP300', 'SP1', 'TGIF1', 'TGIF2', 'MYC', 'CDKN2B', 'PITX2', 'RBX1',
               'CUL1', 'SKP1', 'MAPK1', 'MAPK3', 'IFNG', 'TNF', 'RHOA', 'ROCK1', 
               'PPP2R1B', 'PPP2R1A', 'PPP2CA', 'PPP2CB', 'RPS6KB1', 'RPS6KB2')

###### ANNOTATION ######
organism <- "hg38"
MSigDBFile1 <- file.path(here(),"annotation","hg38.MSIG.gs.rda")
MSigDBFile2 <-file.path(here(),"annotation","hg38.MSIG.ls.rda")
# load(file=file.path(RDataSupdir,paste("Supervised_res_object_edgeR.Rdata",sep="")))
load(MSigDBFile1)
load(MSigDBFile2)
MSIG.ls = hg38.MSIG.ls
MSIG.gs = hg38.MSIG.gs

###### OPTIONS ######
options(stringsAsFactors=FALSE, width=180)
