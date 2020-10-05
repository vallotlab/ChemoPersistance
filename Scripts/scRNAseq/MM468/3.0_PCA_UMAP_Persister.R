library(corrplot)
library(seriation)
library(heatmap.plus)
library(scatterplot3d)
library(here)
library(ChromSCape)
library(devtools)
library(dplyr)
library(DropletUtils)
library(irlba)
library(corrplot)
library(ConsensusClusterPlus)
library(geco.unsupervised)
library(scatterplot3d)
library(scater)
library(Rtsne)
library(ccRemover)
library(colorRamps)
library(geco.supervised)
library(viridis)
library(colorRamps)
library(RColorBrewer)
library(scTools)
library(edgeR)
library(ggplot2)
library(RColorBrewer)
library(genefilter)
library(xtable)
library(WriteXLS)
library(data.table)
library(stringr)
library(limma)
library(edgeR)
library(monocle3)
library(dplyr)
library(WriteXLS)
library(Seurat)
library(gplots)

#Geco packages
library(geco.supervised)
library(geco.RNAseq)
library(geco.utils)
library(geco.visu)
library(scTools)
library(geco.unsupervised)

# Directories -------------------------------------------------------------
maindir= here()
resdir <- file.path(maindir,"output","scRNAseq","MM468","Persister")
QCdir <- file.path(maindir,"output","scRNAseq","QC")

resdir <- file.path(resdir, "Unsupervised") ; if(!file.exists(resdir)){dir.create(resdir)}
resdir_UMAP <- file.path(resdir,"UMAP") ; if(!file.exists(resdir_UMAP)){dir.create(resdir_UMAP)}
resdir_boxplots <- file.path(resdir, "boxplots") ; if(!file.exists(resdir_boxplots)){dir.create(resdir_boxplots)}
resdir_heatmaps = file.path(resdir,"Heatmaps"); if(!dir.exists(resdir_heatmaps)) dir.create(resdir_heatmaps)
resdir_sub = file.path(resdir,"sub500"); if(!dir.exists(resdir_sub)) dir.create(resdir_sub)
RDatadir <- file.path(resdir,"RData") ; if(!file.exists(RDatadir)){dir.create(RDatadir)}

source(file.path(maindir,"Scripts","global_var.R"))

################################################################################
# Parameters --------------------------------------------------------------
################################################################################

## packages Seurat or scater attach an enormous amount of bystander packages, they cannot be loaded in the same session without crashing the DLLs.
##Hierarchical clustering
distHC <- c("distPearson","distCosine","euclidean","maximum","manhattan","canberra","binary","minkowski")[1]
methHC <- c("ward","ward.D","ward.D2","single","complete","average")[2]

##ConsClust
repsCC <- 10
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

################################################################################################
# PCA, UMAP and louvain clustering with Monocle 3  on all cells for the paper------------------
###############################################################################################
load(file.path(RDatadir,"MM468.RData"))

annot_sel <- annot
samples_persister = c(
  "MM468_initial","MM468_5FU3_day50","MM468_5FU3_day77","MM468_5FU3_day202",
  "MM468_5FU5_day67","MM468_5FU5_day171","MM468_5FU6_day33","MM468_5FU6_day214"
)

annot_int <- annot_sel[annot_sel$sample_id %in% samples_persister & annot_sel$doublet_score<10000,]

annot_int$sample_id <- as.character(annot_int$sample_id)
rownames(annot_int) <- annot_int$cell_id

cds <- new_cell_data_set(Signal[,row.names(annot_int)],
                         cell_metadata = annot_int,
                         gene_metadata  = gene_metadata)

#cds <- preprocess_cds(cds,method='PCA', norm_method='log',num_dim=50)
cds <- preprocess_cds(cds,method='PCA', norm_method='size_only',num_dim=50)

cds <- reduce_dimension(cds, reduction_method = 'UMAP')
cds <-  cluster_cells(cds)

 umap_res <- cds@int_colData$reducedDims[[2]] # with R3.6.2 and higher for SummarizedExperiment, otherwise cds@reducedDims[[2]]
 pca_object <- cds@int_colData$reducedDims[[1]]

annot_int$louvain_partition <- cds@clusters[[1]]$partitions
annot_int$louvain_partition <- paste0("C",annot_int$louvain_partition)

save(umap_res,file=file.path(RDatadir,"umap_persister.RData"))
save(pca_object,file=file.path(RDatadir,"pca_persister.RData"))

NormCounts <- t(t(exprs(cds)) /  pData(cds)[, 'Size_Factor'])
persister_LogCounts <- Matrix:: log(NormCounts+1,2)
rm(NormCounts)

save(persister_LogCounts,file=file.path(RDatadir,"persister_LogCounts.RData"))
save(annot_int,gene_metadata,file=file.path(RDatadir,"gene_cell_annot_persister.RData"))
# rownames(NormCounts) <- gene_metadata$Symbol
# colnames(NormCounts) <- annot_sel$cell_id



# ################################################################################
# # Cell cycle scoring with Seurat  ##############################################
# ################################################################################
# 
# library(Seurat)
# 
# # A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# # segregate this list into markers of G2/M phase and markers of S phase
# s.genes <- cc.genes$s.genes
# g2m.genes <- cc.genes$g2m.genes
# 
# # Create our Seurat object and complete the initalization steps
# rownames(Signal) <- gene_metadata$Symbol
# 
# marrow <- CreateSeuratObject(counts = Signal[,row.names(annot_sel)],meta.data=annot_sel)
# marrow <- NormalizeData(marrow)
# marrow <- FindVariableFeatures(marrow, selection.method = "vst")
# marrow <- ScaleData(marrow, features = rownames(marrow))
# marrow <- RunPCA(marrow, features = VariableFeatures(marrow), ndims.print = 1:10, nfeatures.print = 10)
# DimPlot(marrow, reduction = "pca")
# DimHeatmap(marrow, dims = 1:3, cells = 1000, balanced = TRUE)
# ElbowPlot(marrow)
# marrow <- FindNeighbors(marrow, dims = 1:20)
# marrow  <- FindClusters(marrow, resolution = 0.5)
# marrow <- RunUMAP(marrow, dims = 1:20)
# 
# #DimPlot(marrow, reduction = "umap")
# 
# marrow <- CellCycleScoring(marrow, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
# #RidgePlot(marrow, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)
# 
# annot_seurat <- data.frame (cell_id=as.character(names(marrow$Phase)),cell_cycle=as.character((marrow$Phase)))
# 
# annot_sel$cell_cycle <- marrow$Phase[match(annot_sel$cell_id,annot_seurat$cell_id)]
# 
# save(LogCounts,gene_metadata,annot_sel,file=file.path(RDatadir,"forAnalysis.RData"))
# 
################################################################################################
###### Parameters for à façon plots of results on all cells--  #################################
################################################################################################

load("pca_persister.RData")
load("persister_LogCounts.RData")
load("gene_cell_annot_persister.RData")
load("umap_persister.RData")
gc()

#select cells of interest for figure
annot_int$sample_id <- as.character(annot_int$sample_id)

annotCol <- c("sample_id","total_features","cons_BC_lenti","doublet_score","louvain_partition","rRNA","BMP4","VIM","ABCC4","MKI67","KRT14","KLK10","NNMT","TGFB1","TGFBR1","INHBA","INHBB","TGFBR2","TGFBR3")

annotText <- "sample_id"
hcText <- "sample_id"  ## column used in hierarchical clustering # change names from Sample_x to actual sample name

#Common palette for all MM468 experiments
control <- c("#E0E0E0","#BDBDBD","#757575")
persister_color <- c("#DCEDC8","#9CCC65","#4CAF50","#009688")
treated_persister_color <- c("#2196F3", "#4DD0E1" )
treated_GSKJ4_color <- c("#C2185B", "#EC407A" )
res_color <- c("#FFEB3B","#FFC107","#FF9800","#FF5722")
color_MM468 <- c(control[2],persister_color,res_color[c(1:3)], treated_persister_color)

#Adding annotation for expression of genes of interest
annot_int$KRT14 <- persister_LogCounts["KRT14",annot_int$cell_id]; gc()
annot_int$KRT17 <- persister_LogCounts["KRT17",annot_int$cell_id]; 
annot_int$KLK10 <- persister_LogCounts["KLK10",annot_int$cell_id]; gc()
annot_int$MKI67 <- persister_LogCounts["MKI67",annot_int$cell_id]; gc()
annot_int$NNMT <- persister_LogCounts["NNMT",annot_int$cell_id]; gc()
annot_int$TGFB1 <- persister_LogCounts["TGFB1",annot_int$cell_id]; gc()
annot_int$TGFBR1 <- persister_LogCounts["TGFBR1",annot_int$cell_id]; gc()
annot_int$TGFBR3 <- persister_LogCounts["TGFBR3",annot_int$cell_id]; gc()
annot_int$TGFBR2 <- persister_LogCounts["TGFBR2",annot_int$cell_id]; gc()
annot_int$INHBA <- persister_LogCounts["INHBA",annot_int$cell_id]; gc()
annot_int$INHBB <- persister_LogCounts["INHBB",annot_int$cell_id]; gc()
annot_int$BMP4 <- persister_LogCounts["BMP4",annot_int$cell_id]; gc()
annot_int$VIM <- persister_LogCounts["VIM",annot_int$cell_id]; gc()
annot_int$ABCC4 <- persister_LogCounts["ABCC4",annot_int$cell_id]; gc()

ribo <- grepl("^RP[SL]",gene_metadata$Symbol)
annot_int$rRNA <- apply(persister_LogCounts[ribo,],2,mean); gc()

#Generate color annotation & legend
anocol <- geco.annotToCol4(annotS=annot_int[,annotCol],annotT=annot_int,plotLegend=T,plotLegendFile=file.path(resdir,"Annotation_legends.pdf"), categCol=NULL)

#Manual corrections of colors
indice <- which(annotCol=="sample_id")
corres <- data.frame(Sample=c("MM468_initial","MM468_5FU3_day50","MM468_5FU3_day77","MM468_5FU5_day67","MM468_5FU6_day33","MM468_5FU6_day214","MM468_5FU3_day202","MM468_5FU5_day171","MM468_5FU6_persister_day33","MM468_persister_day33"),color=color_MM468)
anocol[,indice] <- as.character(corres$color[match(annot_int$sample_id,corres$Sample)])

indice <- which(annotCol=="lineage_barcode")
anocol[is.na(annot_int$lineage_barcode),indice] <- "white"

corres_cell_cycle <- data.frame(phase=c("G1","G2M","S"),color=c("#ff80b0","#9399ff","#a9fffd"))
indice <- which(annotCol=="cell_cycle")
anocol[,indice] <- as.character(corres_cell_cycle$color[match(annot_int$cell_cycle,corres_cell_cycle$phase)])

corres_cluster <- data.frame(phase=c("C6","C2","C3","C4","C5"),color=c("#b9cced","#6a8caf","#438a5e","#e2979c","#ff9c71"))
indice <- which(annotCol=="clustCol5")
# anocol[,indice] <- as.character(corres_cluster$color[match(annot_sel$clustCol5,corres_cluster$phase)])

png(file.path(resdir,"Legend_cluster_scRNAseq.png"), height=2000,width=1500,res=300)
par(mfrow=c(5,2))
barplot(rep(1,5),col=as.character(corres_cluster$color), cex.names  =0.8,horiz=F,names.arg=c("C1","C2","C3","C4","C5"),las=2)
dev.off()

#Manual corrections of colors

png(file.path(resdir,"Legend_cell_cycle.png"), height=2000,width=1500,res=300)
par(mfrow=c(4,2))
barplot(rep(1,3),col=c("pink2","slateblue3","slategray"), cex.names  =1,horiz=F,names.arg=c("G1","G2M","S"),las=2)
dev.off()

png(file.path(resdir,"Legend_sample_id.png"), height=2000,width=1500,res=300)
par(mfrow=c(5,2))
barplot(rep(1,10),col=as.character(corres$color), cex.names  =0.8,horiz=F,las=2,names.arg=as.character(corres$Sample))
dev.off()

save(anocol,file=file.path(RDatadir,"anocol.RData"))


####################################################################################
###### Plot UMAP results ###########################################################
####################################################################################

for(j in 1:ncol(anocol))
{
  png(file.path(resdir,paste0("UMAP_",colnames(anocol)[j],".png")), height=1350,width=1200,res=300) 
  plot((umap_res[annot_int$cell_id,]), col=alpha(anocol[,j],0.2),pch=20,cex=0.4,main=paste0(colnames(anocol)[j]," perplexity=30"),xlab="component 1",ylab="component 2")
  dev.off()
}
indice <- which(annotCol=="cons_BC_lenti")

png(file.path(resdir,paste0("UMAP_",colnames(anocol)[indice],".png")), height=1150,width=950,res=350) 
plot((umap_res[!is.na(annot_int$cons_BC_lenti),]),
     col=alpha(anocol[!is.na(annot_int$cons_BC_lenti),indice],0.5),
     pch=20,cex=0.4,
     main=paste0(colnames(anocol)[indice]," perplexity=30"),
     xlab="component 1",ylab="component 2",ylim=c(min(umap_res[,2]),
                                                  max(umap_res[,2])))
dev.off()


#count number of barcodes detected per expression group: C10, C2, C3, C4, C6, C8
stats_detection <- annot_int %>% mutate(sample_id=factor(sample_id,levels=c("MM468_initial","MM468_5FU6_day33","MM468_5FU3_day50", "MM468_5FU5_day67","MM468_5FU3_day77","MM468_persister_day33","MM468_5FU6_persister_day33")))%>% group_by(sample_id) %>% summarise(cells=length(cell_id),barcode=length(which(!is.na(cons_BC_lenti))),fraction=length(which(!is.na(cons_BC_lenti)))/length(cell_id))
stats_barcode <- annot_int[!is.na(annot_int$cons_BC_lenti),] %>% mutate(sample_id=factor(sample_id,levels=c("MM468_initial","MM468_5FU6_day33","MM468_5FU3_day50", "MM468_5FU5_day67","MM468_5FU3_day77","MM468_5FU6_persister_day33","MM468_persister_day33")))%>% group_by(sample_id,louvain_partition) %>% summarise(diff_barcode=length(unique(cons_BC_lenti)),barcode=length(cons_BC_lenti),diversity=length(unique(cons_BC_lenti))/length(cons_BC_lenti)) %>% filter(barcode>15) %>% mutate(louvain_partition=factor(louvain_partition,levels=c("C3","C1","C2","C4"))) 

png(file.path(resdir,"Lineage_barcode_detection_scRNA.png"), height=1500,width=1500,res=300)
ggplot(stats_detection, aes(x=sample_id,y=fraction,fill=sample_id,label=barcode)) + geom_bar(position='dodge',stat='identity') +theme_classic() +scale_fill_manual(values=as.character(corres$color[c(1,5,2,4,10,9)])) + geom_text(position = position_dodge(width = 0.8), angle = 90 )
dev.off()

png(file.path(resdir,"Lineage_diversity.png"), height=800,width=1500,res=300)
ggplot(stats_barcode, aes(x=louvain_partition,y=diversity, fill=sample_id)) + geom_bar(stat='identity',position = position_dodge2(width = 0.9, preserve = "single")) +theme_classic() +
  scale_fill_manual(values=as.character(corres$color[c(1,5,2,4,10,9)]))
dev.off()

################################################################################################################
# Hierarchical clustering to group samples in complement to Louvain clustering on a subset of cells ############
################################################################################################################
#take a subset, do ONCE
anocol <- as.data.frame(anocol)
rownames(anocol) <- annot_int$cell_id
# 
annot_subset <- annot_int %>% group_by(sample_id) %>% sample_n(500)
annot_subset <- as.data.frame(annot_subset)
annot_subset$cell_id <- as.character(annot_subset$cell_id)
anocol_subset <- anocol[annot_subset$cell_id,]
# 
# head(anocol_subset)
# anocol_subset <- anocol_subset[,-c(4,6,7)]
# head(anocol_subset)

#Load if second time:

mati <- t(as.data.frame(pca_object)[annot_subset$cell_id,])
hc_PCA <- hclust(as.dist(1 - cor(mati)), method = "ward.D")
mat.so.cor <- mati[,hc_PCA$order]




#add the day to the annot table
# corres_day <- data.frame(Sample=c("MM468_DMSO3_day50","MM468_DMSO5_day67","MM468_5FU3_day50","MM468_5FU3_day77","MM468_5FU5_day67","MM468_5FU6_day33","MM468_5FU3_day202","MM468_5FU5_day171","MM468_5FU6_persister_day33","MM468_persister_day33"),
#                          day=c(50,67,50,77,67,33,202,171,33,33))
# annot_subset$day <- as.integer(corres_day$day[match(annot_subset$sample_id,corres$Sample)])
# 

png(file.path(resdir,"Clustering_correlation_matrix_PCA_subset.png"), height=1500,width=1500,res=300)
geco.hclustAnnotHeatmapPlot(x=cor(mat.so.cor),
                            hc=hc_PCA,
                            hmColors=corColors,
                            anocol=as.matrix(anocol_subset[hc_PCA$order,]),#[,ncol(cc.col):1]
                            xpos=c(0.15,0.9,0.164,0.885),
                            ypos=c(0.1,0.5,0.5,0.6,0.62,0.95),
                            dendro.cex=0.1,
                            xlab.cex=0.2,
                            hmRowNames=FALSE,
                            hmRowNames.cex=0.01
)
dev.off()

#Consensus clustering, do once or load again 
maxKCC <- 10
repsCC <- 100
gc()
consclust <- ConsensusClusterPlus(mati, maxK=maxKCC,plot="png",reps=repsCC, pItem=pItemCC,pFeature=pFeatureCC,title="Consensus clustering", clusterAlg="hc",distance=distCC,innerLinkage=innerLinkageCC, finalLinkage=finalLinkageCC,seed=3.14)
save(consclust,file=file.path(RDatadir,"Consensus_clustering.RData"))

load(file.path(RDatadir,"Consensus_clustering.RData"))
#Add to annotation consensus results
nclust <- 6
clustCol <- "ExpressionGroup"
cc <- consclust[[nclust]]$consensusClass
cc. <- lapply(unique(cc),function(z) names(which(cc==z)))
mat.cc <- geco.groupMat(mati,margin=1,groups=cc.,method="mean")
hcc <- hclust(distPearson(t(mat.cc)),method="ward.D")	
annot_subset[,clustCol] <- paste("C",match(cc,hcc$order),sep="")

# #higher quality consensus graph
# nclust <- 6
# consclust.mat <- consclust[[nclust]]$ml
# hc = hclust(as.dist(1 - consclust.mat), method = finalLinkageCC)
# consclust.mat <- consclust.mat[hc$order,]
# sampnames <- names(consclust[[nclust]]$consensusClass)[consclust[[nclust]]$consensusTree$order]
# png(file.path(resdir,paste("Consensus_matrix_k=",nclust,".png",sep="")))
# heatmap(consclust.mat, Colv=as.dendrogram(hc), Rowv =NA, symm = FALSE, scale = "none", col=colorRampPalette(c("white", "blue"))(100), na.rm = TRUE, labRow = sampnames, labCol = F, mar = c(5, 5),
#         main = paste("consensus matrix k=",nclust, sep = ""))
# dev.off()


# pdf(file.path(resdir,"ICL.pdf"))
# icl <- calcICL(consclust)
# dev.off()
# save(icl,file=file.path(RDatadir,"ICL.RData"))

#Add to annotation hierarchical clustering and consensus results
# nclust <- 6
# #annot_subset <- as.data.frame(annot_subset)
# annot_subset$clustCol6 <- paste("C",match(cutree(hc_PCA,nclust),unique(cutree(hc_PCA,nclust)[hc_PCA$order])),sep="")

#Generate novel color codes
# annotCol <- c("sample_id","ExpressionGroup","total_features","cell_cycle","rRNA","MKI67","KRT17","VIM","TGFB1")
annotCol <- c("sample_id","ExpressionGroup","total_features","rRNA","MKI67","VIM","TGFB1")
anocol_subset <- geco.annotToCol4(annotS=annot_subset[,annotCol],annotT=annot_subset,plotLegend=T,
                                  plotLegendFile=file.path(resdir,"Annotation_legends.pdf"),
                                  categCol=NULL)

#Manual corrections of colors
anocol_subset[,1] <- as.character(corres$color[match(annot_subset$sample_id,corres$Sample)])
indice <- which(annotCol=="cell_cycle")
anocol_subset[,indice] <- as.character(corres_cell_cycle$color[match(annot_subset$cell_cycle,corres_cell_cycle$phase)])

# save(annot_subset2,file=file.path(RDatadir,"annot_subset.RData"))
# save(anocol_subset2,file=file.path(RDatadir,"anocol_subset2.RData"))

png(file.path(resdir,"Clustering_correlation_matrix_PCA_subset_clusterinfo_2.png"), height=1500,width=1500,res=300)
geco.hclustAnnotHeatmapPlot(x=cor(mat.so.cor),
                            hc=hc_PCA,
                            hmColors=corColors,
                            anocol=as.matrix(anocol_subset[hc_PCA$order,c(1:3,6)]),#[,ncol(cc.col):1]
                            xpos=c(0.15,0.9,0.164,0.885),
                            ypos=c(0.1,0.5,0.5,0.6,0.62,0.95),
                            dendro.cex=0.1,
                            xlab.cex=0.2,
                            hmRowNames=FALSE,
                            hmRowNames.cex=0.01
)
dev.off()

load("persister_LogCounts.RData")
LogCounts = persister_LogCounts; gc()
subset_LogCounts <- LogCounts[,annot_subset$cell_id]; gc()
save(subset_LogCounts,hc_PCA,anocol_subset,annot_subset,gene_metadata,file=file.path(RDatadir,"subset_analysis.RData"))

for(j in 1:ncol(anocol_subset))
{
  png(file.path(resdir,paste0("UMAP_subset_",colnames(anocol_subset)[j],".png")), height=1500,width=1500,res=300) 
  
  if(class(annot_subset[1,colnames(anocol_subset)[j]])=="numeric"){
    plot((umap_res[annot_subset$cell_id,]), col=alpha(anocol_subset[,j],0.3),pch=20,cex=0.6,
         main=paste0(colnames(anocol_subset)[j]," min=",round(min(annot_subset[,colnames(anocol_subset)[j]]),digits=3)," max=",round(max(annot_subset[,colnames(anocol_subset)[j]]),digits=3)),
         xlab="component 1",ylab="component 2")} else {
           plot((umap_res[annot_subset$cell_id,]), col=alpha(anocol_subset[,j],0.3),pch=20,cex=0.6,
                main=paste0(colnames(anocol_subset)[j]),
                xlab="component 1",ylab="component 2")
         }
  
  plot((umap_res[annot_subset$cell_id,]), col=alpha(anocol_subset[,j],0.3),pch=20,cex=0.4,main=paste0(colnames(anocol_subset)[j]," perplexity=30"),xlab="component 1",ylab="component 2")
  dev.off()
}


#Plot repartition of Expression groups with time
corres_day <- data.frame(Sample=c("MM468_DMSO3_day50","MM468_DMSO5_day67","MM468_5FU3_day50","MM468_5FU3_day77",
                                  "MM468_5FU5_day67","MM468_5FU6_day33","MM468_5FU3_day202","MM468_5FU5_day171","MM468_5FU6_persister_day33","MM468_persister_day33"),
                          day=c(50,67,50,77,67,33,202,171,33,33))
annot_subset$day <- as.integer(corres_day$day[match(annot_subset$sample_id,corres$Sample)])

annot_subset$col <- anocol_subset[,2]
table(annot_subset$ExpressionGroup,annot_subset$col)
test <- annot_subset %>% mutate(ExpressionGroup=factor(ExpressionGroup,levels = unique(annot_subset$ExpressionGroup))) %>%
  mutate(sample_id=factor(sample_id,levels=c("MM468_5FU3_day50","MM468_5FU5_day67","MM468_5FU6_day33","MM468_5FU6_persister_day33",
                                             "MM468_initial","MM468_persister_day33"))) %>% 
  group_by(sample_id,ExpressionGroup) %>% summarise(freq=n()) 


png(file.path(resdir,"ExpressionGroup_time.png"),width=1500,height=1500,res=300)
ggplot(test) + geom_col(aes(y=freq,x=sample_id,fill=ExpressionGroup)) +scale_fill_manual(values=c("slateblue4","honeydew1","paleturquoise4","gold","mediumseagreen","mediumorchid")) + theme_classic() + theme(axis.text.x = element_text(angle = 90)) 
dev.off()

################################################################################
###Generate heatmaps with top varying genes 
################################################################################
#work on subset_LogCounts 

NbGenes <- 100
exp.sd <- apply(subset_LogCounts,1,sd)
sel <- order(exp.sd,decreasing=T)[1:NbGenes]
mat <- subset_LogCounts[sel,];dim(mat)
tmat <- t(mat)

distHC <- c("distPearson","distCosine","euclidean","maximum","manhattan","canberra","binary","minkowski")[3]
methHC <- c("ward","ward.D","ward.D2","single","complete","average")[3]
if(distHC=="distPearson")	d <- distPearson(tmat)
if(distHC=="distCosine")	d <- distCosine(tmat)
if(distHC %in% c("euclidean","maximum","manhattan","canberra","binary","minkowski"))	d <- dist(tmat,distHC)
hc <- hclust(d,method=methHC)
nclust <- 6
annot_subset$clustCol_reorder <- paste("C",match(cutree(hc,nclust),unique(cutree(hc,nclust)[hc$order])),sep="")
test <- c(which(annot_subset$clustCol_reorder %in% c("C4")), which(annot_subset$clustCol_reorder %in% "C6"),
          which(annot_subset$clustCol_reorder %in% "C5"),which(annot_subset$clustCol_reorder %in% "C2"),
          which(annot_subset$clustCol_reorder %in% "C3"),which(annot_subset$clustCol_reorder %in% "C1"))

hc_reorder <- rotate(hc,test)
mat.so <- as.matrix(mat[,hc_reorder$order])

if(chRangeHM){
  for(i in 1:nrow(mat.so)){
    mat.so[i,] <- geco.changeRange(mat.so[i,],newmin=0,newmax=1)
  }
}
rowClust <- hclust(as.dist(1 - cor(t(mat.so))), method = "ward.D")

png(file.path(resdir,paste0("Re_Clustering_MostVarying_",NbGenes,"_",distHC,"_",methHC,".png")), height=4000,width=4000,res=300)
geco.hclustAnnotHeatmapPlot(x=(mat.so[rowClust$order,]),
                            hc=hc_reorder,
                            hmColors=hmColors,
                            anocol=anocol_subset[hc_reorder$order,c(1:2,6)],#[,ncol(cc.col):1]
                            xpos=c(0.15,0.9,0.164,0.885),
                            ypos=c(0.1,0.5,0.5,0.6,0.62,0.95),
                            dendro.cex=0.01,
                            xlab.cex=1,
                            hmRowNames=TRUE,
                            hmRowNames.cex=0.2
)
dev.off()

if(chRangeHM){
  for(i in 1:nrow(mat.so)){
    mat.so[i,] <- geco.changeRange(mat.so[i,],newmin=0,newmax=1)
  }
}
rowClust <- hclust(as.dist(1 - cor(t(mat.so))), method = "ward.D")

png(file.path(resdir,paste0("Reorder_Clustering_MostVarying_",NbGenes,"_",distHC,"_",methHC,".png")), height=4000,width=4000,res=300)

geco.hclustAnnotHeatmapPlot(x=(mat.so[rowClust$order,]),
                            hc=hc_reorder,
                            hmColors=hmColors,
                            anocol=anocol_subset[hc_reorder$order,c(1:2,6)],#[,ncol(cc.col):1]
                            xpos=c(0.15,0.9,0.164,0.885),
                            ypos=c(0.1,0.5,0.5,0.6,0.62,0.95),
                            dendro.cex=0.01,
                            xlab.cex=1,
                            hmRowNames=TRUE,
                            hmRowNames.cex=0.2
)
dev.off()


################################################################################
###Generate heatmaps with gene lists 
################################################################################
#load "subset_analysis.RData" if starting here
#import my.res et gene lists from diff analysis, logCounts, take top 10
library(gplots)
setwd(RDatadir)

load("../../Supervised/RData/Supervised_res_object_edgeR.Rdata")
resdir <- "../../../Results_MM468/Supervised";if(!file.exists(resdir)){dir.create(resdir)}

RDataSupdir <-  file.path(resdir,"RData");if(!file.exists(RDataSupdir)){dir.create(RDataSupdir)}

# load(file=file.path(RDataSupdir,"Overexpressed_persisterall_vs_DMSO.RData"))
load(file=file.path(RDataSupdir,"Overexpressed_persisterall_vs_DMSO.RData"))
how_many_top <- 50
significPathway <- Overexpressed[Overexpressed$Class %in% c("c2_curated","c5_GO","hallmark"),]

significPathway$Deregulated_genes <- as.character(significPathway$Deregulated_genes)
significPathway <- significPathway[1:how_many_top,]
significPathway_breast <- 	rbind(significPathway[grep(pattern = "BREAST", significPathway$Gene.Set),],significPathway[grep(pattern = "MAMMARY", significPathway$Gene.Set),],significPathway[grep(pattern = "HALLMARK", significPathway$Gene.Set),])

#Calculate average expression for every cells on genes of each top pathways
pathway_mat <- matrix(0,nrow=dim(significPathway)[1],ncol=length(annot_subset$sample_id))		
pathway_names <- c()		

for(i in 1:dim(significPathway)[1]){
  
  gene_list <- unlist(strsplit(significPathway$Deregulated_genes[i],";"))
  mat <- subset_LogCounts[gene_metadata$Symbol %in% gene_list,]
  pathway_mat[i,] <- apply(mat,2,mean) 
  pathway_names[i] <-significPathway$Gene.Set[i]
}		


row.names(pathway_mat) <- pathway_names
tmat <- t(pathway_mat)

distHC <- c("distPearson","distCosine","euclidean","maximum","manhattan","canberra","binary","minkowski")[3]
methHC <- c("ward","ward.D","ward.D2","single","complete","average")[5]
if(distHC=="distPearson")	d <- distPearson(tmat)
if(distHC=="distCosine")	d <- distCosine(tmat)
if(distHC %in% c("euclidean","maximum","manhattan","canberra","binary","minkowski"))	d <- dist(tmat,distHC)
hc <- hclust(d,method=methHC)

mat.so <- pathway_mat[,hc$order]

if(chRangeHM){
  for(i in 1:nrow(mat.so)){
    mat.so[i,] <- geco.changeRange(mat.so[i,],newmin=0,newmax=1)
  }
}

#rowClust <- hclust(as.dist(1 - cor(t(mat.so))), method = "ward.D")
rowClust <- hclust(dist((mat.so)), method = "ward.D")

png(file.path(resdir,paste0("Clustering_",how_many_top,"Pathways_persistersAll_",distHC,"_",methHC,".png")), 
    height=3000,width=2000,res=300)
par(oma=c(2,5,3,7))
geco.hclustAnnotHeatmapPlot(x=(mat.so[rowClust$order,]),
                            hc=hc,
                            hmColors=hmColors,
                            anocol=anocol_subset[hc$order,c(1:2,6)],#[,ncol(cc.col):1]
                            xpos=c(0.15,0.9,0.164,0.885),
                            ypos=c(0.1,0.5,0.5,0.6,0.62,0.95),
                            dendro.cex=0.4,
                            xlab.cex=0.4,
                            hmRowNames=TRUE,
                            hmRowNames.cex=0.4
)
dev.off()




################################################################################
#Plot distribution of correlation scores and vignettes
################################################################################

interm <- data.frame(cor(mati))
interm$first_cell <- rownames(interm)
all_corr <- tidyr::pivot_longer(interm,cols=1:3000,names_to="other_cell",values_to="Corr")
all_corr <- dplyr::distinct(all_corr)
all_corr$first_cell_sample_id <- annot_subset$sample_id[match(all_corr$first_cell,annot_subset$cell_id)]
all_corr$other_cell_sample_id <- annot_subset$sample_id[match(all_corr$other_cell,annot_subset$cell_id)]

#Plot intra-sample corr_scores
intra_corr <- all_corr[all_corr$first_cell_sample_id==all_corr$other_cell_sample_id,]
intra_corr <- as.data.frame(intra_corr)

for(Sample in unique(annot_subset$sample_id)){
  MAT <- intra_corr[intra_corr$first_cell_sample_id==Sample,]
  MAT <- tidyr::pivot_wider(MAT[,c(1:3)],names_from="other_cell",values_from="Corr")
  
  png(file.path(resdir,paste0("Pearson_Correlation_scores_",Sample,".png")), height=1000,width=1000,res=300)
  image((as.matrix(MAT[,-1])),col=rev(viridis(100)))
  dev.off()
  png(file.path(resdir,paste0("Distrib_Pearson_Correlation_scores_",Sample,".png")), height=1000,width=1000,res=300)
  hist(as.matrix(MAT[,-1]),breaks=100,xlim=c(-1,1),main=paste0("Pearson's Correlation Scores ",Sample),xlab="")
  dev.off()
  
}



all_corr$first_cell_cluster_id <- annot_subset$louvain_partition[match(all_corr$first_cell,annot_subset$cell_id)]
all_corr$other_cell_cluster_id <- annot_subset$louvain_partition[match(all_corr$other_cell,annot_subset$cell_id)]

#Plot intra-sample corr_scores
intra_corr_cluster <- all_corr[all_corr$first_cell_cluster_id==all_corr$other_cell_cluster_id,]
intra_corr_cluster <- as.data.frame(intra_corr_cluster)

for(Cluster in unique(annot_subset$louvain_partition)){
  MAT <- intra_corr_cluster[intra_corr_cluster$first_cell_cluster_id==Cluster,]
  MAT <- tidyr::pivot_wider(MAT[,c(1:3)],names_from="other_cell",values_from="Corr")
  
  png(file.path(resdir,paste0("Pearson_Correlation_scores_",Cluster,".png")), height=1000,width=1000,res=300)
  image((as.matrix(MAT[,-1])),col=rev(viridis(100)))
  dev.off()
  png(file.path(resdir,paste0("Distrib_Pearson_Correlation_scores_",Cluster,".png")), height=1000,width=1000,res=300)
  hist(as.matrix(MAT[,-1]),breaks=100,xlim=c(-1,1),main=paste0("Pearson's Correlation Scores ",Sample),xlab="")
  dev.off()
  
}

refill <- corres$color[match(names(table(intra_corr$other_cell_sample_id)), corres$Sample)]

sel <- (intra_corr$first_cell %in% annot_int$cell_id[annot_int$louvain_partition=="C2"] | 
          intra_corr$first_cell %in% annot_int$cell_id[annot_int$louvain_partition=="C3"])

png(file.path(resdir,paste0("DistribPearsonCorrelationScores.png")),height=1000,width=1500,res=300)
ggplot(intra_corr[sel,],aes(x=first_cell_sample_id,y=Corr, fill=first_cell_sample_id)) + 
  geom_violin(alpha=0.8) + theme_classic() + scale_fill_manual(values=as.character(corres$color[c(2:4,9,1)])) +
  stat_summary(fun.y=median, geom="point", size=2, color="black")
dev.off()
#ggplot(intra_corr,aes(x=Corr, fill=first_cell_sample_id)) + geom_histogram(alpha=0.25) + theme_classic(base_size = 25) + scale_fill_manual(values=refill[-c(3,4,10)])

################################################################################
################ Study PC content --  ##########################
################################################################################

# pca_gene_loading <- cds@preprocess_aux$gene_loadings
# 
# variances <- round(cds@preprocess_aux$prop_var_expl,digits=4)
# barplot(variances, names.arg = names(variances), cex.names = 1, ylab = "Proportion of Variance", xlab="Principal components", las=2,col="royalblue",border="darkblue")
# 
# #top 10
# top <- 30
# numberPC <- 20
# 
# for(i in 1:numberPC){
#   pca_genes_top <-rownames(pca_gene_loading)[order(pca_gene_loading[,i],decreasing=TRUE)]
#   print(paste0("Genes with top positive contribution to PC",i))
#   print(pca_genes_top[1:top])
#   
#   pca_genes_down <-rownames(pca_gene_loading)[order(pca_gene_loading[,i])]
#   print(paste0("Genes with top negative contribution to PC",i))
#   print(pca_genes_down[1:top])
# }
# 
# 
# significant.PCs <- 5
# combinations <- as.data.frame(combn(significant.PCs,2,simplify=TRUE))
# extendXlimLeft <- 5 # Extend left Xlim by x percent of min
# extendXlimRight <- 15 # Extend right Xlim by x percent of max
# extendYlimBottom <- 5 # Extend bottom Ylim by x percent of min
# extendYlimTop <- 5 # Extend top Ylim by x percent of max
# TextSize <- 0.4
# 
# 
# for (i in c(1:4)) {
#   #for(j in 1:ncol(anocol))	{
#   for(j in c(1:5)){
#     png(file.path(resdir,paste0(colnames(anocol)[j],"_",sprintf("PCA%sVs%s_plot_2D.png", combinations[1,i], combinations[2,i] ))), height=1500,width=1500,res=300) 
#     xlimits <- c((min(pca_object[,combinations[1,i]])-abs(min(pca_object[,combinations[1,i]])/100*extendXlimLeft)), (max(pca_object[,combinations[1,i]])+abs(max(pca_object[,combinations[1,i]])/100*extendXlimRight)))
#     ylimits <- c((min(pca_object[,combinations[2,i]])-abs(min(pca_object[,combinations[2,i]])/100*extendYlimBottom)), (max(pca_object[,combinations[2,i]])+abs(max(pca_object[,combinations[2,i]])/100*extendYlimTop)))
#     plot(pca_object[,combinations[1,i]],pca_object[,combinations[2,i]],col=alpha(anocol[,j],0.7),xlab=paste(colnames(pca_object)[combinations[1,i]],100*variances[combinations[1,i]]),ylab=paste(colnames(pca_object)[combinations[2,i]],100*variances[combinations[2,i]]),cex=0.8,lwd=1.5, main=colnames(anocol)[j], pch=19, 
#          xlim=xlimits,
#          ylim=ylimits
#     )
#     dev.off()# changed type of point and xlim
#   }
# }


##############################################################################################
###############			ENRICHMENT ANALYSIS on PC component					##############################
##############################################################################################
organism <- c("mm10","hg19","hg38")[3]
annotationDatabases <- c("MSigDB") # "KEGG" or "GO" or "KEGG_GO" or "MSigDB"
EnrichDir <- file.path(resdir,"EnrichmentPCs") ; if(!file.exists(EnrichDir)){dir.create(EnrichDir)}


if (organism == "mm10") {
  MSigDBFile <- "~/Documents/bioinfo/GeCo.Annotation/msigdb_v5.0_GMTs/MSIG_v5_mousified.RData"
  KEGGFile <- "~/Documents/bioinfo/GeCo.Annotation/mm10/gs.KEGG.toSymbol.mm.RData"
  KEGGdefFile <- "~/Documents/bioinfo/GeCo.Annotation/mm10/gs.KEGG.def.RData" # object KEGGdef for annotating path_id with path_name
  GOFile <- "~/Documents/bioinfo/GeCo.Annotation/mm10/gs.GO.toSymbol.mm_2012-10-03.RData" }

if (organism == "hg19") {
  MSigDBFile <- "~/Documents/bioinfo/GeCo.Annotation/msigdb_v5.0_GMTs/MSIG_v5.RData" }
if (organism == "hg38") {
  MSigDBFile <- "~/Google Drive/GeCo.Bioinf_from_MASTER/GeCo.Annotation/msigdb_v5.0_GMTs/MSIG_v5.RData" }


# Load annotation databases files
if (organism == "hg19" | organism == "hg38") { load(MSigDBFile)}
if (organism == "mm10") {load(MSigDBFile); load(KEGGFile) ;load(KEGGdefFile) ; load(GOFile); KEGGdef <- KEGGdef[,1:2]; colnames(KEGGdef) <- c("path_id","path_name")}


annotbase <- "MSigDB"
Results <- data.frame()
database <- MSIG.ls
reflist <- gene_metadata$Symbol

for(i in 1:20) 	{
  
  pca_genes_top <-rownames(pca_gene_loading)[order(pca_gene_loading[,i],decreasing=TRUE)]
  
  
  enrich.test <- geco.enrichmentTest(gene.sets=database,mylist=pca_genes_top[1:top],possibleIds=reflist)
  enrich.test <- data.frame(Gene_set_name=rownames(enrich.test), enrich.test, check.names=FALSE)
  enrich.test <- merge( subset(MSIG.gs, select=-Genes), enrich.test, by.x="Gene.Set", by.y="Gene_set_name", all.y=TRUE, sort=FALSE ) ## Get class of gene set
  
  enrich.test <- enrich.test[order(enrich.test$`p-value`),]
  #ind <- which(enrich.test$`q-value`<= 0.1);if(!length(ind)){ind <- 1:20}
  #Test <- data.frame(module=paste0("module",i),enrich.test[ind,] )
  Test <- data.frame(module=paste0("PC",i,"_positively"),enrich.test )
  Results <- rbind(Results,Test)
  
  
  pca_genes_down <-rownames(pca_gene_loading)[order(pca_gene_loading[,i])]
  
  
  enrich.test <- geco.enrichmentTest(gene.sets=database,mylist=pca_genes_down[1:top],possibleIds=reflist)
  enrich.test <- data.frame(Gene_set_name=rownames(enrich.test), enrich.test, check.names=FALSE)
  enrich.test <- merge( subset(MSIG.gs, select=-Genes), enrich.test, by.x="Gene.Set", by.y="Gene_set_name", all.y=TRUE, sort=FALSE ) ## Get class of gene set
  
  enrich.test <- enrich.test[order(enrich.test$`p-value`),]
  #ind <- which(enrich.test$`q-value`<= 0.1);if(!length(ind)){ind <- 1:20}
  #Test <- data.frame(module=paste0("module",i),enrich.test[ind,] )
  Test <- data.frame(module=paste0("PC",i,"_negatively"),enrich.test )
  Results <- rbind(Results,Test)
  
} # end of module loop

WriteXLS(Results,ExcelFileName = file.path(EnrichDir,paste0("Enrichment_test_PC_1_20.xlsx")), perl = "perl", verbose = FALSE, row.names = FALSE, col.names = TRUE, AdjWidth = TRUE, AutoFilter = TRUE, BoldHeaderRow = TRUE, na = "", FreezeRow = 1, FreezeCol = 1)


###########FIN VENDREDI 24/04 ##################################################
################################################################################
################ Markers based on louvain clustering --  ##########################
################################################################################

##Top markers (n=25 genes studied in detail per partition)
marker_test_res <- top_markers(cds, group_cells_by="clustCol", reference_cells=1000, cores=8)
top_specific_markers <-  marker_test_res %>% dplyr::filter(fraction_expressing >= 0.30) %>% group_by(cell_group) %>% top_n(5, specificity)
#top_specific_markers <-  marker_test_res %>% dplyr::filter(fraction_expressing >= 0.80,specificity>=0.3) %>% group_by(cell_group)%>% dplyr::slice(which.max(marker_score))
top_specific_marker_ids <-  unique(top_specific_markers %>% pull(gene_id))

plot_genes_by_group(cds,
                    top_specific_marker_ids,
                    group_cells_by="partition",
                    ordering_type="maximal_on_diag",
                    max.size=3)

plot_genes_by_group(cds,
                    top_specific_marker_ids,
                    group_cells_by="partition",
                    ordering_type="cluster_row_col",
                    max.size=3)

##Subset 
cds_subset <- choose_cells(cds)
cds_subset = cluster_cells(cds_subset, resolution=1e-3)
plot_cells(cds_subset, color_cells_by="cluster")
pr_graph_test_res <- graph_test(cds_subset, neighbor_graph="knn", cores=8)
pr_deg_ids <- row.names(subset(pr_graph_test_res, morans_I > 0.01 & q_value < 0.01))
gene_module_df <- find_gene_modules(cds_subset[pr_deg_ids,], resolution=1e-3)

cell_group_df <- tibble::tibble(cell=row.names(colData(cds_subset)), 
                                cell_group=clusters(cds_subset)[colnames(cds_subset)])
agg_mat <- aggregate_gene_expression(cds_subset, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
colnames(agg_mat) <- stringr::str_c("Cluster ", colnames(agg_mat))

pheatmap::pheatmap(agg_mat, cluster_rows=TRUE, cluster_cols=TRUE,
                   scale="column", clustering_method="ward.D2",
                   fontsize=6)

##Detailed marker analysis for all genes in the datasets
pr_graph_test_res <- graph_test(cds, neighbor_graph="knn", cores=8)
pr_deg_ids <- row.names(subset(pr_graph_test_res, morans_I > 0.1 & q_value < 0.01))
#pr_deg_ids <- row.names(subset(pr_graph_test_res, q_value < 0.01))


gene_module_df <- find_gene_modules(cds[pr_deg_ids,], resolution=1e-2)
cell_group_df <- tibble::tibble(cell=row.names(colData(cds)), 
                                cell_group=annot_sel$clustCol)
agg_mat <- aggregate_gene_expression(cds, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
colnames(agg_mat) <- stringr::str_c("Cluster ", colnames(agg_mat))

save(list = c(pr_graph_test_res,pr_deg_ids,gene_module_df,cell_group_df  ))
test <- pheatmap::pheatmap(agg_mat, cluster_rows=TRUE, cluster_cols=TRUE,
                   scale="column", clustering_method="ward.D2",
                   fontsize=6)

plot_cells(cds, genes=gene_module_df, 
           show_trajectory_graph=FALSE, 
           label_cell_groups=FALSE)

##############################################################################################
###############			ENRICHMENT ANALYSIS on gene modules					##############################
##############################################################################################
organism <- c("mm10","hg19","hg38")[3]
annotationDatabases <- c("MSigDB") # "KEGG" or "GO" or "KEGG_GO" or "MSigDB"
EnrichDir <- file.path(resdir,"EnrichmentModules") ; if(!file.exists(EnrichDir)){dir.create(EnrichDir)}


if (organism == "mm10") {
  MSigDBFile <- "~/Documents/bioinfo/GeCo.Annotation/msigdb_v5.0_GMTs/MSIG_v5_mousified.RData"
  KEGGFile <- "~/Documents/bioinfo/GeCo.Annotation/mm10/gs.KEGG.toSymbol.mm.RData"
  KEGGdefFile <- "~/Documents/bioinfo/GeCo.Annotation/mm10/gs.KEGG.def.RData" # object KEGGdef for annotating path_id with path_name
  GOFile <- "~/Documents/bioinfo/GeCo.Annotation/mm10/gs.GO.toSymbol.mm_2012-10-03.RData" }

if (organism == "hg19") {
  MSigDBFile <- "~/Documents/bioinfo/GeCo.Annotation/msigdb_v5.0_GMTs/MSIG_v5.RData" }
if (organism == "hg38") {
  MSigDBFile <- "~/Documents/bioinfo/GeCo.Annotation/msigdb_v5.0_GMTs/MSIG_v5.RData" }


# Load annotation databases files
if (organism == "hg19" | organism == "hg38") { load(MSigDBFile)}
if (organism == "mm10") {load(MSigDBFile); load(KEGGFile) ;load(KEGGdefFile) ; load(GOFile); KEGGdef <- KEGGdef[,1:2]; colnames(KEGGdef) <- c("path_id","path_name")}


  annotbase <- "MSigDB"
  Results <- data.frame()
  database <- MSIG.ls
  reflist <- gene_metadata$Symbol
  
  for(i in 1:dim(table(gene_module_df$module))) 	{
    genes_int <- gene_module_df$id[gene_module_df$module==i]
    

      enrich.test <- geco.enrichmentTest(gene.sets=database,mylist=genes_int,possibleIds=reflist)
      enrich.test <- data.frame(Gene_set_name=rownames(enrich.test), enrich.test, check.names=FALSE)
      enrich.test <- merge( subset(MSIG.gs, select=-Genes), enrich.test, by.x="Gene.Set", by.y="Gene_set_name", all.y=TRUE, sort=FALSE ) ## Get class of gene set
    
      enrich.test <- enrich.test[order(enrich.test$`p-value`),]
      #ind <- which(enrich.test$`q-value`<= 0.1);if(!length(ind)){ind <- 1:20}
      #Test <- data.frame(module=paste0("module",i),enrich.test[ind,] )
      Test <- data.frame(module=paste0("module",i),enrich.test )
      Results <- rbind(Results,Test)
    
   
    
  } # end of module loop
  WriteXLS(Results,ExcelFileName = file.path(EnrichDir,paste0("Enrichment_test_Modules_I0.1_Clusterhierarchique.xlsx")), perl = "perl", verbose = FALSE, row.names = FALSE, col.names = TRUE, AdjWidth = TRUE, AutoFilter = TRUE, BoldHeaderRow = TRUE, na = "", FreezeRow = 1, FreezeCol = 1)
  
  test <- as.matrix(table(annot_sel$louvain_partition,annot_sel$sample_id))
  
  WriteXLS(test,ExcelFileName = file.path(resdir,paste0("CellCounts.xlsx")), perl = "perl", verbose = FALSE, row.names = FALSE, col.names = TRUE, AdjWidth = TRUE, AutoFilter = TRUE, BoldHeaderRow = TRUE, na = "", FreezeRow = 1, FreezeCol = 1)


##### Integrative heatmap of pathway enrichment###############
  
#pathway selection:
  
significPathway <-   unique(Results$Gene.Set[Results$q.value<0.001])
  
  significPathway_breast <- significPathway[grep(pattern = "BREAST",significPathway)]
  
  pathwayInterest <- c("SMID_BREAST_CANCER_BASAL_UP","LIM_MAMMARY_STEM_CELL_UP","CHARAFE_BREAST_CANCER_LUMINAL_VS_MESENCHYMAL_DN","HALLMARK_TNFA_SIGNALING_VIA_NFKB","HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION","CHARAFE_BREAST_CANCER_LUMINAL_VS_BASAL_DN","HALLMARK_MTORC1_SIGNALING","HALLMARK_INTERFERON_GAMMA_RESPONSE","HALLMARK_HYPOXIA")
  significPathway_mammary <- significPathway[grep(pattern = "MAMMARY",significPathway)]
  significPathway_hallmark <- significPathway[grep(pattern = "HALLMARK",significPathway)]
  #order module: 
  hc_row <-test$tree_row
  
  significPathway_interest <- significPathway[grep(pattern = "SMID_BREAST",significPathway)]
  significPathway_interest <- significPathway[grep(pattern = "SMID_BREAST",significPathway)]
  
#Figure breast:
  annot <- Results[Results$Gene.Set %in% pathwayInterest,colnames(Results) %in% c("module","q.value","Gene.Set")]
  annot$logq <- -log10(annot$q.value)
  annotCol <- c("q.value","logq")
  anocol <- geco.annotToCol3(annotS = annot[,annotCol], annotT = annot, maxnumcateg = 6)
  
  annot$value <- anocol[,2]
  annot. <- annot[,c(1,2,5)]   %>% spread(Gene.Set,value)
  rownames(annot.) <- annot.$module
  annot. <- annot.[hc_row$order,-c(1)]
  
  geco.imageCol(as.matrix(annot.))

   
    par(mar=c(2,5,3,2))
  #ordered by drivers
  annot. <- annot[indsamp,] %>% arrange(desc(ERG),desc(ETV1),AR,PTEN,TP53,CDK12,BRCA2,FOXA1,KMT2C,KMT2D,ATM,JAK1,CHD1,COL5A3,NCOR1,SPEN)
  anocol. <- anocol[annot.$Sample,]
  geco.imageCol(anocol.)
  #ordered by RNA cluster
  annot. <- annot[indsamp,] %>% arrange(RNA_cluster,desc(ERG),desc(ETV1),AR,PTEN,TP53,CDK12,BRCA2,FOXA1,KMT2C,KMT2D,ATM,JAK1,CHD1,COL5A3,NCOR1,SPEN)
  anocol. <- anocol[annot.$Sample,]
  geco.imageCol(anocol.)
  #ordered by NE_biopsy
  annot. <- annot[indsamp,] %>% arrange(desc(NE_biopsy),desc(ERG),desc(ETV1),AR,PTEN,TP53,CDK12,BRCA2,FOXA1,KMT2C,KMT2D,ATM,JAK1,CHD1,COL5A3,NCOR1,SPEN)
  anocol. <- anocol[annot.$Sample,]
  geco.imageCol(anocol.)
  #ordered by ARI_pretreated
  annot. <- annot[indsamp,] %>% arrange(desc(ARI_pretreated),desc(ERG),desc(ETV1),AR,PTEN,TP53,CDK12,BRCA2,FOXA1,KMT2C,KMT2D,ATM,JAK1,CHD1,COL5A3,NCOR1,SPEN)
  anocol. <- anocol[annot.$Sample,]
  geco.imageCol(anocol.)
  #ordered by M1_dia
  annot. <- annot[indsamp,] %>% arrange(desc(M1_dia),desc(ERG),desc(ETV1),AR,PTEN,TP53,CDK12,BRCA2,FOXA1,KMT2C,KMT2D,ATM,JAK1,CHD1,COL5A3,NCOR1,SPEN)
  anocol. <- anocol[annot.$Sample,]
  geco.imageCol(anocol.)
  dev.off()
  
  
  