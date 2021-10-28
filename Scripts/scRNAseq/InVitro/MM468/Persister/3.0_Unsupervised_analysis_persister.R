
library(here)
# Directories -------------------------------------------------------------
maindir= here()
resdir <- file.path(maindir,"output","scRNAseq","MM468","Persister")
QCdir <- file.path(maindir,"output","scRNAseq","QC")

resdir <- file.path(resdir, "Unsupervised") ; if(!file.exists(resdir)){dir.create(resdir)}
resdir_UMAP <- file.path(resdir,"UMAP") ; if(!file.exists(resdir_UMAP)){dir.create(resdir_UMAP)}
resdir_boxplots <- file.path(resdir,"boxplots") ; if(!file.exists(resdir_boxplots)){dir.create(resdir_boxplots)}
resdir_heatmaps = file.path(resdir,"Heatmaps"); if(!dir.exists(resdir_heatmaps)) dir.create(resdir_heatmaps)
resdir_sub = file.path(resdir,"sub500"); if(!dir.exists(resdir_sub)) dir.create(resdir_sub)
RDatadir <- file.path(resdir,"RData") ; if(!file.exists(RDatadir)){dir.create(RDatadir)}

source(file.path(maindir,"Scripts","global_var.R"))

annotCol <- c("sample_id","total_features","rRNA","louvain_partition","cons_BC_lenti","CDH2","TWIST1","TGFB1")
# Select initial population, 4 'persister' states (early) and 3 'resistant' states (late)
sample_persisters_study = c(
  "MM468_chemonaive","MM468_5FU1_day33","MM468_5FU3_day50","MM468_5FU2_day67",
  "MM468_5FU3_day77", "MM468_5FU2_day171","MM468_5FU3_day202",
  "MM468_5FU1_day214"
)

#Common palette for all MM468 experiments
control <- c("#E0E0E0","#BDBDBD","#757575")
persister_color <- c("#009688","#DCEDC8","#4CAF50","#9CCC65")
res_color <- c("#FF9800","#FFC107","#FFEB3B","#FF5722")
color_MM468 <- c(control[2],persister_color,res_color[c(1:3)])

corres_sample <- data.frame(
  Sample=sample_persisters_study,
  color=color_MM468)

# If re-computing from scratch - set RECOMPUTE to TRUE else if you just want to
# plot the UMAPS, set RECOMPUTE to FALSE 
RECOMPUTE = FALSE

if(RECOMPUTE){
  
  ################################################################################################
  # PCA, UMAP and louvain clustering with Monocle 3  on all cells for the paper------------------
  ################################################################################################
  load(file.path(RDatadir,"MM468.RData"))
  annot_int <- annot

  annot_int <- annot_int[annot_int$sample_id %in% sample_persisters_study & annot_int$doublet_score<10000,]
  annot_int$sample_id <- as.character(annot_int$sample_id)
  rownames(annot_int) <- annot_int$cell_id
  
  cds <- new_cell_data_set(Signal[row.names(gene_metadata),row.names(annot_int)],
                           cell_metadata = annot_int,
                           gene_metadata  = gene_metadata)
  
  gc()
  cds <- preprocess_cds(cds,method='PCA', norm_method='size_only',num_dim=50)
  
  cds <- reduce_dimension(cds, reduction_method = 'UMAP')
  cds <-  cluster_cells(cds)
  
  umap_res <- cds@int_colData$reducedDims[[2]] # with R3.6.2 and higher for SummarizedExperiment, otherwise cds@reducedDims[[2]]
  pca_object <- cds@int_colData$reducedDims[[1]]
  
  annot_int$louvain_partition <- cds@clusters[[1]]$partitions
  annot_int$louvain_partition <- paste0("C",annot_int$louvain_partition)
  
  #visual check before save
  plot_cells(cds, color_cells_by="sample_id", group_cells_by="partition",label_cell_groups = F)
  #plot_cells(cds,genes="TGFB1")
  
  save(umap_res,file=file.path(RDatadir,"umap_persister.RData"))
  save(pca_object,file=file.path(RDatadir,"pca_persister.RData"))
  save(cds,file=file.path(RDatadir,"cds_persister.RData"))
  
  NormCounts <- t(t(exprs(cds)) /  pData(cds)[, 'Size_Factor'])
  LogCounts <- log(NormCounts+1,2)
  rm(NormCounts)
  gc()
  save(LogCounts,file=file.path(RDatadir,"LogCounts.RData"))
  RawCounts = exprs(cds)
  save(RawCounts,file=file.path(RDatadir,"RawCounts.RData"))
  
  # ################################################################################
  # # Cell cycle scoring with Seurat  ##############################################
  # ################################################################################
  
  # # A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
  # # segregate this list into markers of G2/M phase and markers of S phase
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  
  # Create our Seurat object and complete the initalization steps
  rownames(Signal) <- gene_metadata$Symbol
  marrow <- CreateSeuratObject(counts = Signal[,row.names(annot_int)],meta.data=annot_int)
  marrow <- NormalizeData(marrow)
  marrow <- FindVariableFeatures(marrow, selection.method = "vst")
  marrow <- ScaleData(marrow, features = rownames(marrow))
  marrow <- RunPCA(marrow, features = VariableFeatures(marrow), ndims.print = 1:10, nfeatures.print = 10)
  # DimPlot(marrow, reduction = "pca")
  # DimHeatmap(marrow, dims = 1:3, cells = 1000, balanced = TRUE)
  # ElbowPlot(marrow)
  marrow <- FindNeighbors(marrow, dims = 1:20)
  marrow  <- FindClusters(marrow, resolution = 0.5)
  marrow <- RunUMAP(marrow, dims = 1:20)
  
  #DimPlot(marrow, reduction = "umap")
  
  marrow <- CellCycleScoring(marrow, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  #RidgePlot(marrow, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)
  
  annot_seurat <- data.frame (cell_id=as.character(names(marrow$Phase)),cell_cycle=as.character((marrow$Phase)))
  
  annot_int$cell_cycle <- marrow$Phase[match(annot_int$cell_id,annot_seurat$cell_id)]
  save(annot_int,gene_metadata,file=file.path(RDatadir,"persister_gene_cell_annot.RData"))
  
  ################################################################################################
  ###### Coloring and annotations  #################################
  ################################################################################################
  load(file.path(RDatadir,"LogCounts.RData"))
  load(file.path(RDatadir,"persister_gene_cell_annot.RData"))
  load(file.path(RDatadir,"umap_persister.RData"))
  load(file.path(RDatadir,"pca_persister.RData"))

  rownames(annot_int) <- annot_int$cell_id
  annot_int$sample_id <- as.character(annot_int$sample_id)
  
  annotCol <- c("sample_id","total_features","rRNA","louvain_partition","cons_BC_lenti","CDH2","TWIST1","TGFB1")
  
  annotText <- "sample_id"
  hcText <- "sample_id"  ## column used in hierarchical clustering # change names from Sample_x to actual sample name
  
  ribo <- grepl("^RP[SL]", rownames(LogCounts))
  annot_int$rRNA <- apply(LogCounts[ribo,],2,mean)
  gc()
  set.seed(7)
  
  annot_int$CDH2 <- LogCounts[which(rownames(LogCounts) == "CDH2"), ]
  annot_int$TWIST1 <- LogCounts[which(rownames(LogCounts) == "TWIST1"),]
  annot_int$TGFB1 <- LogCounts[which(rownames(LogCounts) == "TGFB1"),]
  annot_int$LAMB1 <- LogCounts[which(rownames(LogCounts) == "LAMB1"),]
  
  anocol <- geco.unsupervised::geco.annotToCol4(
    annotS=annot_int[,annotCol],annotT=annot_int,plotLegend=T,
    plotLegendFile=file.path(resdir,"Annotation_legends.pdf"), scale_q = "inferno")
  rownames(anocol) = annot_int$cell_id

  anocol[,"sample_id"] <- as.character(corres_sample$color[match(annot_int$sample_id,corres_sample$Sample)])
  
  corres_cell_cycle <- data.frame(phase=c("G1","G2M","S"),color=c("#e896aaff","#553fc2ff","#5d6c7dff"))
  indice <- which(annotCol=="cell_cycle")
  anocol[,indice] <- as.character(corres_cell_cycle$color[match(annot_int$cell_cycle,corres_cell_cycle$phase)])
  
  corres_cluster <- data.frame(cluster=c("C10","C2","C3","C4","C6","C8"),
                               color=c("#b9cced","#438a5e","#ffeadb","#ade498","#ff9c71","#ff847c"))
  corres_cluster$cluster <- as.character(corres_cluster$cluster)
  corres_cluster$color <- as.character(corres_cluster$color)
  indice <- which(annotCol=="louvain_partition")
  anocol[,indice] <- as.character(corres_cluster$color[match(annot_int$louvain_partition,corres_cluster$cluster)])
  
  png(file.path(resdir,"Legend_sample_scRNAseq.png"), height=2000,width=1500,res=300)
  barplot(rep(1,8),col=as.character(corres_sample$color), cex.names  =0.8,horiz=F,names.arg=as.character(corres_sample$Sample),las=2)
  dev.off()
  
  png(file.path(resdir,"Legend_cluster_scRNAseq.png"), height=2000,width=1500,res=300)
  barplot(rep(1,5),col=corres_cluster$color, cex.names  =0.8,horiz=F,
          names.arg=c("C1","C2","C3","C4","C5"),las=2)
  dev.off()
  
  png(file.path(resdir_boxplots,"Boxplot_rRNA.png"),width=1500,height=1500,res=300)
  boxplot(annot_int$rRNA~annot_int$sample_id,las=2)
  dev.off()
  
  save(anocol,file=file.path(RDatadir,"anocol_tokeep.RData"))
  
}

####################################################################################
###### Plot UMAP results ###########################################################
####################################################################################
load(file.path(RDatadir,"anocol.RData"))
load(file.path(RDatadir,"umap_persister.RData"))
load(file.path(RDatadir,"persister_gene_cell_annot.RData"))

umap_res <- umap_res[annot_int$cell_id,]

for(j in 1:ncol(anocol))
{
  png(file.path(resdir_UMAP,paste0("UMAP_",colnames(anocol)[j],".png")), height=1350,width=1200,res=300)
  if(class(annot_int[1,colnames(anocol)[j]])=="numeric"){
  plot((umap_res[annot_int$cell_id,]), col=alpha(anocol[annot_int$cell_id,j],0.2),pch=20,
       cex=0.4,main=paste0(colnames(anocol)[j],
                           " min=",round(min(annot_int[,colnames(anocol)[j]]),digits=3),
                           " max=",round(max(annot_int[,colnames(anocol)[j]]),digits=3)),
       xlab="component 1",ylab="component 2")
  } else{
    plot((umap_res[annot_int$cell_id,]), col=alpha(anocol[,j],0.2),pch=20,
         cex=0.4,main=paste0(colnames(anocol)[j]),
         xlab="component 1",ylab="component 2")
  }
  dev.off()
}

indice <- which(annotCol=="cons_BC_lenti")
if(length(indice) >0 ){
  png(file.path(resdir_UMAP,paste0("UMAP_",colnames(anocol)[indice],".png")), height=1150,width=950,res=350)
  plot((umap_res[!is.na(annot_int$cons_BC_lenti),]),
       col=alpha(anocol[!is.na(annot_int$cons_BC_lenti),indice],0.5),
       pch=20,cex=0.4,main=paste0(colnames(anocol)[indice]," perplexity=30"),xlab="component 1",ylab="component 2",ylim=c(min(umap_res[,2]),max(umap_res[,2])))
  dev.off()
  
  png(file.path(resdir_UMAP,paste0("exp3_UMAP_",colnames(anocol)[indice],".png")), height=1150,width=950,res=350)
  indice <- which(annotCol=="cons_BC_lenti")
  plot((umap_res[!is.na(annot_int$cons_BC_lenti) & annot_int$sample_id %in% c("MM468_chemonaive","MM468_5FU3_day50","MM468_5FU3_day77","MM468_5FU3_day202"),]), col=alpha(anocol[!is.na(annot_int$cons_BC_lenti) & annot_int$sample_id %in% c("MM468_chemonaive","MM468_5FU3_day50","MM468_5FU3_day77","MM468_5FU3_day202"),indice],0.8),pch=20,cex=0.4,main=paste0(colnames(anocol)[indice]," perplexity=30"),
       xlab="component 1",ylab="component 2",ylim=c(min(umap_res[,2]),max(umap_res[,2])), xlim=c(min(umap_res[,1]),max(umap_res[,1])))
  dev.off()
  
  png(file.path(resdir_UMAP,paste0("exp5_UMAP_",colnames(anocol)[indice],".png")), height=1150,width=950,res=350)
  indice <- which(annotCol=="cons_BC_lenti")
  plot((umap_res[!is.na(annot_int$cons_BC_lenti) & annot_int$sample_id %in% c("MM468_chemonaive","MM468_5FU5_day67","MM468_5FU5_day171"),]), col=alpha(anocol[!is.na(annot_int$cons_BC_lenti) & annot_int$sample_id %in% c("MM468_chemonaive","MM468_5FU5_day67","MM468_5FU5_day171"),indice],0.8),pch=20,cex=0.4,main=paste0(colnames(anocol)[indice]," perplexity=30"),
       xlab="component 1",ylab="component 2",ylim=c(min(umap_res[,2]),max(umap_res[,2])), xlim=c(min(umap_res[,1]),max(umap_res[,1])))
  dev.off()
  
  png(file.path(resdir_UMAP,paste0("exp6_UMAP_",colnames(anocol)[indice],".png")), height=1150,width=950,res=350)
  indice <- which(annotCol=="cons_BC_lenti")
  plot((umap_res[!is.na(annot_int$cons_BC_lenti) & annot_int$sample_id %in% c("MM468_chemonaive","MM468_5FU6_day33","MM468_5FU6_day214"),]), col=alpha(anocol[!is.na(annot_int$cons_BC_lenti) & annot_int$sample_id %in% c("MM468_chemonaive","MM468_5FU6_day33","MM468_5FU6_day214"),indice],0.8),pch=20,cex=0.4,main=paste0(colnames(anocol)[indice]," perplexity=30"),
       xlab="component 1",ylab="component 2",ylim=c(min(umap_res[,2]),max(umap_res[,2])), xlim=c(min(umap_res[,1]),max(umap_res[,1])))
  dev.off()
  
  #count number of barcodes detected per expression group: C10, C2, C3, C4, C6, C8
  stats_detection <- annot_int %>% mutate(
    sample_id=factor(sample_id,levels=c(sample_persisters_study))) %>%
    group_by(sample_id) %>% 
    summarise(
      cells=length(cell_id),barcode=length(which(!is.na(cons_BC_lenti))),
      fraction=length(which(!is.na(cons_BC_lenti)))/length(cell_id)
    )
  
  stats_barcode <- annot_int[!is.na(annot_int$cons_BC_lenti),] %>%
    mutate(sample_id=factor(sample_id,levels=sample_persisters_study)) %>%
    group_by(sample_id,louvain_partition) %>%
    summarise(diff_barcode=length(unique(cons_BC_lenti)),
              barcode=length(which(!is.na(cons_BC_lenti))),
              diversity=length(unique(cons_BC_lenti))/length(which(!is.na(cons_BC_lenti)))) %>% 
    filter(barcode>10) 
  
  
  png(file.path(resdir_boxplots,"Lineage_barcode_detection_scRNA.png"), height=1500,width=1500,res=300)
  ggplot(stats_detection, aes(x=sample_id,y=fraction,fill=sample_id,label=barcode)) + 
    geom_bar(position='dodge',stat='identity') +theme_classic() + 
    scale_fill_manual(values=as.character(corres_sample$color)) + 
    geom_text(position = position_dodge(width = 0.8), angle = 90 )
  dev.off()
  
  png(file.path(resdir_boxplots,"Lineage_diversity.png"), height=800,width=1500,res=300)
  ggplot(stats_barcode, aes(x=louvain_partition,y=diversity, fill=sample_id)) +
    geom_bar(stat='identity',position = position_dodge2(width = 0.9, preserve = "single")) + 
    theme_classic() +scale_fill_manual(values=as.character(corres_sample$color))
  dev.off()
}

WriteXLS(as.data.frame(stats_barcode),ExcelFileName=file.path(resdir,"NumberUniqueLineage_sample_cluster.xls"))
contingency <- as.data.frame.matrix(table(annot_int$sample_id,annot_int$louvain_partition))
WriteXLS((contingency),ExcelFileName=file.path(resdir,"NumberCell_cluster_sampleid.xls"),row.names = T)

################################################################################################################
# Hierarchical clustering to group samples in complement to Louvain clustering on a subset of cells ############
################################################################################################################
anocol <- as.data.frame(anocol)
rownames(anocol) <- annot_int$cell_id

set.seed(47)
annot_subset <- annot_int %>% group_by(sample_id) %>% sample_n(500)
annot_subset <- as.data.frame(annot_subset)
annot_subset$cell_id <- as.character(annot_subset$cell_id)
anocol_subset <- anocol[annot_subset$cell_id,]

mati <- t(as.data.frame(pca_object)[annot_subset$cell_id,])
hc_PCA <- hclust(as.dist(1 - cor(mati)), method = "ward.D")
mat.so.cor <- mati[,hc_PCA$order]

png(file.path(resdir_sub,"Clustering_correlation_matrix_PCA_subset.png"), height=1500,width=1500,res=300)
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

subset_LogCounts <- LogCounts[,annot_subset$cell_id]
save(subset_LogCounts,hc_PCA,anocol_subset,annot_subset,gene_metadata,file=file.path(RDatadir,"subset_analysis.RData"))

##############################################################################################
############################### Consensus Hierarchical Clustering  ###########################
##############################################################################################
load(file.path(RDatadir,"subset_analysis.RData"))
maxKCC <- 6
repsCC <- 50
consclust <- ConsensusClusterPlus(
  mati, maxK=maxKCC,plot="png",reps=repsCC, pItem=pItemCC,
  pFeature=pFeatureCC, title=file.path(resdir_heatmaps,"Consensus_clustering"),
  clusterAlg="hc", distance=distCC, innerLinkage=innerLinkageCC, finalLinkage=finalLinkageCC,seed=3.14)
save(consclust,file=file.path(RDatadir,"Consensus_clustering.RData"))
#---------------------------------------------------------------------------------

#Plot repartition of Expression groups with time
corres_day <- data.frame(Sample=sample_persisters_study,
                         day=gsub(".*day","",sample_persisters_study))
corres_day[which(corres_day$Sample=="MM468_chemonaive"),"day"] = 0
corres_day$day = as.numeric(corres_day$day)
annot_subset$day <- as.integer(corres_day$day[match(annot_subset$sample_id,corres_sample$Sample)])

annot_subset$col <- anocol_subset[,"louvain_partition"]
table(annot_subset$louvain_partition,annot_subset$col)
test <- annot_subset %>% mutate(
  louvain_partition=factor(louvain_partition, levels = corres_cluster$cluster)) %>%
  mutate(sample_id=factor(
    sample_id,levels=sample_persisters_study)) %>%
  group_by(sample_id,louvain_partition) %>% summarise(freq=n()) 
test = test[-which(is.na(test$louvain_partition)),]

png(file.path(resdir_boxplots,"Louvain_partition_repartition.png"),width=1500,height=1500,res=300)
ggplot(test) + geom_col(aes(y=freq,x=sample_id,fill=louvain_partition)) + 
  scale_fill_manual(values=corres_cluster$color) + 
  theme_classic() + theme(axis.text.x = element_text(angle = 90)) 
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
methHC <- c("ward","ward.D","ward.D2","single","complete","average")[1]
if(distHC=="distPearson")	d <- distPearson(tmat)
if(distHC=="distCosine")	d <- distCosine(tmat)
if(distHC %in% c("euclidean","maximum","manhattan","canberra","binary","minkowski"))	d <- dist(tmat,distHC)
hc <- hclust(d,method=methHC)

mat.so <- as.matrix(mat[,hc$order])

if(chRangeHM){
  for(i in 1:nrow(mat.so)){
    mat.so[i,] <- geco.changeRange(mat.so[i,],newmin=0,newmax=1)
  }
}
rowClust <- hclust(as.dist(1 - cor(t(mat.so))), method = "ward.D")

png(file.path(resdir_heatmaps,paste0("Clustering_MostVarying_",NbGenes,"_",distHC,"_",methHC,".png")), height=4000,width=4000,res=300)
geco.hclustAnnotHeatmapPlot(x=(mat.so[rowClust$order,]),
                            hc=hc,
                            hmColors=hmColors,
                            anocol=as.matrix(anocol_subset[hc$order,]),#[,ncol(cc.col):1]
                            xpos=c(0.15,0.9,0.164,0.885),
                            ypos=c(0.1,0.5,0.5,0.6,0.62,0.95),
                            dendro.cex=0.01,
                            xlab.cex=0.5,
                            hmRowNames=TRUE,
                            hmRowNames.cex=0.2
)
dev.off()

#reorder clustering
nclust <- 6
annot_subset$clustCol_reorder <- paste("C",match(cutree(hc,nclust),unique(cutree(hc,nclust)[hc$order])),sep="")
test <- c(which(annot_subset$clustCol_reorder %in% c("C4")), which(annot_subset$clustCol_reorder %in% "C6"),
          which(annot_subset$clustCol_reorder %in% "C5"),which(annot_subset$clustCol_reorder %in% "C2"),
          which(annot_subset$clustCol_reorder %in% "C3"),which(annot_subset$clustCol_reorder %in% "C1"))

hc_reorder <- dendextend::rotate(hc,test)

mat.so <- as.matrix(mat[,hc_reorder$order])

if(chRangeHM){
  for(i in 1:nrow(mat.so)){
    mat.so[i,] <- geco.changeRange(mat.so[i,],newmin=0,newmax=1)
  }
}
rowClust <- hclust(as.dist(1 - cor(t(mat.so))), method = "ward.D")

png(file.path(resdir_heatmaps,paste0("Reorder_Clustering_MostVarying_",NbGenes,"_",distHC,"_",methHC,".png")), height=4000,width=4000,res=300)

geco.hclustAnnotHeatmapPlot(x=(mat.so[rowClust$order,]),
                            hc=hc_reorder,
                            hmColors=hmColors,
                            anocol=as.matrix(anocol_subset[hc_reorder$order,]),
                            xpos=c(0.15,0.9,0.164,0.885),
                            ypos=c(0.1,0.5,0.5,0.6,0.62,0.95),
                            dendro.cex=0.01,
                            xlab.cex=0.5,
                            hmRowNames=TRUE,
                            hmRowNames.cex=0.2
)
dev.off()



################################################################################
#Plot distribution of correlation scores and vignettes
################################################################################
interm <- data.frame(cor(mati))
interm$first_cell <- rownames(interm)
all_corr <- tidyr::pivot_longer(interm,cols=1:4000,names_to="other_cell",values_to="Corr")
all_corr <- dplyr::distinct(all_corr)
all_corr$first_cell_ExpressionGroup <- annot_subset$louvain_partition[match(all_corr$first_cell,annot_subset$cell_id)]
all_corr$other_cell_ExpressionGroup <- annot_subset$louvain_partition[match(all_corr$other_cell,annot_subset$cell_id)]

#Plot intra-group corr_scores
intra_corr <- all_corr[all_corr$first_cell_ExpressionGroup==all_corr$other_cell_ExpressionGroup,]
intra_corr <- as.data.frame(intra_corr)

mat.so.cor <- mati[,hc_PCA$order]

taille_groupe <- as.numeric(table(annot_subset$ExpressionGroup))

IntraCorr_test <- data.frame(expressionGroup=corres_cluster$cluster,pvalue=0)

#ref sample for wilcox testing, C10 cluster
MAT <- intra_corr[intra_corr$first_cell_ExpressionGroup=="C10",]
MAT <- tidyr::pivot_wider(MAT[,c(1:3)],names_from="other_cell",values_from="Corr")
MAT_ref <- as.matrix(MAT[,-1])

#to keep same color scale, viridis(100) for -1 to 1 correlation scores.
for(i in 1:6){
  
  Group <- as.character(corres_cluster$cluster[i])
  Col_clust <- as.character(corres_cluster$color[i])
  MAT <- intra_corr[intra_corr$first_cell_ExpressionGroup==Group,]
  MAT <- tidyr::pivot_wider(MAT[,c(1:3)],names_from="other_cell",values_from="Corr")
  test <- wilcox.test(as.matrix(MAT[,-1]),MAT_ref)
  IntraCorr_test$pvalue[i] <- test$p.value
  
  png(file.path(resdir_heatmaps,paste0("Pearson_Correlation_Inferno_scores_",Group,".png")), height=1150,width=1000,res=300)
  MAT <- mat.so.cor[,colnames(mat.so.cor) %in% annot_subset$cell_id[annot_subset$louvain_partition==Group]]
  debut <- round(100-(1-min(cor(MAT)))/0.02,0)
  
  image(cor(MAT),col=(inferno(100)[debut:100]))
  dev.off()
}

WriteXLS(IntraCorr_test,file.path(resdir_heatmaps,"IntraCorr_test.xls"),row.names = T)

intra_corr <- intra_corr %>% 
  mutate(first_cell_ExpressionGroup=factor(first_cell_ExpressionGroup,levels=c("C10","C2","C4","C3","C6","C8")))

png(file.path(resdir_boxplots,"ExpressionGroup_intracorrelation.png"),width=3000,height=1500,res=300)
ggplot(arrange(intra_corr,first_cell_ExpressionGroup),aes(x=Corr, fill=first_cell_ExpressionGroup)) +
  geom_density(alpha=0.6) + theme_classic(base_size = 25) + 
  scale_fill_manual(values=as.character(corres_cluster$color))
dev.off()


sel <- (intra_corr$other_cell_ExpressionGroup %in% corres_cluster$cluster)

png(file.path(resdir_boxplots,paste0("ExpressionGroup_intracorrelation_violin.png")),height=1000,width=1500,res=300)
ggplot(intra_corr[sel,],aes(x=first_cell_ExpressionGroup,y=Corr, fill=first_cell_ExpressionGroup)) + 
  geom_violin(alpha=0.8) + theme_classic() + scale_fill_manual(values=as.character(corres_cluster$color)) +
  stat_summary(fun.y=median, geom="point", size=2, color="black")
dev.off()


##Lineage diversity
mat_index <- annot_int %>% group_by(sample_id,louvain_partition) %>% summarise(index=length(unique(lineage_barcode))/length(which(!is.na(lineage_barcode))),number=length(which(!is.na(lineage_barcode)))) 
mat_index <- mat_index[mat_index$number>10,]
WriteXLS(as.data.frame(mat_index),ExcelFileName=file.path(resdir,"NumberUniqueLineage_sample_cluster.xls"))
contingency <- as.data.frame.matrix(table(annot_int$sample_id,annot_int$louvain_partition))
WriteXLS((contingency),ExcelFileName=file.path(resdir,"NumberCell_cluster_sampleid.xls"),row.names = T)

#####################################################################################
#  Differences between persisting chemonaive & non persisting chemonaive ############
#####################################################################################

annot_barcoded = annot_int[which(!is.na(annot_int$cons_BC_lenti)),]

annot_barcoded$cell_id = paste0(annot_barcoded$sample_id,"_", annot_barcoded$BC_10x)

dim(annot_barcoded)

umap_barcoded = umap_res[which(!is.na(annot_int$cons_BC_lenti)),]
rownames(umap_barcoded) = annot_barcoded$cell_id


chemonaive_2 = new.env()
load(file.path(RDatadir,"RData_perSample","MM468_chemonaive_2.RData"), envir = chemonaive_2)
annot_chemonaive_2 = chemonaive_2$metadata.MM468_chemonaive_2
annot_chemonaive_2 = annot_chemonaive_2[which(!is.na(annot_chemonaive_2$cons_BC_lenti)),]
annot_chemonaive_2$sample_id = "MM468_chemonaive"
dim(annot_chemonaive_2)

annot_barcoded = rbind(annot_barcoded[,c(1:6,9)],annot_chemonaive_2[,1:7])

barcodes_chemonaive = (annot_barcoded$cons_BC_lenti[which(annot_barcoded$sample_id %in% c("MM468_chemonaive"))])
cells_chemonaive = unique(annot_barcoded$cell_id[which(annot_barcoded$sample_id %in% c("MM468_chemonaive"))])
barcodes_persister = (annot_barcoded$cons_BC_lenti[which(!annot_barcoded$sample_id %in% c("MM468_chemonaive"))])
cells_persister = unique(annot_barcoded$cell_id[which(!annot_barcoded$sample_id %in% c("MM468_chemonaive"))])

persisting_chemonaive = intersect(barcodes_chemonaive, barcodes_persister)
cells_persisting_chemonaive = annot_barcoded$cell_id[which(annot_barcoded$cons_BC_lenti %in% persisting_chemonaive &
                                                             annot_barcoded$sample_id == "MM468_chemonaive")]
cells_non_persisting_chemonaive = annot_barcoded$cell_id[which(!(annot_barcoded$cons_BC_lenti %in% persisting_chemonaive) &
                                                                 annot_barcoded$sample_id == "MM468_chemonaive")]

non_chemonaive_persisters = setdiff(barcodes_persister,barcodes_chemonaive)



png(file.path(resdir_UMAP,paste0("UMAP_persisting_chemonaives.png")), height=1500,width=1400,res=350) 
plot(umap_barcoded[which(rownames(umap_barcoded) %in% cells_chemonaive),],
     col=ifelse(barcodes_chemonaive %in% persisting_chemonaive, "orange", "grey50"),
     pch=20,cex=0.4,main="Chemonaive persisting barcodes",xlab="component 1",ylab="component 2",ylim=c(min(umap_res[,2]),max(umap_res[,2])),
     xlim = c(min(umap_barcoded[match(c(cells_chemonaive, cells_persister),rownames(umap_barcoded)),1], na.rm = T),
              max(umap_barcoded[match(c(cells_chemonaive, cells_persister),rownames(umap_barcoded)),1], na.rm = T)))
dev.off()

space = umap_barcoded[which(rownames(umap_barcoded) %in% cells_chemonaive),]
cells_small_cluster =names(which(space[,1] < -6.6 & space[,2] > -2.5))
cells_big_cluster =names(which(space[,1] < -6.6 & space[,2] < -2.5))
length(cells_small_cluster)
length(intersect(cells_small_cluster, cells_non_persisting_chemonaive) )
length(intersect(cells_small_cluster, cells_persisting_chemonaive) )

length(cells_big_cluster)
length(intersect(cells_big_cluster, cells_non_persisting_chemonaive) )
length(intersect(cells_big_cluster, cells_persisting_chemonaive) )


png(file.path(resdir_UMAP,paste0("UMAP_persisting_chemonaives_and_persister.png")), height=1500,width=1400,res=350) 
plot(umap_barcoded[match(c(cells_chemonaive, cells_persister),rownames(umap_barcoded)),],
     col=c(ifelse(barcodes_chemonaive %in% persisting_chemonaive, "orange", "grey50"),
           ifelse(barcodes_persister %in% persisting_chemonaive, "red", "grey50")),
     pch=20,cex=0.4,main="Chemonaive & persister common barcodes",xlab="component 1",ylab="component 2",ylim=c(min(umap_res[,2]),max(umap_res[,2])))
dev.off()

### Differential analysis between persisting & non persisting cells :
# COMPARISON SETUP

## Type of analysis
algoType <- c("edgeR") ## Limma or DESeq or a vector of both or edgeR for single-cell
useBlockFactor <- FALSE ## TRUE or FALSE
annotationDatabases <- c("MSigDB") # "KEGG" or "GO" or "KEGG_GO" or "MSigDB"

mygps <- list(
  'p_chemonaive' = cells_persisting_chemonaive
)

myrefs <- list(
  'non_p_chemonaive' = cells_non_persisting_chemonaive
) # Sets reference (can be 1 or more samples)


refs <- names(myrefs) 
groups <- names(mygps) 

##raw counts for edgeR
load(file.path(RDatadir,"RawCounts.RData"))
load(file.path(RDatadir,"gene_cell_annot.RData"))
feature <- gene_metadata


# select genes which are detected in at least 1% of cells
Percent=0.01
threshold = Percent*dim(RawCounts)[2]

BinaryCounts = RawCounts
BinaryCounts[BinaryCounts>0]=1
gc()

N_above_0 = rowSums(BinaryCounts)
gc()
sel = which(N_above_0>threshold)

RawCounts <- RawCounts[sel,]
gc()

# Add chemonaive_2 cells to RawCount
RawCounts_chemonaive_2 = chemonaive_2$counts.MM468_chemonaive_2
common_genes = intersect(rownames(RawCounts), rownames(RawCounts_chemonaive_2))
RawCounts_chemonaive_2 = RawCounts_chemonaive_2[match(common_genes, rownames(RawCounts_chemonaive_2)),]
RawCounts = RawCounts[match(common_genes, rownames(RawCounts)),]
feature = data.frame("Symbol" = common_genes, "gene_short_name" = common_genes)
rownames(feature) = common_genes
RawCounts = cbind(RawCounts, RawCounts_chemonaive_2)
gc()

colnames(RawCounts) = gsub("MM468_chemonaive_2","MM468_chemonaive", colnames(RawCounts))
annot_barcoded$cell_id = gsub("MM468_chemonaive_2","MM468_chemonaive", annot_barcoded$cell_id)

myrefs$non_p_chemonaive = gsub("MM468_chemonaive_2","MM468_chemonaive", myrefs$non_p_chemonaive)
mygps$p_chemonaive = gsub("MM468_chemonaive_2","MM468_chemonaive", mygps$p_chemonaive)

RawCounts = RawCounts[,which(colnames(RawCounts) %in% annot_barcoded$cell_id)]

# EdgeR DA
my.res <- geco.CompareedgeRGLM (dataMat=RawCounts,
                                annot=annot_barcoded[which(annot_barcoded$sample_id == "MM468_chemonaive"),],
                                ref=myrefs,
                                groups=mygps,
                                featureTab=feature
)

log2FC_thresholds = log2(c(2))
Signif_threshold = 0.01
tabdir = file.path(resdir, "Tables"); if(!dir.exists(tabdir)) dir.create(tabdir)

for(log2FC_threshold in log2FC_thresholds){
  under_res = my.res %>% dplyr::filter(log2FC.p_chemonaive < -log2FC_threshold & qval.p_chemonaive < Signif_threshold) %>% 
    dplyr::arrange(qval.p_chemonaive) %>% dplyr::select(Symbol,log2FC.p_chemonaive,qval.p_chemonaive)
  over_res = my.res %>% dplyr::filter(log2FC.p_chemonaive > log2FC_threshold & qval.p_chemonaive < Signif_threshold) %>% 
    dplyr::arrange(qval.p_chemonaive) %>% dplyr::select(Symbol,log2FC.p_chemonaive,qval.p_chemonaive)
  
  WriteXLS(c("over_res","under_res","my.res"),
           ExcelFileName = file.path(
             tabdir,paste0("Differential_analysis_Limma_logFC_",
                           round(log2FC_threshold,2),"_p_chemonaive_vs_non_p_chemonaive.xlsx")),
           SheetNames = c(
             paste0("Over_",round(log2FC_threshold,2),"_",Signif_threshold,"_n",
                    nrow(over_res)),paste0(
                      "Under_-",round(log2FC_threshold,2),"_",Signif_threshold,
                      "_n",nrow(under_res)),"All"),
           perl = "perl", verbose = FALSE, row.names = FALSE,
           col.names = TRUE, AdjWidth = T, AutoFilter = TRUE,
           BoldHeaderRow = TRUE, na = "", FreezeRow = 1, FreezeCol = 1)
  
  summaryTab <- geco.summaryCompareedgeR(restab=my.res,
                                         ref=myrefs,
                                         groups=mygps,
                                         qval.th=Signif_threshold,
                                         fc.th=log2FC_threshold,
                                         plotdir=tabdir
  )
  
  write.table(summaryTab, file.path(
    tabdir,paste0("Number_differentially_expressed_genes_per_group_logFC_",
                  round(log2FC_threshold,2),".csv")),row.names=F,quote=F,sep=";")
} ## end of if edgeR
save(my.res, file=file.path(tabdir,paste("Supervised_res_object_edgeR.Rdata",sep="")))
load(file=file.path(tabdir,paste("Supervised_res_object_edgeR.Rdata",sep="")))
print("summaryedgeRCompare Done")


################################################################################
# Find Differences between resisting persister & non resisting persisters
################################################################################

annot_barcoded = annot_int[which(!is.na(annot_int$cons_BC_lenti)),]
annot_barcoded = annot_barcoded[-which(annot_barcoded$sample_id == "MM468_chemonaive"),]
dim(annot_barcoded)
umap_res = umap_res[match(annot_barcoded$cell_id, rownames(umap_res)),]

persister_samples = c("MM468_5FU1_day33","MM468_5FU2_day67","MM468_5FU3_day50","MM468_5FU3_day77")

barcodes_persister = (annot_barcoded$cons_BC_lenti[which(annot_barcoded$sample_id %in% persister_samples)])
cells_persister = unique(annot_barcoded$cell_id[which(annot_barcoded$sample_id %in% persister_samples)])
barcodes_resistant = (annot_barcoded$cons_BC_lenti[which(!annot_barcoded$sample_id %in% persister_samples)])
cells_resistant = unique(annot_barcoded$cell_id[which(!annot_barcoded$sample_id %in% persister_samples)])

resisting_persister = intersect(barcodes_persister, barcodes_resistant)
cells_resisting_persister = annot_barcoded$cell_id[which(annot_barcoded$cons_BC_lenti %in% resisting_persister &
                                                             annot_barcoded$sample_id %in% persister_samples)]
cells_non_resisting_persister = annot_barcoded$cell_id[which(!(annot_barcoded$cons_BC_lenti %in% resisting_persister) &
                                                                 annot_barcoded$sample_id %in% persister_samples)]

non_persister_resistants = setdiff(barcodes_resistant,barcodes_persister)

png(file.path(resdir_UMAP,paste0("UMAP_resisting_persisters.png")), height=1500,width=1400,res=350) 
plot(umap_res[which(rownames(umap_res) %in% cells_persister),],
     col=ifelse(barcodes_persister %in% resisting_persister, "orange", "grey50"),
     pch=20,cex=0.4,main="resisting persister barcodes",xlab="component 1",ylab="component 2",
     ylim=c(min(umap_res[,2]),max(umap_res[,2])))
dev.off()

png(file.path(resdir_UMAP,paste0("UMAP_resisting_persisters_and_resistant.png")), height=1500,width=1400,res=350) 
plot(umap_res[match(c(cells_persister, cells_resistant),rownames(umap_res)),],
     col=c(ifelse(barcodes_persister %in% resisting_persister, "orange", "grey50"),
           ifelse(barcodes_resistant %in% resisting_persister, "red", "grey50")),
     pch=20,cex=0.4,main="persister & resistant common barcodes",xlab="component 1",ylab="component 2",ylim=c(min(umap_res[,2]),max(umap_res[,2])))
dev.off()

### Differential analysis between persisting & non persisting cells :
# COMPARISON SETUP

log2FC_thresholds <- log2(3)
Signif_threshold <- 0.01

## Type of analysis
algoType <- c("edgeR") ## Limma or DESeq or a vector of both or edgeR for single-cell
useBlockFactor <- FALSE ## TRUE or FALSE
annotationDatabases <- c("MSigDB") # "KEGG" or "GO" or "KEGG_GO" or "MSigDB"

mygps <- list(
  'r_persister' = cells_resisting_persister
)

myrefs <- list(
  'non_r_persister' = cells_non_resisting_persister
) # Sets reference (can be 1 or more samples)


refs <- names(myrefs) 
groups <- names(mygps) 

##raw counts for edgeR
load(file.path(RDatadir,"RawCounts.RData"))
feature <- gene_metadata

sel = which(colnames(RawCounts) %in% paste0(annot_barcoded$sample_id,"_",annot_barcoded$BC_10x))

RawCounts = RawCounts[,sel]
gc()

# select genes which are detected in at least 1% of cells
Percent=0.01
sel <- which(apply(RawCounts,1,function(x) length(which(x>0)))>Percent*dim(RawCounts)[2])
gc()
RawCounts <- RawCounts[sel,];feature <- feature[sel,]
gc()

colnames(RawCounts) =  annot_barcoded$cell_id[match(colnames(RawCounts), paste0(annot_barcoded$sample_id,"_",annot_barcoded$BC_10x))]
# EdgeR DA
my.res <- geco.CompareedgeRGLM (dataMat=RawCounts,
                                annot=annot_barcoded[which(annot_barcoded$sample_id %in% persister_samples),],
                                ref=myrefs,
                                groups=mygps,
                                featureTab=feature
)

log2FC_thresholds = log2(c(2))
tabdir = file.path(resdir, "Tables"); if(!dir.exists(tabdir)) dir.create(tabdir)

for(log2FC_threshold in log2FC_thresholds){
  under_res = my.res %>% dplyr::filter(log2FC.r_persister < -log2FC_threshold & qval.r_persister < Signif_threshold) %>% 
    dplyr::arrange(qval.r_persister) %>% dplyr::select(Symbol,log2FC.r_persister,qval.r_persister)
  over_res = my.res %>% dplyr::filter(log2FC.r_persister > log2FC_threshold & qval.r_persister < Signif_threshold) %>% 
    dplyr::arrange(qval.r_persister) %>% dplyr::select(Symbol,log2FC.r_persister,qval.r_persister)
  
  WriteXLS(c("over_res","under_res","my.res"),
           ExcelFileName = file.path(
             tabdir,paste0("DA_r_persister_vs_non_r_persister_",
                           round(log2FC_threshold,2),".xlsx")),
           SheetNames = c(
             paste0("Over_",round(log2FC_threshold,2),"_",Signif_threshold,"_n",
                    nrow(over_res)),paste0(
                      "Under_-",round(log2FC_threshold,2),"_",Signif_threshold,
                      "_n",nrow(under_res)),"All"),
           perl = "perl", verbose = FALSE, row.names = FALSE,
           col.names = TRUE, AdjWidth = T, AutoFilter = TRUE,
           BoldHeaderRow = TRUE, na = "", FreezeRow = 1, FreezeCol = 1)
  
  summaryTab <- geco.summaryCompareedgeR(restab=my.res,
                                         ref=myrefs,
                                         groups=mygps,
                                         qval.th=Signif_threshold,
                                         fc.th=log2FC_threshold,
                                         plotdir=tabdir
  )
  
  write.table(summaryTab, file.path(
    tabdir,paste0("Number_differentially_expressed_genes_per_group_logFC_",
                  round(log2FC_threshold,2),".csv")),row.names=F,quote=F,sep=";")
} ## end of if edgeR
print("summaryedgeRCompare Done")


annotbase <- "MSigDB" 
database <- MSIG.ls ##MSigDB

Overexpressed  <- Underexpressed <- data.frame()

#my.res$Gene<- sub("hg19_","",my.res$Gene)

#rownames(my.res) <- my.res$Gene
reflist <- unique(my.res$Symbol);length(reflist)

for(i in 1:length(groups)) 	{
  for(log2FC_threshold in log2FC_thresholds){
    gp <- groups[i]
    if (length(refs)>1) {ref <- refs[i]} else {ref <- refs[1]}
    print(paste0("Processing ",gp, " vs ", ref, " _ ",annotbase ))
    
    signific <- which(my.res[,paste("qval",gp,sep=".")] <= Signif_threshold & abs(my.res[,paste("log2FC",gp,sep=".")]) > log2FC_threshold)
    over <- which(my.res[,paste("qval",gp,sep=".")] <= Signif_threshold & my.res[,paste("log2FC",gp,sep=".")] > log2FC_threshold)
    under <- which(my.res[,paste("qval",gp,sep=".")] <= Signif_threshold & my.res[,paste("log2FC",gp,sep=".")] < -log2FC_threshold)
    print(paste0("significant = ", length(signific))) ; print(paste0("over = ", length(over))) ; print(paste0("under = ", length(under)))
    
    
    if(length(over)){
      enrich.test <- geco.enrichmentTest(gene.sets=database,mylist=my.res$Symbol[over],possibleIds=reflist)
      enrich.test <- data.frame(Gene_set_name=rownames(enrich.test), enrich.test, check.names=FALSE)
      enrich.test <- merge( subset(MSIG.gs, select=-Genes), enrich.test, by.x="Gene.Set", by.y="Gene_set_name", all.y=TRUE, sort=FALSE ) ## Get class of gene set
      enrich.test <- enrich.test[order(enrich.test$`p-value`),]
      enrich.test <- enrich.test[order(enrich.test$`p-value`),]
      ind <- which(enrich.test$`q-value`<= 0.1);if(!length(ind)){ind <- 1:20}
      Overexpressed  <- enrich.test[ind,]		}
      
    if(length(under)){
        enrich.test <- geco.enrichmentTest(gene.sets=database,mylist=my.res$Symbol[under],possibleIds=reflist)
        enrich.test <- data.frame(Gene_set_name=rownames(enrich.test), enrich.test, check.names=FALSE)
        enrich.test <- merge( subset(MSIG.gs, select=-Genes), enrich.test, by.x="Gene.Set", by.y="Gene_set_name", all.y=TRUE, sort=FALSE ) ## Get class of gene set
        enrich.test <- enrich.test[order(enrich.test$`p-value`),]
        ind <- which(enrich.test$`q-value`<= 0.1);if(!length(ind)){ind <- 1:20}
        Underexpressed <- enrich.test[ind,]		}
      
      WriteXLS(
        c("Overexpressed", "Underexpressed"),
        ExcelFileName = file.path(tabdir,
                                  paste0("Enrichment_test_",gp,"_vs_",ref,"_",annotbase,
                                         "_logFC",round(log2FC_threshold,2),".xlsx")), 
        SheetNames = c( paste0("Overexp_in_", gp), paste0("Underexp_in_", gp) ),
        perl = "perl", verbose = FALSE, row.names = FALSE, col.names = TRUE,
        AdjWidth = TRUE, AutoFilter = TRUE, BoldHeaderRow = TRUE, na = "",
        FreezeRow = 1, FreezeCol = 1)
      
      save(Overexpressed,Underexpressed, file = file.path(tabdir,paste0("Enrichment_test_",gp,"_vs_",ref,"_",annotbase,"_logFC",round(log2FC_threshold,2),".RData")))
    }
  }