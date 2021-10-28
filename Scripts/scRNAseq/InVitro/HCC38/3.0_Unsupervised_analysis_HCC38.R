library(here)
# Directories -------------------------------------------------------------
maindir= here()
resdir <- file.path(maindir,"output","scRNAseq","HCC38","Persister")
QCdir <- file.path(maindir,"output","scRNAseq","QC")

resdir <- file.path(resdir, "Unsupervised") ; if(!file.exists(resdir)){dir.create(resdir)}
resdir_UMAP <- file.path(resdir,"UMAP") ; if(!file.exists(resdir_UMAP)){dir.create(resdir_UMAP)}
resdir_boxplots <- file.path(resdir,"boxplots") ; if(!file.exists(resdir_boxplots)){dir.create(resdir_boxplots)}
resdir_heatmaps = file.path(resdir,"Heatmaps"); if(!dir.exists(resdir_heatmaps)) dir.create(resdir_heatmaps)
resdir_sub = file.path(resdir,"sub500"); if(!dir.exists(resdir_sub)) dir.create(resdir_sub)
RDatadir <- file.path(resdir,"RData") ; if(!file.exists(RDatadir)){dir.create(RDatadir)}

source(file.path(maindir,"Scripts","global_var.R"))

annotCol <- c("sample_id","total_features","rRNA","louvain_partition","CDH2","TWIST1","TGFB1")
# Select initial population, 4 'persister' states (early) and 3 'resistant' states (late)
sample_persisters_study = c(
  "HCC38_chemonaive","HCC38_persister"
)

#Common palette for all HCC38 experiments
control <- c("#E0E0E0","#BDBDBD","#757575")
persister_color <- c("#009688","#DCEDC8","#4CAF50","#9CCC65")
res_color <- c("#FF9800","#FFC107","#FFEB3B","#FF5722")
color_HCC38 <- c(control[2],persister_color[1])

corres_sample <- data.frame(
  Sample=sample_persisters_study,
  color=color_HCC38)

# If re-computing from scratch - set RECOMPUTE to TRUE else if you just want to
# plot the UMAPS, set RECOMPUTE to FALSE 
RECOMPUTE = TRUE

if(RECOMPUTE){
  load(file.path(RDatadir,"HCC38.RData"))
  annot_int <- annot
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
  
  
  ################################################################################################
  # PCA, UMAP and louvain clustering with Monocle 3  on all cells for the paper------------------
  ##############################################################################################
  library(monocle3)
  annot_int <- annot_int[annot_int$sample_id %in% sample_persisters_study & annot_int$doublet_score<10000,]
  annot_int$sample_id <- as.character(annot_int$sample_id)
  rownames(annot_int) <- annot_int$cell_id
  
  cds <- new_cell_data_set(Signal[row.names(gene_metadata),row.names(annot_int)],
                           cell_metadata = annot_int,
                           gene_metadata  = gene_metadata)
  
  gc()
  cds <- preprocess_cds(cds,method='PCA', norm_method='size_only',num_dim=50)
  
  cds <- reduce_dimension(cds, reduction_method = 'UMAP')
  cds <-  cluster_cells(cds, resolution = 0.5, k = 250)
  
  pdf(file.path(resdir_UMAP,"monocle_UMAPs.pdf"))
  monocle3::plot_cells(cds, color_cells_by = "sample_id",reduction_method = "UMAP", 
                       group_label_size = 4, cell_stroke =2 )
  monocle3::plot_cells(cds, color_cells_by = "total_counts",reduction_method = "UMAP", 
                       group_label_size = 4, cell_stroke =2 )
  monocle3::plot_cells(cds, color_cells_by = "total_features",reduction_method = "UMAP", 
                       group_label_size = 4, cell_stroke =2 )
  monocle3::plot_cells(cds, color_cells_by = "doublet_score",reduction_method = "UMAP", 
                       group_label_size = 4, cell_stroke =2 )
  monocle3::plot_cells(cds, color_cells_by = "cluster",reduction_method = "UMAP", 
                       group_label_size = 4, cell_stroke =2 )
  dev.off()
  
  umap_res <- cds@int_colData$reducedDims[[2]] # with R3.6.2 and higher for SummarizedExperiment, otherwise cds@reducedDims[[2]]
  pca_object <- cds@int_colData$reducedDims[[1]]
  
  annot_int$louvain_partition <- cds@clusters[[1]]$partitions
  annot_int$louvain_partition <- paste0("C",annot_int$louvain_partition)
  annot_int$louvain_cluster <- cds@clusters[[1]]$clusters
  annot_int$louvain_cluster <- paste0("C",annot_int$louvain_cluster)
  
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
  
  annotCol <- c("sample_id","total_features","rRNA","louvain_cluster","CDH2","TWIST1","TGFB1")
  
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
  
  anocol <- geco.annotToCol4(
    annotS=annot_int[,annotCol],annotT=annot_int,plotLegend=T,
    plotLegendFile=file.path(resdir,"Annotation_legends.pdf"), scale_q = "inferno")
  rownames(anocol) = annot_int$cell_id

  anocol[,"sample_id"] <- as.character(corres_sample$color[match(annot_int$sample_id,corres_sample$Sample)])
  
  corres_cell_cycle <- data.frame(phase=c("G1","G2M","S"),color=c("#e896aaff","#553fc2ff","#5d6c7dff"))
  indice <- which(annotCol=="cell_cycle")
  anocol[,indice] <- as.character(corres_cell_cycle$color[match(annot_int$cell_cycle,corres_cell_cycle$phase)])

  png(file.path(resdir,"Legend_sample_scRNAseq.png"), height=2000,width=1500,res=300)
  barplot(rep(1,2),col=as.character(corres_sample$color), cex.names  =0.8,horiz=F,names.arg=as.character(corres_sample$Sample),las=2)
  dev.off()
  
  png(file.path(resdir,"Legend_cluster_scRNAseq.png"), height=2000,width=1500,res=300)
  barplot(rep(1,2),col=corres_cluster$color, cex.names  =0.8,horiz=F,
          names.arg=c("C1","C2"),las=2)
  dev.off()
  
  png(file.path(resdir_boxplots,"Boxplot_rRNA.png"),width=1500,height=1500,res=300)
  boxplot(annot_int$rRNA~annot_int$sample_id,las=2)
  dev.off()
  
  save(anocol,file=file.path(RDatadir,"anocol.RData"))
  
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

################################################################################################################
# Hierarchical clustering to group samples in complement to Louvain clustering on a subset of cells ############
################################################################################################################
anocol <- as.data.frame(anocol)
rownames(anocol) <- annot_int$cell_id

set.seed(47)
annot_subset <- annot_int
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
