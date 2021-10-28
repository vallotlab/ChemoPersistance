
library(here)
# Directories -------------------------------------------------------------
maindir= here()
resdir <- file.path(maindir,"output","scRNAseq","MM468","UNC_5FU")
QCdir <- file.path(maindir,"output","scRNAseq","QC")

resdir <- file.path(resdir, "Unsupervised") ; if(!file.exists(resdir)){dir.create(resdir)}
resdir_UMAP <- file.path(resdir,"UMAP") ; if(!file.exists(resdir_UMAP)){dir.create(resdir_UMAP)}
resdir_boxplots <- file.path(resdir,"boxplots") ; if(!file.exists(resdir_boxplots)){dir.create(resdir_boxplots)}
resdir_heatmaps = file.path(resdir,"Heatmaps"); if(!dir.exists(resdir_heatmaps)) dir.create(resdir_heatmaps)
resdir_sub = file.path(resdir,"sub500"); if(!dir.exists(resdir_sub)) dir.create(resdir_sub)
resdir_sublineage = file.path(resdir,"sublineage"); if(!dir.exists(resdir_sub)) dir.create(resdir_sub)

RDatadir <- file.path(resdir,"RData") ; if(!file.exists(RDatadir)){dir.create(RDatadir)}

source(file.path(maindir,"Scripts","global_var.R"))
# source("~/Desktop/geco.annotToCol4_Pacôme.R")
# be careful with geco.annotTocol4 version

annotCol <- c("sample_id","total_features","louvain_partition","cons_BC_lenti","CDH2","TWIST1","TGFB1","KRT14")
# Select initial population, 4 'persister' states (early) and 3 'resistant' states (late)
sample_UNC_5FU_study = c(
  "MM468_chemonaive", "MM468_chemonaive_2", "MM468_5FU7_day21","MM468_5FU1_day33","MM468_5FU3_day50","MM468_5FU2_day67", "MM468_5FU7_UNC_day21", "MM468_5FU1_UNC_day33"
)


#Common palette for all MM468 experiments
control <- c("#E0E0E0","#BDBDBD","#757575")
pers_color <- c("#009688","#DCEDC8","#4CAF50","#9CCC65","#009053ff")
res_color <- c("#FF9800","#FFC107","#FFEB3B","#FF5722")
treated_UNC_color <- c("#2196F3", "#00028cff" )
treated_GSKJ4_color <- c("#C2185B", "#EC407A" )

color_MM468 <- c(control[2:1],pers_color[c(5,1:3)],treated_UNC_color[2:1])

corres_sample <- data.frame(
  Sample=sample_UNC_5FU_study,
  color=color_MM468)

barplot(rep(1, 8), cex.names = 0.5, col = color_MM468, names.arg = corres_sample$Sample)
# If re-computing from scratch - set RECOMPUTE to TRUE else if you just want to
# plot the UMAPS, set RECOMPUTE to FALSE 
RECOMPUTE = TRUE

if(RECOMPUTE){
  
  ################################################################################################
  # PCA, UMAP and louvain clustering with Monocle 3  on all cells for the paper------------------
  ################################################################################################
  load(file.path(RDatadir,"MM468_UNC_5FU_All.RData"))
  #load(file.path(RDatadir,"UNC_5FU_gene_cell_annot.RData"))
  annot_int <- annot
  
  annot_int <- annot_int[annot_int$doublet_score<10000,]
  annot_int$sample_id <- as.character(annot_int$sample_id)
  
  rownames(annot_int) <- annot_int$cell_id
  library(monocle3)
  cds <- new_cell_data_set(Signal[row.names(gene_metadata),row.names(annot_int)],
                           cell_metadata = annot_int,
                           gene_metadata  = gene_metadata)
  
  gc()
  cds <- preprocess_cds(cds,method='PCA', norm_method='size_only',num_dim=50)
  
  # cds@colData$batch = ifelse(cds@colData$sample_id %in% c("MM468_5FU7_day21","MM468_5FU7_UNC_day21"),"batch_day21","other_batch")
  # cds@colData$batch = as.factor(cds@colData$batch)
  # table(cds@colData$batch, cds@colData$sample_id)
  # cds = align_cds(cds, num_dim = 50, alignment_group = "batch")

  cds <- reduce_dimension(cds, reduction_method = 'UMAP')
  cds <-  cluster_cells(cds)
  
  markers = monocle3::top_markers(cds, group_cells_by = "sample_id")
  
  plot_cells(cds, genes = "CAMK2N1")
  
  pdf(file.path(resdir_UMAP, "markers_5FU7_UNC_day21.pdf"))
  for(i in markers$gene_id[which(markers$cell_group == "MM468_5FU7_UNC_day21")]){
    print(plot_cells(cds, genes = i,  group_cells_by = "partition"))
  }
  dev.off()
  
  pdf(file.path(resdir_UMAP, "markers_5FU7_day21.pdf"))
  for(i in markers$gene_id[which(markers$cell_group == "MM468_5FU7_day21")]){
    print(plot_cells(cds, genes = i,  group_cells_by = "partition"))
  }
  dev.off()
  
  markers_cluster = monocle3::top_markers(cds, group_cells_by = "cluster")
  
  pdf(file.path(resdir_UMAP, "markers_5FU7_day21_cluster8.pdf"))
  for(i in markers_cluster$gene_id[which(markers_cluster$cell_group == 8)]){
    print(plot_cells(cds, genes = i,  group_cells_by = "cluster"))
  }
  dev.off()
  
  plot_cells(cds, genes = "CAMK2N1")
  
  umap_res <- cds@int_colData@listData$reducedDims@listData$UMAP[,1:2] # with R3.6.2 and higher for SummarizedExperiment, otherwise cds@reducedDims[[2]]
  pca_object <- cds@int_colData$reducedDims[[1]]
  
  annot_int$louvain_partition <- cds@clusters[[1]]$partitions
  annot_int$louvain_partition <- paste0("C",annot_int$louvain_partition)
  
  #visual check before save
  plot_cells(cds, color_cells_by="sample_id", group_cells_by="partition",label_cell_groups = F)
  plot_cells(cds, color_cells_by="partition", group_cells_by="partition",label_cell_groups = F)
  
  
  save(umap_res,file=file.path(RDatadir,"umap_UNC_5FU.RData"))
  save(pca_object,file=file.path(RDatadir,"pca_UNC_5FU.RData"))
  
  RawCounts = exprs(cds)
  size_factor = pData(cds)[, 'Size_Factor']
  rm(cds)
  gc()
  save(RawCounts,file=file.path(RDatadir,"UNC_5FU_RawCounts.RData"))
  
  NormCounts <- t(t(RawCounts) /  size_factor)
  rm(RawCounts); rm(Signal); 
  gc()
  NormCounts <- as.matrix(NormCounts) + 1
  gc()
  UNC_5FU_LogCounts <- log(NormCounts,2)
  rm(NormCounts)
  gc()
  UNC_5FU_LogCounts = as(UNC_5FU_LogCounts,"dgCMatrix")
  save(UNC_5FU_LogCounts,file=file.path(RDatadir,"UNC_5FU_LogCounts.RData"))
  gc()
  
  save(annot_int, gene_metadata,file=file.path(RDatadir,"UNC_5FU_gene_cell_annot_louvain.RData"))
  
  ################################################################################################
  ###### Parameters for à façon plots of results on all cells--  #################################
  ################################################################################################
  #Transient loading of persister annot_int and anocol to retreive reference colors for BC
  load(file.path(maindir, "output","scRNAseq","MM468","Persister","Unsupervised","Rdata","anocol.RData"))
  load(file.path(maindir, "output","scRNAseq","MM468","Persister","Unsupervised","Rdata","persister_gene_cell_annot.RData"))

  #load(file.path(maindir,"output","scRNAseq","MM468_PDX","common_over_genes_pers_vs_unt_log2FC1.58.RData"))
  load(file.path(RDatadir,"UNC_5FU_LogCounts.RData"))
  load(file.path(RDatadir,"umap_UNC_5FU.RData"))
  load(file.path(RDatadir,"pca_UNC_5FU.RData"))
  load(file.path(RDatadir,"UNC_5FU_gene_cell_annot_louvain.RData"))
  
  annot_int$sample_id <- as.character(annot_int$sample_id)
  
  annotText <- "sample_id"
  hcText <- "sample_id"  ## column used in hierarchical clustering # change names from Sample_x to actual sample name
  
  # ribo <- grepl("^RP[SL]",gene_metadata$Symbol)
  # annot_int$rRNA <- apply(UNC_LogCounts[ribo,],2,mean)
  gc()
  set.seed(2047)
  
  annot_int$CDH2 <- UNC_5FU_LogCounts["CDH2",]
  annot_int$TWIST1 <- UNC_5FU_LogCounts["TWIST1",]
  annot_int$TGFB1 <- UNC_5FU_LogCounts["TGFB1",]
  annot_int$KRT14 <- UNC_5FU_LogCounts["KRT14",]
  
  anocol <- geco.annotToCol4(
    annotS=annot_int,annotT=annot_int,plotLegend=T,
    plotLegendFile=file.path(resdir,"Annotation_legends.pdf"), scale_q = "inferno")
  
  
  anocol[,"sample_id"] <- as.character(corres_sample$color[match(annot_int$sample_id,corres_sample$Sample)])
  
  corres_cluster <- data.frame(
    cluster=unique(annot_int$louvain_partition),
    color= palette.colors(8))
  
  anocol[,"louvain_partition"] <- as.character(corres_cluster$color[match(annot_int$louvain_partition,corres_cluster$cluster)])
  
  corres_lenti_BC <- data.frame(
    lenti_BC = annot_int$cons_BC_lenti,
    cell_id = annot_int$cell_id,
    color = anocol[,"cons_BC_lenti"]
  )
  
  corres_lenti_BC = corres_lenti_BC[which(corres_lenti_BC$cell_id %in% annot_int$cell_id),]
  in_Pers = which(rownames(anocol) %in% as.character(corres_lenti_BC$cell_id))
  
  anocol[in_Pers,"cons_BC_lenti"] =  as.character(corres_lenti_BC$color[match(rownames(anocol)[in_Pers],corres_lenti_BC$cell_id)])
  
  
  png(file.path(resdir,"Legend_sample_scRNAseq.png"), height=2000,width=1500,res=300)
  barplot(rep(1,8),col=as.character(corres_sample$color), cex.names  =0.8,horiz=F,names.arg=corres_sample$Sample,las=2)
  dev.off()
  
  png(file.path(resdir,"Legend_cluster_scRNAseq.png"), height=2000,width=1500,res=300)
  barplot(rep(1,8),col=as.character(corres_cluster$color), cex.names  =0.8,horiz=F,names.arg=corres_cluster$cluster,las=2)
  dev.off()
  
  # png(file.path(resdir,"Legend_cluster_scRNAseq.png"), height=2000,width=1500,res=300)
  # barplot(rep(1,5),col=unique(anocol[,"louvain_partition"]), cex.names  =0.8,horiz=F,
  #         names.arg=c("C1","C2","C3","C4","C5"),las=2)
  # dev.off()
  
  png(file.path(resdir_boxplots,"Boxplot_rRNA.png"),width=1500,height=1500,res=300)
  boxplot(annot_int$rRNA~annot_int$sample_id,las=2)
  dev.off()
  
  save(anocol,file=file.path(RDatadir,"anocol.RData"))
  save(annot_int,file=file.path(RDatadir,"annot_int.RData"))
}

####################################################################################
###### Plot UMAP results ###########################################################
####################################################################################
load(file.path(RDatadir,"annot_int.RData"))
load(file.path(RDatadir,"anocol.RData"))
load(file.path(RDatadir,"pca_UNC_5FU.RData"))
load(file.path(RDatadir,"umap_UNC_5FU.RData"))

umap_res <- umap_res[annot_int$cell_id,]

for(j in 1:ncol(anocol))
{
  png(file.path(resdir_UMAP,paste0("UMAP_",colnames(anocol)[j],".png")), height=1500,width=1400,res=350)
  if(class(annot_int[1,colnames(anocol)[j]])=="numeric"){
    plot((umap_res[annot_int$cell_id,]), col=alpha(anocol[,j],0.15),pch=20,
         cex=0.6,main=paste0(colnames(anocol)[j],
                             " min=",round(min(annot_int[,colnames(anocol)[j]]),digits=3),
                             " max=",round(max(annot_int[,colnames(anocol)[j]]),digits=3)),
         xlab="component 1",ylab="component 2")
  } else{
    plot((umap_res[annot_int$cell_id,]), col=alpha(anocol[,j],0.15),pch=20,
         cex=0.6,main=paste0(colnames(anocol)[j]),
         xlab="component 1",ylab="component 2")
  }
  dev.off()
}


indice <- which(colnames(anocol)=="cons_BC_lenti")
if(length(indice) >0 ){
  png(file.path(resdir_UMAP,paste0("UMAP_",colnames(anocol)[indice],".png")), height=1500,width=1400,res=350) 
  plot((umap_res[!is.na(annot_int$cons_BC_lenti),]),
       col=alpha(anocol[!is.na(annot_int$cons_BC_lenti),indice],0.5),
       pch=20,cex=0.4,main=paste0(colnames(anocol)[indice]," perplexity=30"),xlab="component 1",ylab="component 2",ylim=c(min(umap_res[,2]),max(umap_res[,2])))
  dev.off()
  
  png(file.path(resdir_UMAP,paste0("exp3_UMAP_",colnames(anocol)[indice],".png")),  height=1150,width=950,res=350)
  plot((umap_res[!is.na(annot_int$cons_BC_lenti) & annot_int$sample_id %in% c("MM468_chemonaive","MM468_chemonaive_2","MM468_5FU3_day50"),]), col=alpha(anocol[!is.na(annot_int$cons_BC_lenti) & annot_int$sample_id %in% c("MM468_chemonaive","MM468_5FU3_day50","MM468_5FU3_day77","MM468_5FU3_day202"),indice],0.8),pch=20,cex=0.4,main=paste0(colnames(anocol)[indice]," perplexity=30"),
       xlab="component 1",ylab="component 2",ylim=c(min(umap_res[,2]),max(umap_res[,2])), xlim=c(min(umap_res[,1]),max(umap_res[,1])))
  dev.off()
  
  png(file.path(resdir_UMAP,paste0("exp2_UMAP_",colnames(anocol)[indice],".png")),  height=1150,width=950,res=350) 
  plot((umap_res[!is.na(annot_int$cons_BC_lenti) & annot_int$sample_id %in% c("MM468_chemonaive","MM468_chemonaive_2","MM468_5FU2_day67"),]), col=alpha(anocol[!is.na(annot_int$cons_BC_lenti) & annot_int$sample_id %in% c("MM468_chemonaive","MM468_5FU5_day67","MM468_5FU5_day171"),indice],0.8),pch=20,cex=0.4,main=paste0(colnames(anocol)[indice]," perplexity=30"),
       xlab="component 1",ylab="component 2",ylim=c(min(umap_res[,2]),max(umap_res[,2])), xlim=c(min(umap_res[,1]),max(umap_res[,1])))
  dev.off()
  
  png(file.path(resdir_UMAP,paste0("exp1_UMAP_",colnames(anocol)[indice],".png")), height=1150,width=950,res=350) 
  plot((umap_res[!is.na(annot_int$cons_BC_lenti) & annot_int$sample_id %in% c("MM468_chemonaive","MM468_chemonaive_2","MM468_5FU1_day33","MM468_5FU1_UNC_day33"),]), col=alpha(anocol[!is.na(annot_int$cons_BC_lenti) & annot_int$sample_id %in% c("MM468_chemonaive","MM468_5FU6_day33","MM468_5FU6_day214"),indice],0.8),pch=20,cex=0.4,main=paste0(colnames(anocol)[indice]," perplexity=30"),
       xlab="component 1",ylab="component 2",ylim=c(min(umap_res[,2]),max(umap_res[,2])), xlim=c(min(umap_res[,1]),max(umap_res[,1])))
  dev.off()
  
  png(file.path(resdir_UMAP,paste0("exp7_UMAP_",colnames(anocol)[indice],".png")), height=1150,width=950,res=350) 
  plot((umap_res[!is.na(annot_int$cons_BC_lenti) & annot_int$sample_id %in% c("MM468_chemonaive","MM468_chemonaive_2","MM468_5FU7_day21","MM468_5FU7_UNC_day21"),]), col=alpha(anocol[!is.na(annot_int$cons_BC_lenti) & annot_int$sample_id %in% c("MM468_chemonaive","MM468_5FU6_day33","MM468_5FU6_day214"),indice],0.8),pch=20,cex=0.4,main=paste0(colnames(anocol)[indice]," perplexity=30"),
       xlab="component 1",ylab="component 2",ylim=c(min(umap_res[,2]),max(umap_res[,2])), xlim=c(min(umap_res[,1]),max(umap_res[,1])))
  dev.off()
  
  #count number of barcodes detected per expression group: C10, C2, C3, C4, C6, C8
  stats_detection <- annot_int %>% mutate(
    sample_id=factor(sample_id,levels=c(sample_UNC_5FU_study))) %>%
    group_by(sample_id) %>% 
    summarise(
      cells=length(cell_id),barcode=length(which(!is.na(cons_BC_lenti))),
      fraction=length(which(!is.na(cons_BC_lenti)))/length(cell_id)
    )
  
  stats_barcode <- annot_int[!is.na(annot_int$cons_BC_lenti),] %>%
    mutate(sample_id=factor(sample_id,levels=sample_UNC_5FU_study)) %>%
    group_by(sample_id,louvain_partition) %>%
    summarise(diff_barcode=length(unique(cons_BC_lenti)),
              barcode=length(cons_BC_lenti),
              diversity=length(unique(cons_BC_lenti))/length(cons_BC_lenti)) %>% 
    filter(barcode>20) 
  
  
  png(file.path(resdir_boxplots,"Lineage_barcode_detection_scRNA.png"), height=1500,width=1500,res=300)
  ggplot(stats_detection, aes(x=sample_id,y=fraction,fill=sample_id,label=barcode)) + 
    geom_bar(position='dodge',stat='identity') +theme_classic() + 
    scale_fill_manual(values=as.character(corres_sample$color)) + 
    geom_text(position = position_dodge(width = 0.8), angle = 90 ) +
    theme(axis.text.x = element_text(angle = 90))
  dev.off()
  
  png(file.path(resdir_boxplots,"Lineage_diversity.png"), height=800,width=1500,res=300)
  ggplot(stats_barcode, aes(x=louvain_partition,y=diversity, fill=sample_id)) +
    geom_bar(stat='identity',position = position_dodge2(width = 0.9, preserve = "single")) + 
    theme_classic() +scale_fill_manual(values=as.character(corres_sample$color))
  dev.off()
}


#Tables
#Lineage diversity
WriteXLS(stats_barcode,ExcelFileName=file.path(resdir,"NumberUniqueLineage_sample_cluster.xls"))
contingency <- as.data.frame.matrix(table(annot_int$sample_id,annot_int$louvain_partition))
WriteXLS((contingency),ExcelFileName=file.path(resdir,"NumberCell_cluster_sampleid.xls"),row.names = T)



################################################################################################################
# Hierarchical clustering to group samples in complement to Louvain clustering on a subset of cells ############
################################################################################################################
#take a subset, do ONCE
load(file.path(RDatadir,"annot_int.RData"))
load(file.path(RDatadir,"anocol.RData"))
load(file.path(RDatadir,"pca_UNC_5FU.RData"))
load(file.path(RDatadir,"umap_UNC_5FU.RData"))

anocol <- as.data.frame(anocol)
rownames(anocol) <- annot_int$cell_id

set.seed(47)
#annot_subset <- annot_int %>% group_by(sample_id) %>% sample_n(500)
annot_subset <- annot_int[!is.na(annot_int$cons_BC_lenti),]

# annot_int = as.data.frame(annot_subset)
pca_object <- pca_object[annot_subset$cell_id,]
anocol_subset <- anocol[annot_subset$cell_id,]

set.seed(47)
mati <- t(as.data.frame(pca_object))
hc_PCA <- hclust(as.dist(1 - cor(mati)), method = "ward.D")
mat.so.cor <- mati[,hc_PCA$order]

annot_subset$hierarchical_cluster = paste0("C", cutree(hc_PCA,k = 4))
corres_hcluster <- data.frame(
  cluster=unique(annot_subset$hierarchical_cluster),
  color=brewer.pal(4, "Set2"))

png(file.path(resdir_sublineage,"Legend_cluster_scRNAseq.png"), height=2000,width=1500,res=300)
barplot(rep(1,4),col=as.character(corres_hcluster$color), cex.names  =0.8,horiz=F,names.arg=corres_hcluster$cluster,las=2)
dev.off()

anocol_subset$hierarchical_cluster <- corres_hcluster$color[match(annot_subset$hierarchical_cluster,corres_hcluster$cluster)]
#anocol. = annotToCol2(annot_int) 
#anocol[,"hierarchical_cluster"] = anocol.[,"hierarchical_cluster"]

png(file.path(resdir_sublineage,"Clustering_correlation_matrix_PCA_subset.png"), height=1500,width=1500,res=300)
geco.hclustAnnotHeatmapPlot(x=cor(mat.so.cor),
                            hc=hc_PCA,
                            hmColors=corColors,
                            anocol=as.matrix(anocol_subset[hc_PCA$order,c(3,4,7,15,20)]),#[,ncol(cc.col):1]
                            xpos=c(0.15,0.9,0.164,0.885),
                            ypos=c(0.1,0.5,0.5,0.6,0.62,0.95),
                            dendro.cex=0.1,
                            xlab.cex=0.2,
                            hmRowNames=FALSE,
                            hmRowNames.cex=0.01
)
dev.off()


umap_res <- umap_res[annot_subset$cell_id,]

for(j in 1:ncol(anocol_subset))
{
  png(file.path(resdir_sublineage,paste0("UMAP_subset_",colnames(anocol_subset)[j],".png")), height=1500,width=1400,res=350) 
  if(class(annot_int[1,colnames(anocol)[j]])=="numeric"){
    plot((umap_res[annot_subset$cell_id,]), col=alpha(anocol_subset[,j],0.5),pch=20,
         cex=0.6,main=paste0(colnames(anocol)[j],
                             " min=",round(min(annot_subset[,colnames(anocol)[j]]),digits=3),
                             " max=",round(max(annot_subset[,colnames(anocol)[j]]),digits=3)),
         xlab="component 1",ylab="component 2")
  } else{
    plot((umap_res[annot_subset$cell_id,]), col=alpha(anocol_subset[,j],0.5),pch=20,
         cex=0.6,main=paste0(colnames(anocol_subset)[j]),
         xlab="component 1",ylab="component 2")
  }
  dev.off()
}

stats_detection <- annot_subset %>% mutate(
  sample_id=factor(sample_id,levels=c(sample_UNC_5FU_study))) %>%
  group_by(sample_id) %>% 
  summarise(
    cells=length(cell_id),barcode=length(which(!is.na(cons_BC_lenti))),
    fraction=length(which(!is.na(cons_BC_lenti)))/length(cell_id)
  )

stats_barcode_louvain <- annot_subset[!is.na(annot_subset$cons_BC_lenti),] %>%
  mutate(sample_id=factor(sample_id,levels=sample_UNC_5FU_study)) %>%
  group_by(sample_id,louvain_partition) %>%
  summarise(diff_barcode=length(unique(cons_BC_lenti)),
            barcode=length(cons_BC_lenti),
            diversity=length(unique(cons_BC_lenti))/length(cons_BC_lenti)) %>% 
  filter(barcode>20) 

stats_barcode_hc <- annot_subset[!is.na(annot_subset$cons_BC_lenti),] %>%
  mutate(sample_id=factor(sample_id,levels=sample_UNC_5FU_study)) %>%
  group_by(sample_id,hierarchical_cluster) %>%
  summarise(diff_barcode=length(unique(cons_BC_lenti)),
            barcode=length(cons_BC_lenti),
            diversity=length(unique(cons_BC_lenti))/length(cons_BC_lenti)) %>% 
  filter(barcode>20) 


png(file.path(resdir_sub,"Lineage_barcode_detection_scRNA.png"), height=1500,width=1500,res=300)
ggplot(stats_detection, aes(x=sample_id,y=fraction,fill=sample_id,label=barcode)) + 
  geom_bar(position='dodge',stat='identity') +theme_classic() + 
  scale_fill_manual(values=as.character(corres_sample$color)) + 
  geom_text(position = position_dodge(width = 0.8), angle = 90 ) +
  theme(axis.text.x = element_text(angle = 90))
dev.off()

png(file.path(resdir_sub,"Lineage_diversity_louvain.png"), height=800,width=1500,res=300)
ggplot(stats_barcode_louvain, aes(x=louvain_partition,y=diversity, fill=sample_id)) +
  geom_bar(stat='identity',position = position_dodge2(width = 0.9, preserve = "single")) + 
  theme_classic() +scale_fill_manual(values=as.character(corres_sample$color))
dev.off()

png(file.path(resdir_sub,"Lineage_diversity_hierarchical.png"), height=800,width=1500,res=300)
ggplot(stats_barcode_hc, aes(x=hierarchical_cluster,y=diversity, fill=sample_id)) +
  geom_bar(stat='identity',position = position_dodge2(width = 0.9, preserve = "single")) + 
  theme_classic() +scale_fill_manual(values=as.character(corres_sample$color))
dev.off()


WriteXLS(stats_barcode_hc,ExcelFileName=file.path(resdir_sub,"NumberUniqueLineage_sample_cluster_hc.xls"))
contingency <- as.data.frame.matrix(table(annot_int$sample_id,annot_int$louvain_partition))
WriteXLS((contingency),ExcelFileName=file.path(resdir,"NumberCell_cluster_sampleid.xls"),row.names = T)


#Plot repartition of Expression groups with time
corres_day <- data.frame(Sample=sample_UNC_5FU_study,
                         day=gsub(".*day","",sample_UNC_5FU_study))
corres_day[which(corres_day$Sample=="MM468_chemonaive"),"day"] = 0
corres_day[which(corres_day$Sample=="MM468_chemonaive_2"),"day"] = 0
corres_day$day = as.numeric(corres_day$day)
annot_int$day <- as.integer(corres_day$day[match(annot_int$sample_id,corres_sample$Sample)])

annot_int$col <- annot_int[,"louvain_partition"]
table(annot_int$louvain_partition,annot_int$col)
test <- annot_int %>% mutate(
  louvain_partition=factor(louvain_partition, levels = paste0("C",1:9))) %>%
  mutate(sample_id=factor(
    sample_id,levels=sample_UNC_5FU_study)) %>%
  group_by(sample_id,louvain_partition) %>% summarise(freq=n()) 
# test = test[-which(is.na(test$louvain_partition)),]

png(file.path(resdir_boxplots,"Louvain_partition_repartition.png"),width=1500,height=1500,res=300)
ggplot(test) + geom_col(aes(y=freq,x=sample_id,fill=louvain_partition)) + 
  scale_fill_manual(values=unique(anocol[,"louvain_partition"])) + 
  theme_classic() + theme(axis.text.x = element_text(angle = 90)) 
dev.off()

################################################################################
###Generate heatmaps with top varying genes 
################################################################################
#work on subset_LogCounts 
NbGenes <- 100
exp.sd <- apply(UNC_5FU_LogCounts,1,sd)
sel <- order(exp.sd,decreasing=T)[1:NbGenes]
mat <- UNC_5FU_LogCounts[sel,];dim(mat)
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

annot_int$hierarchical_cluster = paste0("C", cutree(hc,k = 7))
anocol. = annotToCol2(annot_int) 
anocol[,"hierarchical_cluster"] = anocol.[,"hierarchical_cluster"]

png(file.path(resdir_heatmaps,paste0("Clustering_MostVarying_",NbGenes,"_",distHC,"_",methHC,".png")), height=4000,width=4000,res=300)
geco.hclustAnnotHeatmapPlot(x=(mat.so[rowClust$order,]),
                            hc=hc,
                            hmColors=hmColors,
                            anocol=as.matrix(anocol[hc$order,c(1,2,4,6)]),#[,ncol(cc.col):1]
                            xpos=c(0.15,0.9,0.164,0.885),
                            ypos=c(0.1,0.5,0.5,0.6,0.62,0.95),
                            dendro.cex=0.01,
                            xlab.cex=0.5,
                            hmRowNames=TRUE,
                            hmRowNames.cex=0.5
)
dev.off()


##########################
# Plot UMAPs #############
load(file=file.path(RDatadir,"anocol.RData"))
load(file=file.path(RDatadir,"UNC_LogCounts.RData"))
load(file=file.path(RDatadir,"umap_UNC.RData"))
load(file.path(RDataSupdir,"Supervised_res_object_edgeR.Rdata"))

diff_genes <-  c("KRT17", "INHBB", "KRT5","KRT8","KRT14","BMP6","CDKN2B","TGFBR2",
                                        "TGFBR3","TGFB2","INHBA","INHBB","SMAD2","TGFB1",
                                        "FOSL1","NNMT","KRT14","KLK10","KLK5","TAGLN","NNMT",
                                        "ELN","TGFBR3","MIF", "ABCC4","ABCC3","CD24","VIM","CDH2",
                                        "CDH1","ABCA5","LBH")

pcaText <- FALSE
annotText <- "sample_id"



annotCol = unique(c(colnames(annot_int),diff_genes))


for(i in annotCol){
  if(i %in% rownames(UNC_5FU_LogCounts))
    annot_int[,i] <- UNC_5FU_LogCounts[which(rownames(UNC_5FU_LogCounts)==i),]
}
gc()

anocol2 <- geco.unsupervised::geco.annotToCol4(annotS=annot_int[,annotCol],plotLegend = F, scale_q = "inferno")
save(anocol2,annot_int,file=file.path(RDatadir,"annot_anocol_final.RData")) 

for(i in setdiff(colnames(anocol2),colnames(anocol)))
{
  j = which(colnames(anocol2) == i)
  png(file.path(resdir_UMAP,paste0("UMAP_",colnames(anocol2)[j],".png")), height=1500,width=1500,res=300)
  if(class(annot_int[1,colnames(anocol2)[j]])=="numeric"){
    plot((umap_res), col=alpha(anocol2[,j],0.3),pch=20,cex=0.6,
         main=paste0(colnames(anocol2)[j]," min=",round(min(annot_int[,colnames(anocol2)[j]]),digits=3)," max=",
                     round(max(annot_int[,colnames(anocol2)[j]]),digits=3)),
         xlab="component 1",ylab="component 2")} else {
           plot((umap_res), col=alpha(anocol2[,j],0.3),pch=20,cex=0.6,
                main=paste0(colnames(anocol2)[j]),
                xlab="component 1",ylab="component 2")
           if(colnames(anocol2)[j]=="sample_id"){
             #Plot UNC with add plot since there is too low cell number
             pers = which(annot_int$sample_id=="HBCx95_persister_6souris")
             points(umap_res[pers,], col=alpha(anocol2[pers,j],0.5),pch=20,cex=0.6)
             
           }
         }
  
  
  dev.off()
}

#Lineage diversity
mat_index <- annot_int %>% group_by(sample_id,louvain_partition) %>% summarise(index=length(unique(cons_BC_lenti))/length(which(!is.na(cons_BC_lenti))),number=length(which(!is.na(cons_BC_lenti)))) 
mat_index <- mat_index[mat_index$number>10,]
WriteXLS(as.data.frame(mat_index),ExcelFileName=file.path(resdir,"NumberUniqueLineage_sample_cluster.xls"))
contingency <- as.data.frame.matrix(table(annot_int$sample_id,annot_int$louvain_partition))
WriteXLS((contingency),ExcelFileName=file.path(resdir,"NumberCell_cluster_sampleid.xls"),row.names = T)

  