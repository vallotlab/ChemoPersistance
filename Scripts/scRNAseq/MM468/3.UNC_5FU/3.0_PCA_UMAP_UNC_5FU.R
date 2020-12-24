
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
RDatadir <- file.path(resdir,"RData") ; if(!file.exists(RDatadir)){dir.create(RDatadir)}

source(file.path(maindir,"Scripts","global_var.R"))

################################################################################################
# PCA, UMAP and louvain clustering with Monocle 3  on all cells for the paper------------------
################################################################################################
load(file.path(maindir,"output","scRNAseq","MM468","Persister","Unsupervised","RData","MM468.RData"))
load(file.path(RDatadir,"UNC_5FU_gene_cell_annot.RData"))
annot_int <- annot

table(annot_int$sample_id)
# Select initial population, 4 'persister' states (early) and 3 'resistant' states (late)
sample_UNC_5FU_study = c(
  "MM468_chemonaive","MM468_5FU6_day33","MM468_5FU3_day50","MM468_5FU5_day67",
  "MM468_5FU3_day77", "MM468_5FU6_UNC_day33"
  )

annot_int <- annot_int[annot_int$sample_id %in% sample_UNC_5FU_study & annot_int$doublet_score<10000,]
annot_int$sample_id <- as.character(annot_int$sample_id)

# set.seed(47)
# annot_subset <- annot_int %>% group_by(sample_id) %>% sample_n(1000)
# annot_int = as.data.frame(annot_subset)

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

save(umap_res,file=file.path(RDatadir,"umap_UNC.RData"))
save(pca_object,file=file.path(RDatadir,"pca_UNC.RData"))

RawCounts = exprs(cds)
size_factor = pData(cds)[, 'Size_Factor']
rm(cds)
gc()
save(RawCounts,file=file.path(RDatadir,"UNC_RawCounts.RData"))

NormCounts <- t(t(RawCounts) /  size_factor)
rm(RawCounts); rm(Signal); gc()
NormCounts <- as.matrix(NormCounts) + 1
gc()
UNC_LogCounts <- log(NormCounts,2)
rm(NormCounts)
gc()
UNC_LogCounts = as(UNC_LogCounts,"dgCMatrix")
save(UNC_LogCounts,file=file.path(RDatadir,"UNC_LogCounts.RData"))
gc()

################################################################################################
###### Parameters for à façon plots of results on all cells--  #################################
################################################################################################

load(file.path(maindir,"output","scRNAseq","MM468_PDX","common_over_genes_pers_vs_unt_log2FC1.58.RData"))
load(file.path(RDatadir,"UNC_LogCounts.RData"))
load(file.path(RDatadir,"umap_UNC.RData"))
load(file.path(RDatadir,"pca_UNC.RData"))
load(file.path(RDatadir,"annot_int.RData"))

annot_int$sample_id <- as.character(annot_int$sample_id)
annotCol <- c("sample_id","total_features","rRNA","louvain_partition","cons_BC_lenti")

annotText <- "sample_id"
hcText <- "sample_id"  ## column used in hierarchical clustering # change names from Sample_x to actual sample name

ribo <- grepl("^RP[SL]",gene_metadata$Symbol)
annot_int$rRNA <- apply(UNC_LogCounts[ribo,],2,mean)
gc()
set.seed(2047)
anocol <- geco.unsupervised::geco.annotToCol4(
  annotS=annot_int[,annotCol],annotT=annot_int,plotLegend=T,
  plotLegendFile=file.path(resdir,"Annotation_legends.pdf"), scale_q = "inferno")

#Common palette for all MM468 experiments
control <- c("#E0E0E0","#BDBDBD","#757575")
pers_color <- c("#009688","#DCEDC8","#4CAF50","#9CCC65")
res_color <- c("#FF9800","#FFC107","#FFEB3B","#FF5722")
treated_UNC_color <- c("#2196F3", "#4DD0E1" )
treated_GSKJ4_color <- c("#C2185B", "#EC407A" )

color_MM468 <- c(control[2],pers_color[1:4],treated_UNC_color[1])

corres_sample <- data.frame(
  Sample=sample_UNC_5FU_study,
  color=color_MM468)

anocol[,"sample_id"] <- as.character(corres_sample$color[match(annot_int$sample_id,corres_sample$Sample)])

Persister = new.env()
load(file.path(resdir,"..","..","Persister","Unsupervised","RData","anocol.RData"),Persister)
load(file.path(resdir,"..","..","Persister","Unsupervised","RData","gene_cell_annot_persister.RData"),Persister)
corres_lenti_BC <- data.frame(
  lenti_BC = Persister$annot_int$cons_BC_lenti,
  cell_id = Persister$annot_int$cell_id,
  color = Persister$anocol[,"cons_BC_lenti"]
)
corres_lenti_BC = corres_lenti_BC[which(corres_lenti_BC$cell_id %in% annot_int$cell_id),]
in_Pers = which(rownames(anocol) %in% corres_lenti_BC$cell_id)

anocol[in_Pers,"cons_BC_lenti"] =  corres_lenti_BC$color[match(rownames(anocol)[in_Pers],corres_lenti_BC$cell_id)]


png(file.path(resdir,"Legend_sample_scRNAseq.png"), height=2000,width=1500,res=300)
barplot(rep(1,6),col=corres_sample$color, cex.names  =0.8,horiz=F,names.arg=corres_sample$Sample,las=2)
dev.off()

png(file.path(resdir,"Legend_cluster_scRNAseq.png"), height=2000,width=1500,res=300)
barplot(rep(1,5),col=unique(anocol[,"louvain_partition"]), cex.names  =0.8,horiz=F,
        names.arg=c("C1","C2","C3","C4","C5"),las=2)
dev.off()

png(file.path(resdir_boxplots,"Boxplot_rRNA.png"),width=1500,height=1500,res=300)
boxplot(annot_int$rRNA~annot_int$sample_id,las=2)
dev.off()

save(anocol,file=file.path(RDatadir,"anocol.RData"))
save(annot_int,file=file.path(RDatadir,"annot_int.RData"))

####################################################################################
###### Plot UMAP results ###########################################################
####################################################################################
load(file.path(RDatadir,"annot_int.RData"))
load(file.path(RDatadir,"anocol.RData"))

umap_res <- umap_res[annot_int$cell_id,]

for(j in 1:ncol(anocol))
{
  png(file.path(resdir_UMAP,paste0("UMAP_",colnames(anocol)[j],".png")), height=1500,width=1500,res=300) 
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

indice <- which(annotCol=="cons_BC_lenti")
if(length(indice) >0 ){
  png(file.path(resdir_UMAP,paste0("UMAP_",colnames(anocol)[indice],".png")), height=1150,width=950,res=350) 
  plot((umap_res[!is.na(annot_int$cons_BC_lenti),]),
       col=alpha(anocol[!is.na(annot_int$cons_BC_lenti),indice],0.5),
       pch=20,cex=0.4,main=paste0(colnames(anocol)[indice]," perplexity=30"),xlab="component 1",ylab="component 2",ylim=c(min(umap_res[,2]),max(umap_res[,2])))
  dev.off()
  
  png(file.path(resdir_UMAP,paste0("exp3_UMAP_",colnames(anocol)[indice],".png")),  height=1150,width=950,res=350)
  plot((umap_res[!is.na(annot_int$cons_BC_lenti) & annot_int$sample_id %in% c("MM468_chemonaive","MM468_5FU3_day50","MM468_5FU3_day77","MM468_5FU3_day202"),]), col=alpha(anocol[!is.na(annot_int$cons_BC_lenti) & annot_int$sample_id %in% c("MM468_chemonaive","MM468_5FU3_day50","MM468_5FU3_day77","MM468_5FU3_day202"),indice],0.8),pch=20,cex=0.4,main=paste0(colnames(anocol)[indice]," perplexity=30"),
       xlab="component 1",ylab="component 2",ylim=c(min(umap_res[,2]),max(umap_res[,2])), xlim=c(min(umap_res[,1]),max(umap_res[,1])))
  dev.off()
  
  png(file.path(resdir_UMAP,paste0("exp5_UMAP_",colnames(anocol)[indice],".png")),  height=1150,width=950,res=350) 
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
    filter(barcode>15) 
  
  
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

################################################################################################################
# Hierarchical clustering to group samples in complement to Louvain clustering on a subset of cells ############
################################################################################################################
#take a subset, do ONCE
anocol <- as.data.frame(anocol)
rownames(anocol) <- annot_int$cell_id

set.seed(47)
mati <- t(as.data.frame(pca_object))
hc_PCA <- hclust(as.dist(1 - cor(mati)), method = "ward.D")
mat.so.cor <- mati[,hc_PCA$order]

annot_int$hierarchical_cluster = paste0("C", cutree(hc_PCA,k = 7))
anocol. = annotToCol2(annot_int) 
anocol[,"hierarchical_cluster"] = anocol.[,"hierarchical_cluster"]

png(file.path(resdir_heatmaps,"Clustering_correlation_matrix_PCA_subset.png"), height=1500,width=1500,res=300)
geco.hclustAnnotHeatmapPlot(x=cor(mat.so.cor),
                            hc=hc_PCA,
                            hmColors=corColors,
                            anocol=as.matrix(anocol[hc_PCA$order,]),#[,ncol(cc.col):1]
                            xpos=c(0.15,0.9,0.164,0.885),
                            ypos=c(0.1,0.5,0.5,0.6,0.62,0.95),
                            dendro.cex=0.1,
                            xlab.cex=0.2,
                            hmRowNames=FALSE,
                            hmRowNames.cex=0.01
)
dev.off()

#Plot repartition of Expression groups with time
corres_day <- data.frame(Sample=sample_UNC_5FU_study,
                         day=gsub(".*day","",sample_UNC_5FU_study))
corres_day[which(corres_day$Sample=="MM468_chemonaive"),"day"] = 0
corres_day$day = as.numeric(corres_day$day)
annot_int$day <- as.integer(corres_day$day[match(annot_int$sample_id,corres_sample$Sample)])

annot_int$col <- annot_int[,"louvain_partition"]
table(annot_int$louvain_partition,annot_int$col)
test <- annot_int %>% mutate(
  louvain_partition=factor(louvain_partition, levels = paste0("C",1:4))) %>%
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
exp.sd <- apply(UNC_LogCounts,1,sd)
sel <- order(exp.sd,decreasing=T)[1:NbGenes]
mat <- UNC_LogCounts[sel,];dim(mat)
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

  