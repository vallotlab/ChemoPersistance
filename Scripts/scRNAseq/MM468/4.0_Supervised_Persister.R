
##############################################################################################
###############			LIBRARIES AND FUNCTIONS					##############################
##############################################################################################
options(stringsAsFactors=FALSE, width=180)
organism <- "hg38"
GencodeVersion <- ""

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
library(gplots)
library(ggplot2)
library(rgl)
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
library(dendextend)

#Geco packages
library(geco.supervised)
library(geco.RNAseq)
library(geco.utils)
library(geco.visu)
library(scTools)
library(geco.unsupervised)

### PARAMETERS	
####################################################################
useClusterInfoFromUnsupp <- TRUE ## TRUE , FALSE

####################################################################
### DIRECTORIES and FILES 
####################################################################

###### PATHS
options(width=180)

maindir= here()
resdir <- file.path(maindir,"output","scRNAseq","MM468","Persister","Supervised")
tabdir <- file.path(resdir,"Tables");if(!file.exists(tabdir)){dir.create(tabdir)}
plotdir <- file.path(resdir,"Plots");if(!file.exists(plotdir)){dir.create(plotdir)}
EnrichDir <- file.path(tabdir,"Gene_set_analysis");if(!file.exists(EnrichDir)){dir.create(EnrichDir)}
RDataSupdir <-  file.path(resdir,"RData");if(!file.exists(RDataSupdir)){dir.create(RDataSupdir)}
RDatadir <- file.path(file.path(maindir,"output","scRNAseq","MM468","Persister","Unsupervised","RData"))

## Import annotation file
load(file.path(RDatadir,"subset_analysis.RData"))
load(file.path(RDatadir,"MM468.RData"))
metadata <- as.data.frame(annot_subset)

####################################################################
### COMPARISON SETUP  
####################################################################
log2FC_thresholds <- log2(c(2,3,4))
Signif_threshold <- 0.01

## Type of analysis
algoType <- c("edgeR") ## Limma or DESeq or a vector of both or edgeR for single-cell
useBlockFactor <- FALSE ## TRUE or FALSE
#BlockFactor <- "Sex"

annotationDatabases <- c("MSigDB") # "KEGG" or "GO" or "KEGG_GO" or "MSigDB"

table(metadata[,c("sample_id","louvain_partition")])
mygps <- list(
'C2_pers'=metadata[which(metadata$louvain_partition %in% c("C2") ),"cell_id"]
)
myrefs <- list(
'DMSO'=metadata[which(metadata$sample_id %in% c("MM468_initial") ),"cell_id"]
) # Sets reference (can be 1 or more samples)

refs <- names(myrefs) 
groups <- names(mygps) 

##############################################################################################
###############			LOAD FILES			        			##############################
##############################################################################################
MSigDBFile1 <- file.path(maindir,"annotation","hg38.MSIG.gs.rda")
MSigDBFile2 <-file.path(maindir,"annotation","hg38.MSIG.ls.rda")

# load(file=file.path(RDataSupdir,paste("Supervised_res_object_edgeR.Rdata",sep="")))
load(MSigDBFile1)
load(MSigDBFile2)
MSIG.ls = hg38.MSIG.ls
MSIG.gs = hg38.MSIG.gs

##raw counts for edgeR
RawCounts <- Signal[,metadata$cell_id]
feature <- gene_metadata
rm(Signal); gc()

##############################################################################################
###############			START COMPARISONS						##############################
##############################################################################################

##############################################################################################
###############			0.PRE-PROCESS				##############################
##############################################################################################
# select genes which are detected in at least 1% of cells
Percent=0.01
sel <- which(apply(RawCounts,1,function(x) length(which(x>0)))>Percent*dim(RawCounts)[2])
RawCounts <- RawCounts[sel,];feature <- feature[sel,]
gc()

##############################################################################################
###############			1.Wilcox                     			##############################
##############################################################################################


if ("wilcox" %in% algoType) {
    
    #en entrée matrice de compte non normalisée, ou voir scran pour normalisation 
    my.res <- geco.CompareWilcox (dataMat=Counts,
                                  annot=metadata,
                                  ref=myrefs,
                                  groups=mygps,
                                  featureTab=feature
    )
    
    write.table(my.res, file.path(tabdir,"Supervised_analysis_Wilcox.csv"),row.names=F,quote=F,sep=";",dec="," )
    save(my.res, file=file.path(RDataSupdir,paste("Supervised_res_Wilcox_object.Rdata",sep="")))
    
    summaryTab <- geco.summaryCompareWilcox(restab=my.res,
                                            ref=myrefs,
                                            groups=mygps,
                                            qval.th=Signif_threshold,
                                            fc.th=log2FC_threshold,
                                            plotdir=plotdir
    )
    
    summaryTab <- data.frame("."=row.names(summaryTab), summaryTab, check.names=FALSE)		
    write.table(summaryTab, file.path(tabdir,"Number_differentially_bound_peaks_per_group.csv"),row.names=F,quote=F,sep=";",dec=","  )
    print("summaryWilcoxRCompare Done")
    
} ## end of if wilcox

##############################################################################################
###############			2.EdgeR                     			##############################
##############################################################################################
if ("edgeR" %in% algoType) {
    
    
    
    #en entrée matrice de compte non normalisée, ou voir scran pour normalisation 
    my.res <- geco.CompareedgeRGLM (dataMat=RawCounts,
                                    annot=metadata,
                                    ref=myrefs,
                                    groups=mygps,
                                    featureTab=feature
    )
    
    for(log2FC_threshold in log2FC_thresholds){
        under_res = my.res %>% dplyr::filter(log2FC.C2_pers < -log2FC_threshold & qval.C2_pers < Signif_threshold) %>% 
            dplyr::arrange(qval.C2_pers) %>% select(Symbol,log2FC.C2_pers,qval.C2_pers)
        over_res = my.res %>% dplyr::filter(log2FC.C2_pers > log2FC_threshold & qval.C2_pers < Signif_threshold) %>% 
            dplyr::arrange(qval.C2_pers) %>% select(Symbol,log2FC.C2_pers,qval.C2_pers)
        
        WriteXLS(c("over_res","under_res","my.res"),
                 ExcelFileName = file.path(
                     tabdir,paste0("Differential_analysis_Limma_logFC_",
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
                                               plotdir=plotdir
        )
        
        write.table(summaryTab, file.path(
            tabdir,paste0("Number_differentially_expressed_genes_per_group_logFC_",
                          round(log2FC_threshold,2),".csv")),row.names=F,quote=F,sep=";")
    } ## end of if edgeR
    save(my.res, file=file.path(RDataSupdir,paste("Supervised_res_object_edgeR.Rdata",sep="")))
    print("summaryedgeRCompare Done")
}

##############################################################################################
###############			3.ENRICHMENT ANALYSIS					##############################
##############################################################################################
load(file.path(RDataSupdir,paste("Supervised_res_object_edgeR.Rdata",sep="")))
# Load annotation databases files
rownames(my.res) <- my.res$Symbol

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
        ExcelFileName = file.path(EnrichDir,
                                  paste0("Enrichment_test_",gp,"_vs_",ref,"_",annotbase,
                                         "_logFC",round(log2FC_threshold,2),".xlsx")), 
        SheetNames = c( paste0("Overexp_in_", gp), paste0("Underexp_in_", gp) ),
        perl = "perl", verbose = FALSE, row.names = FALSE, col.names = TRUE,
        AdjWidth = TRUE, AutoFilter = TRUE, BoldHeaderRow = TRUE, na = "",
        FreezeRow = 1, FreezeCol = 1)
    
    save(Overexpressed,Underexpressed, file = file.path(RDataSupdir,paste0("Enrichment_test_",gp,"_vs_",ref,"_",annotbase,"_logFC",round(log2FC_threshold,2),".RData")))
    }
}

##########################
# Plot UMAPs #############
load(file=file.path(RDatadir,"anocol.RData"))
load(file=file.path(RDatadir,"persister_LogCounts.RData"))
load(file=file.path(RDatadir,"gene_cell_annot_persister.RData"))

diff_genes <-  c("NNMT","TAGLN","INHBA","KRT6A","KRT6B","KRT16","KRT14","KLK10",
                 "KLK5","KLK7")
diff_genes_figures <-  c("KRT17", "INHBB", "KRT5","KRT8","BMP6","CDKN2B","TGFBR2",
                         "TGFBR3","TGFB2","INHBA","INHBB","SMAD2","TGFB1",
                         "FOSL1","NNMT","KRT14","KLK10","KLK5","TAGLN","NNMT",
                         "ELN","TGFBR3","MIF")


pcaText <- FALSE
annotText <- "sample_id"

gene_lists <- c("KEGG_TGF_BETA_SIGNALING_PATHWAY",
                "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
                "HALLMARK_APOPTOSIS","HALLMARK_TNFA","LEI_MYB_TARGETS")

annotCol = unique(c(colnames(annot_int), gene_lists,diff_genes,diff_genes_figures))
for(log2FC_threshold in log2FC_thresholds){
    overgenes = my.res$Symbol[which(my.res$log2FC.C2_pers > log2FC_threshold & 
                                        my.res$qval.C2_pers < Signif_threshold)]
    for(i in gene_lists) {
        genes <- intersect(overgenes, unlist(str_split(MSIG.gs[which(MSIG.gs$Gene.Set==i),3],",")))
        annot_int$test <- apply(persister_LogCounts[which(rownames(persister_LogCounts) %in% genes),],2,mean)
        colnames(annot_int)[dim(annot_int)[2]] <- paste0(i,"_logFC",round(log2FC_threshold,2))
    }
}

for(i in annotCol){
    if(i %in% rownames(persister_LogCounts))
        annot_int[,i] <- persister_LogCounts[which(rownames(persister_LogCounts)==i),]
}
gc()

if(file.exists(file.path(here(), "output","scRNAseq","common_over_genes_pers_vs_unt.RData"))){
    load(file.path(here(), "output","scRNAseq","common_over_genes_pers_vs_unt.RData"))
    for(i in common_overexpressed_genes){
        annot_int[,i] = persister_LogCounts[i,which(rownames(persister_LogCounts) %in% annot_int$cell_id)]
    }
    annotCol = unique(c(annotCol, common_overexpressed_genes))
}

anocol2 <- geco.unsupervised::geco.annotToCol4(annotS=annot_int,annotT=annot_int,
                                               plotLegend = F)


png(file.path(plotdir,"Boxplot_KEGG_TGFB_FC2.png"),width=1500,height=1500,res=300)
boxplot(annot_int$KEGG_TGF_BETA_SIGNALING_PATHWAY_logFC1~annot_int$sample_id,las=2)
dev.off()

png(file.path(plotdir,"Boxplot_KEGG_TGFB_FC1.png"),width=1500,height=1500,res=300)
boxplot(annot_int$KEGG_TGF_BETA_SIGNALING_PATHWAY_logFC1~annot_int$sample_id,las=2)
dev.off()

png(file.path(resdir_boxplots,"Boxplot_INHBA.png"),width=1500,height=1500,res=300)
boxplot(annot_int$INHBB~annot_int$sample_id,las=2)
dev.off()

save(anocol2,annot_int,file=file.path(RDataSupdir,"annot_anocol_final.RData")) 

################################################################################
###Generate heatmaps with gene lists 
################################################################################
#load "subset_analysis.RData" if starting here
#import my.res et gene lists from diff analysis, logCounts, take top 10

load(file.path(file.path("~/Desktop/scRNAseq_data_local/Results_MM468/Supervised/RData/Supervised_res_object_edgeR.Rdata"))
resdir <- "~/Desktop/scRNAseq_data_local/Results_MM468/Supervised";if(!file.exists(resdir)){dir.create(resdir)}

RDataSupdir <-  file.path(resdir,"RData");if(!file.exists(RDataSupdir)){dir.create(RDataSupdir)}

load(file=file.path(RDataSupdir,"Overexpressed_persisterall_vs_DMSO.RData"))
how_many_top <- 15
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
methHC <- c("ward","ward.D","ward.D2","single","complete","average")[2]
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

png(file.path(resdir,paste0("Clustering_",how_many_top,"Pathways_PersistersAll_",distHC,"_",methHC,".png")), 
    height=3000,width=2000,res=300)
par(oma=c(2,5,3,7))
geco.hclustAnnotHeatmapPlot(x=(mat.so[rowClust$order,]),
                            hc=hc,
                            hmColors=hmColors,
                            anocol=as.matrix(anocol_subset[hc$order,c(1:2,6)]),#[,ncol(cc.col):1]
                            xpos=c(0.15,0.9,0.164,0.885),
                            ypos=c(0.1,0.5,0.5,0.6,0.62,0.95),
                            dendro.cex=0.4,
                            xlab.cex=0.4,
                            hmRowNames=TRUE,
                            hmRowNames.cex=0.4
)
dev.off()




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


