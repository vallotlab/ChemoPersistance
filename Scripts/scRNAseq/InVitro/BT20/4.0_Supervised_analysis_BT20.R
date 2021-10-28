library(here)

maindir= here()
resdir <- file.path(maindir,"output","scRNAseq","BT20","Persister","Supervised")
if(!dir.exists(resdir)) dir.create(resdir)
tabdir <- file.path(resdir,"Tables");if(!file.exists(tabdir)){dir.create(tabdir)}
plotdir <- file.path(resdir,"Plots");if(!file.exists(plotdir)){dir.create(plotdir)}

resdir_boxplots <-  file.path(resdir,"Boxplots");if(!file.exists(resdir_boxplots)){dir.create(resdir_boxplots)}
resdir_UMAPs <-  file.path(resdir,"UMAP");if(!file.exists(resdir_UMAPs)){dir.create(resdir_UMAPs)}
resdir_heatmaps <-  file.path(resdir,"Heatmaps");if(!file.exists(resdir_heatmaps)){dir.create(resdir_heatmaps)}

EnrichDir <- file.path(tabdir,"Gene_set_analysis");if(!file.exists(EnrichDir)){dir.create(EnrichDir)}
RDataSupdir <-  file.path(resdir,"RData");if(!file.exists(RDataSupdir)){dir.create(RDataSupdir)}
RDatadir <- file.path(file.path(maindir,"output","scRNAseq","BT20","Persister","Unsupervised","RData"))

source(file.path(maindir,"Scripts","global_var.R"))

## Import annotation file
load(file.path(RDatadir,"BT20.RData"))
metadata <- as.data.frame(annot_int)

####################################################################
### COMPARISON SETUP  
####################################################################
# log2FC_thresholds <- log2(3)
Signif_threshold <- 0.01

## Type of analysis
algoType <- c("edgeR") ## Limma or DESeq or a vector of both or edgeR for single-cell
useBlockFactor <- FALSE ## TRUE or FALSE
annotationDatabases <- c("MSigDB") # "KEGG" or "GO" or "KEGG_GO" or "MSigDB"

table(metadata[,c("sample_id","louvain_partition")])
mygps <- list(
    'persister'=metadata[which(metadata$sample_id %in% c("BT20_persister") ),"cell_id"]
)
myrefs <- list(
    'chemonaive'=metadata[which(metadata$sample_id %in% c("BT20_chemonaive") ),"cell_id"]
) # Sets reference (can be 1 or more samples)

refs <- names(myrefs) 
groups <- names(mygps) 

##############################################################################################
###############			LOAD FILES			    ##############################
##############################################################################################

##raw counts for edgeR
RawCounts <- Signal[,metadata$cell_id]
feature <- gene_metadata
rm(Signal); gc()

##############################################################################################
###############			START COMPARISONS		 ##############################
##############################################################################################

##############################################################################################
###############			0.PRE-PROCESS			 ##############################
##############################################################################################
# select genes which are detected in at least 1% of cells
log2FC_thresholds = log2(3)
Percent=0.01
sel <- which(apply(RawCounts,1,function(x) length(which(x>0)))>Percent*dim(RawCounts)[2])
RawCounts <- RawCounts[sel,];feature <- feature[sel,]
gc()

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
        under_res = my.res %>% dplyr::filter(log2FC.persister < -log2FC_threshold & qval.persister < Signif_threshold) %>% 
            dplyr::arrange(qval.persister) %>% dplyr::select(Symbol,log2FC.persister,qval.persister)
        over_res = my.res %>% dplyr::filter(log2FC.persister > log2FC_threshold & qval.persister < Signif_threshold) %>% 
            dplyr::arrange(qval.persister) %>% dplyr::select(Symbol,log2FC.persister,qval.persister)
        
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
            #ind <- which(enrich.test$`q-value`<= 0.1);if(!length(ind)){ind <- 1:20}
            #Overexpressed  <- enrich.test[ind,]		}
            Overexpressed  <- enrich.test
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
}
##########################
# Plot UMAPs #############
load(file=file.path(RDatadir,"anocol.RData"))
load(file=file.path(RDatadir,"persister_gene_cell_annot.RData"))
load(file=file.path(RDatadir,"umap_persister.RData"))
load(file=file.path(RDatadir,"LogCounts.RData"))
load(file.path(RDataSupdir,"Supervised_res_object_edgeR.Rdata"))

diff_genes <-  c("NNMT","TAGLN","INHBA","KRT6A","KRT6B","KRT16","KRT14","KLK10",
                 "KLK5","KLK7")
diff_genes_figures <-  c("KRT17", "INHBB", "KRT5","KRT8","BMP6","CDKN2B","TGFBR2",
                         "TGFBR3","TGFB2","INHBA","INHBB","SMAD2","TGFB1",
                         "FOSL1","NNMT","KRT14","KLK10","KLK5","TAGLN","NNMT",
                         "ELN","TGFBR3","MIF", "ABCC4","ABCC3","CD24","VIM","CDH2",
                         "CDH1","ABCA5","LBH")

pcaText <- FALSE
annotText <- "sample_id"

gene_lists <- c("KEGG_TGF_BETA_SIGNALING_PATHWAY",
                "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
                "HALLMARK_APOPTOSIS","HALLMARK_TNFA","LEI_MYB_TARGETS")

annotCol = unique(c(colnames(annot_int),diff_genes,diff_genes_figures))

for(log2FC_threshold in log2FC_thresholds){
    overgenes = my.res$Symbol[which(my.res$log2FC.persister > log2FC_threshold & 
                                        my.res$qval.persister < Signif_threshold)]
    for(i in gene_lists) {
        genes <- intersect(overgenes, unlist(str_split(MSIG.gs[which(MSIG.gs$Gene.Set==i),3],",")))
        annot_int$test <- apply(LogCounts[which(rownames(LogCounts) %in% genes),],2,mean)
        colnames(annot_int)[dim(annot_int)[2]] <- paste0(i,"_logFC",round(log2FC_threshold,2))
        annotCol = c(annotCol,colnames(annot_int)[dim(annot_int)[2]])
    }
}

for(i in diff_genes_figures){
    if(i %in% rownames(LogCounts))
        annot_int[,i] <- LogCounts[which(rownames(LogCounts)==i),]
}
gc()

anocol2 <- geco.unsupervised::geco.annotToCol4(annotS=annot_int[,intersect(annotCol,colnames(annot_int))],plotLegend = F, scale_q = "inferno")
for(i in intersect(colnames(anocol),colnames(anocol2))) anocol2[,i] = anocol[,i]

png(file.path(resdir_boxplots,"Boxplot_KEGG_TGFB_FC2.png"),width=1500,height=1500,res=300)
boxplot(annot_int$KEGG_TGF_BETA_SIGNALING_PATHWAY_logFC1.58~annot_int$sample_id,las=2)
dev.off()

png(file.path(resdir_boxplots,"Boxplot_KEGG_TGFB_FC1.png"),width=1500,height=1500,res=300)
boxplot(annot_int$KEGG_TGF_BETA_SIGNALING_PATHWAY_logFC1.58~annot_int$sample_id,las=2)
dev.off()

png(file.path(resdir_boxplots,"Boxplot_INHBA.png"),width=1500,height=1500,res=300)
boxplot(annot_int$INHBB~annot_int$sample_id,las=2)
dev.off()

save(anocol2,annot_int,file=file.path(RDataSupdir,"annot_anocol_final.RData")) 

for(i in setdiff(colnames(anocol2),colnames(anocol)))
{
    j = which(colnames(anocol2) == i)
    png(file.path(resdir_UMAPs,paste0("UMAP_",colnames(anocol2)[j],".png")), height=1350,width=1200,res=300)
    if(class(annot_int[1,colnames(anocol2)[j]])=="numeric"){
        plot((umap_res), col=alpha(anocol2[,j],0.3),pch=20,cex=0.6,
             main=paste0(colnames(anocol2)[j]," min=",round(min(annot_int[,colnames(anocol2)[j]]),digits=3)," max=",
                         round(max(annot_int[,colnames(anocol2)[j]]),digits=3)),
             xlab="component 1",ylab="component 2")} else {
                 plot((umap_res), col=alpha(anocol2[,j],0.3),pch=20,cex=0.6,
                      main=paste0(colnames(anocol2)[j]),
                      xlab="component 1",ylab="component 2")
                 if(colnames(anocol2)[j]=="sample_id"){
                     #Plot persister with add plot since there is too low cell number
                     pers = which(annot_int$sample_id=="HBCx95_persister_6souris")
                     points(umap_res[pers,], col=alpha(anocol2[pers,j],0.5),pch=20,cex=0.6)
                     
                 }
             }
    
    
    dev.off()
}

################################################################################
###Generate heatmaps with gene lists 
################################################################################
#load "subset_analysis.RData" if starting here
#import my.res et gene lists from diff analysis, logCounts, take top 10

how_many_top <- 15
for(log2FC_threshold in log2FC_thresholds){
    load(file.path(RDataSupdir,paste0("Enrichment_test_",groups,"_vs_",refs,"_",
                                      annotbase,"_logFC",round(log2FC_threshold,2),".RData")))
    
    significPathway <- Overexpressed[Overexpressed$Class %in% c("c2_curated","c5_GO","hallmark"),]
    significPathway$Deregulated_genes <- as.character(significPathway$Deregulated_genes)
    significPathway <- significPathway[1:how_many_top,]
    significPathway_breast <- rbind(
        significPathway[grep(pattern = "BREAST", significPathway$Gene.Set),],
        significPathway[grep(pattern = "MAMMARY", significPathway$Gene.Set),],
        significPathway[grep(pattern = "HALLMARK", significPathway$Gene.Set),])
    
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
    
    rowClust <- hclust(dist((mat.so)), method = "ward.D")
    png(file.path(plotdir,paste0(
        "Clustering_",how_many_top,"Pathways_PersistersAll_",distHC,"_",methHC,
        "log2FC",round(log2FC_threshold,2),".png")), 
        height=2000,width=3000,res=300)
    par(oma=c(2,5,3,7))
    geco.hclustAnnotHeatmapPlot(x=(mat.so[rowClust$order,]),
                                hc=hc,
                                hmColors=hmColors,
                                anocol=as.matrix(anocol_subset[hc$order,c(1,2,5)]),#[,ncol(cc.col):1]
                                xpos=c(0.15,0.9,0.164,0.885),
                                ypos=c(0.1,0.5,0.5,0.6,0.62,0.95),
                                dendro.cex=0.4,
                                xlab.cex=0.4,
                                hmRowNames=TRUE,
                                hmRowNames.cex=0.4
    )
    dev.off()
}


##### Integrative heatmap of pathway enrichment###############


