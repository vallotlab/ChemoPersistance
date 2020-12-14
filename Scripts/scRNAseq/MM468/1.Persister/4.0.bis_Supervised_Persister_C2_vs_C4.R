library(here)

maindir= here()
resdir <- file.path(maindir,"output","scRNAseq","MM468","Persister","Supervised")
tabdir <- file.path(resdir,"Tables");if(!file.exists(tabdir)){dir.create(tabdir)}
plotdir <- file.path(resdir,"Plots");if(!file.exists(plotdir)){dir.create(plotdir)}

resdir_boxplots <-  file.path(resdir,"Boxplots");if(!file.exists(resdir_boxplots)){dir.create(resdir_boxplots)}
resdir_UMAPs <-  file.path(resdir,"UMAP");if(!file.exists(resdir_UMAPs)){dir.create(resdir_UMAPs)}
resdir_heatmaps <-  file.path(resdir,"Heatmaps");if(!file.exists(resdir_heatmaps)){dir.create(resdir_heatmaps)}

EnrichDir <- file.path(tabdir,"Gene_set_analysis");if(!file.exists(EnrichDir)){dir.create(EnrichDir)}
RDataSupdir <-  file.path(resdir,"RData");if(!file.exists(RDataSupdir)){dir.create(RDataSupdir)}
RDatadir <- file.path(file.path(maindir,"output","scRNAseq","MM468","Persister","Unsupervised","RData"))

source(file.path(maindir,"Scripts","global_var.R"))

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
annotationDatabases <- c("MSigDB") # "KEGG" or "GO" or "KEGG_GO" or "MSigDB"

table(metadata[,c("sample_id","louvain_partition")])
mygps <- list(
    'C2_pers'=metadata[which(metadata$louvain_partition %in% c("C2") ),"cell_id"]
)
myrefs <- list(
    'C4_pers'=metadata[which(metadata$louvain_partition %in% c("C4") ),"cell_id"]
) # Sets reference (can be 1 or more samples)

refs <- names(myrefs) 
groups <- names(mygps) 

##############################################################################################
###############			LOAD FILES			        			##############################
##############################################################################################

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
            dplyr::arrange(qval.C2_pers) %>% dplyr::select(Symbol,log2FC.C2_pers,qval.C2_pers)
        over_res = my.res %>% dplyr::filter(log2FC.C2_pers > log2FC_threshold & qval.C2_pers < Signif_threshold) %>% 
            dplyr::arrange(qval.C2_pers) %>% dplyr::select(Symbol,log2FC.C2_pers,qval.C2_pers)
        
        WriteXLS(c("over_res","under_res","my.res"),
                 ExcelFileName = file.path(
                     tabdir,paste0("Differential_analysis_Limma_logFC_C2_vs_C4_",
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
            tabdir,paste0("Number_differentially_expressed_genes_per_group_logFC_C2_vs_C4",
                          round(log2FC_threshold,2),".csv")),row.names=F,quote=F,sep=";")
    } ## end of if edgeR
    save(my.res, file=file.path(RDataSupdir,paste("Supervised_res_object_edgeR_C2_vs_C4.Rdata",sep="")))
    print("summaryedgeRCompare Done")
}

##############################################################################################
###############			3.ENRICHMENT ANALYSIS					##############################
##############################################################################################
load(file.path(RDataSupdir,paste("Supervised_res_object_edgeR_C2_vs_C4.Rdata",sep="")))
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


