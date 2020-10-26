
##############################################################################################
###############			LIBRARIES AND FUNCTIONS					##############################
##############################################################################################
options(stringsAsFactors=FALSE, width=180)
organism <- "hg38"
GencodeVersion <- ""
library(geco.supervised);  library(edgeR); library(ggplot2); library(rgl); library(RColorBrewer); library(genefilter); library(xtable); library(geco.RNAseq); library(WriteXLS); library(data.table); library(stringr); library(limma); library(edgeR);library(dplyr) ####################################################################
library(scTools); library(WriteXLS)
#library(tidyverse)
### PARAMETERS	
####################################################################


res_folder_name <- "All_clusters_vs_other_clusters_edgeR"
useClusterInfoFromUnsupp <- TRUE ## TRUE , FALSE
####################################################################
### DIRECTORIES and FILES 
####################################################################

######VARIABLES
options(width=180)
maindir <- "/media/pprompsy/LaCie/InstitutCurie/Z_server/Manuscripts/2020_ChemoTolerance/Raw_analysis/scRNAseq/Results_MM468/"
resdir <- "/media/pprompsy/LaCie/InstitutCurie/Z_server/Manuscripts/2020_ChemoTolerance/Raw_analysis/scRNAseq/Results_MM468/Supervised_UNC";if(!file.exists(resdir)){dir.create(resdir)}
tabdir <- file.path(resdir,"Tables");if(!file.exists(tabdir)){dir.create(tabdir)}
plotdir <- file.path(resdir,"Plots");if(!file.exists(plotdir)){dir.create(plotdir)}
EnrichDir <- file.path(tabdir,"Gene_set_analysis");if(!file.exists(EnrichDir)){dir.create(EnrichDir)}
RDataSupdir <-  file.path(resdir,"RData");if(!file.exists(RDataSupdir)){dir.create(RDataSupdir)}
RDatadir <- "/media/pprompsy/LaCie/InstitutCurie/Z_server/Manuscripts/2020_ChemoTolerance/Raw_analysis/scRNAseq/Results_MM468/Unsupervised_UNC/RData/"
## Import annotation file
load(file.path(RDatadir,"subset_analysis.RData"))
#load(file.path(RDatadir,"MM468.RData"))
metadata <- as.data.frame(annot_subset)

####################################################################
### COMPARISON SETUP  
####################################################################
FC_threshold <- 2
Signif_threshold <- 0.01
## Type of analysis
algoType <- c("edgeR") ## Limma or DESeq or a vector of both or edgeR for single-cell
useBlockFactor <- FALSE ## TRUE or FALSE
#BlockFactor <- "Sex"
annotationDatabases <- c("MSigDB") # "KEGG" or "GO" or "KEGG_GO" or "MSigDB"

table(metadata[,c("sample_id","louvain_partition")])
mygps <- list(
  'C1_MM468_UNC_day33_a'=metadata[which(metadata$louvain_partition == "C1" ),"cell_id"],
  'C1_MM468_UNC_day33_b'=metadata[which(metadata$louvain_partition == "C1" ),"cell_id"],
  'C2_persister_5FU_UNC'=metadata[which(metadata$louvain_partition == "C2"),"cell_id"],
  'MM468_5FU_UNC'=metadata[which(metadata$sample_id == "MM468_5FU6_UNC_day33"),"cell_id"]
)

myrefs <- list(
'C3_DMSO_a'=metadata[which(metadata$louvain_partition == "C3" ),"cell_id"],
'C2_persister_5FU_UNC'=metadata[which(metadata$louvain_partition == "C2"),"cell_id"],
'C3_DMSO_b'=metadata[which(metadata$louvain_partition == "C3" ),"cell_id"],
'persister_non_UNC'=metadata[which(metadata$sample_id %in% c("MM468_5FU3_day50","MM468_5FU3_day77","MM468_5FU5_day67","MM468_5FU6_day33") ),"cell_id"]
) # Sets reference (can be 1 or more samples)


refs <- names(myrefs) 
groups <- names(mygps) 


##############################################################################################
###############			LOAD FILES			        			##############################
##############################################################################################
if (organism == "mm10") {
MSigDBFile <- "~/Documents/bioinfo/GeCo.Annotation/msigdb_v5.0_GMTs/MSIG_v5_mousified.RData"
KEGGFile <- "~/Documents/bioinfo/GeCo.Annotation/mm10/gs.KEGG.toSymbol.mm.RData"
KEGGdefFile <- "~/Documents/bioinfo/GeCo.Annotation/mm10/gs.KEGG.def.RData" # object KEGGdef for annotating path_id with path_name
GOFile <- "~/Documents/bioinfo/GeCo.Annotation/mm10/gs.GO.toSymbol.mm_2012-10-03.RData" }

if (organism == "hg19") {
MSigDBFile <- "~/Documents/bioinfo/GeCo.Annotation/msigdb_v5.0_GMTs/MSIG_v5.RData" }
if (organism == "hg38") {
MSigDBFile1 <- file.path(maindir,"annotation","hg38.MSIG.gs.rda")
MSigDBFile2 <- file.path(maindir,"annotation","hg38.MSIG.ls.rda")
}



# load(file.path(RDatadir,"../../Unsupervised_persister/RData/MM468.RData") )
#load(file.path(RDatadir,"NormCounts.RData"))

##raw counts for edgeR
# RawCounts <- Signal[,metadata$cell_id] #Matrice de Signal brute MM468.RData
feature <- gene_metadata
# rm(Signal); gc()
##alternative for normalized data
# Counts <- as.matrix(NormCounts)
# dim(Counts)
# feature <- gene_metadata[rownames(Counts),]
# dim(feature)

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
# save(RawCounts,file = file.path(RDatadir,"RawCounts_UNC.RData"))
load(file.path(RDatadir,"RawCounts_UNC.RData"))

dim(RawCounts)
dim(feature)
# sel <- which(apply(Counts,1,function(x) length(which(x>0)))>Percent*dim(Counts)[2])
# Counts <- Counts[sel,];feature <- feature[sel,]
# dim(RawCounts)
# dim(feature)



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
                                          fc.th=FC_threshold,
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
  
  #write.table(my.res, file.path(tabdir,"Supervised_analysis_EdgeR.csv"),row.names=F,quote=F,sep=";" )
  WriteXLS(my.res, ExcelFileName = file.path(tabdir,"Supervised_analysis_EdgeR.xls"), perl = "perl",
           verbose = FALSE, row.names = FALSE, col.names = TRUE, AdjWidth = TRUE, AutoFilter = TRUE,
           BoldHeaderRow = TRUE, na = "", FreezeRow = 1, FreezeCol = 1)
  
  save(my.res, file=file.path(RDataSupdir,paste("Supervised_res_object_edgeR.Rdata",sep="")))
  
  summaryTab <- geco.summaryCompareedgeR(restab=my.res,
                                         ref=myrefs,
                                         groups=mygps,
                                         qval.th=Signif_threshold,
                                         fc.th=FC_threshold,
                                         plotdir=plotdir
  )
  
  # summaryTab <- data.frame("."=row.names(summaryTab), summaryTab, check.names=FALSE)
  # write.table(summaryTab, file.path(tabdir,"Number_differentially_expressed_genes_per_group.csv"),row.names=F,quote=F,sep=";"  )
  print("summaryedgeRCompare Done")
  
} ## end of if edgeR



##############################################################################################
###############			3.ENRICHMENT ANALYSIS					##############################
##############################################################################################
# Load annotation databases files
if (organism == "hg19" | organism == "hg38") { 
  load(MSigDBFile1)
  load(MSigDBFile2)
  MSIG.ls = hg38.MSIG.ls
  MSIG.gs = hg38.MSIG.gs
  }
if (organism == "mm10") {load(MSigDBFile); load(KEGGFile) ;load(KEGGdefFile) ; load(GOFile); KEGGdef <- KEGGdef[,1:2]; colnames(KEGGdef) <- c("path_id","path_name")}

rownames(my.res) <- my.res$Symbol

	annotbase <- "MSigDB" 
	database <- MSIG.ls ##MSigDB
	
	Overexpressed  <- Underexpressed <- data.frame()
  
	#my.res$Gene<- sub("hg19_","",my.res$Gene)
	
	#rownames(my.res) <- my.res$Gene
	reflist <- unique(my.res$Symbol);length(reflist)

	for(i in 1:length(groups)) 	{
		gp <- groups[i]
  		if (length(refs)>1) {ref <- refs[i]} else {ref <- refs[1]}
		print(paste0("Processing ",gp, " vs ", ref, " _ ",annotbase ))

			signific <- which(my.res[,paste("qval",gp,sep=".")] <= Signif_threshold & abs(my.res[,paste("log2FC",gp,sep=".")]) > FC_threshold)
			over <- which(my.res[,paste("qval",gp,sep=".")] <= Signif_threshold & my.res[,paste("log2FC",gp,sep=".")] > FC_threshold)
			under <- which(my.res[,paste("qval",gp,sep=".")] <= Signif_threshold & my.res[,paste("log2FC",gp,sep=".")] < -FC_threshold)
			print(paste0("significant = ", length(signific))) ; print(paste0("over = ", length(over))) ; print(paste0("under = ", length(under)))
			
			
			 if(length(over)){
			   enrich.test <- geco.enrichmentTest(gene.sets=database,mylist=my.res$Symbol[over],possibleIds=reflist)
			   enrich.test <- data.frame(Gene_set_name=rownames(enrich.test), enrich.test, check.names=FALSE)
			   enrich.test <- merge( subset(MSIG.gs, select=-Genes), enrich.test, by.x="Gene.Set", by.y="Gene_set_name", all.y=TRUE, sort=FALSE ) ## Get class of gene set
			   
			   enrich.test <- enrich.test[order(enrich.test$`p-value`),]
			   
# 			if (annotbase=="MSigDB") {
#             enrich.test <- merge( subset(MSIG.gs, select=-Genes), enrich.test, by.x="Gene.Set", by.y="Gene_set_name", all.y=TRUE, sort=FALSE ) ## Get class of gene set
#             enrich.test <- enrich.test %>% dplyr::rename(Gene_set_type = Class, Gene_set_name = Gene.Set) %>% 
#                                 dplyr::mutate(Gene_set_type = str_replace(Gene_set_type, "c[0-9]_", "")) %>% ## clean class type
#                                 mutate(Gene_set_type = str_to_title(Gene_set_type) ) %>% ## uppercase first letter
#                                 mutate(Gene_set_type =  dplyr::str_replace(Gene_set_type, "Go", "GO")) ## teat special "GO" case
#             }
			#if (annotbase=="KEGG") {enrich.test <-merge( KEGGdef, enrich.test, by.x="id",by.y="Gene_set_name", all.y=TRUE )}
			enrich.test <- enrich.test[order(enrich.test$`p-value`),]
			ind <- which(enrich.test$`q-value`<= 0.1);if(!length(ind)){ind <- 1:20}
				Overexpressed  <- enrich.test[ind,]		}
			
			 if(length(under)){
			enrich.test <- geco.enrichmentTest(gene.sets=database,mylist=my.res$Symbol[under],possibleIds=reflist)
            enrich.test <- data.frame(Gene_set_name=rownames(enrich.test), enrich.test, check.names=FALSE)
            enrich.test <- merge( subset(MSIG.gs, select=-Genes), enrich.test, by.x="Gene.Set", by.y="Gene_set_name", all.y=TRUE, sort=FALSE ) ## Get class of gene set
            
# 			if (annotbase=="MSigDB") {
#             enrich.test <- merge( subset(MSIG.gs, select=-Genes), enrich.test, by.x="Gene.Set", by.y="Gene_set_name", all.y=TRUE, sort=FALSE ) ## Get class of gene set
#             enrich.test <- enrich.test %>% 
#                                 rename(Gene_set_type = Class, Gene_set_name = Gene.Set) %>% 
#                                 mutate(Gene_set_type =  str_replace(Gene_set_type, "c[0-9]_", "")) %>% ## clean class type
#                                 mutate(Gene_set_type = str_to_title(Gene_set_type) ) %>% ## uppercase first letter
#                                 mutate(Gene_set_type =  str_replace(Gene_set_type, "Go", "GO")) ## teat special "GO" case
#             }
			# if (annotbase=="KEGG") {enrich.test <-merge( KEGGdef, enrich.test, by="Gene_set_name", all.y=TRUE )}
			enrich.test <- enrich.test[order(enrich.test$`p-value`),]
			ind <- which(enrich.test$`q-value`<= 0.1);if(!length(ind)){ind <- 1:20}
				Underexpressed <- enrich.test[ind,]		}
				
	WriteXLS(c("Overexpressed", "Underexpressed"), ExcelFileName = file.path(EnrichDir,paste0("Enrichment_test_",gp,"_vs_",ref,"_",annotbase,".xlsx")),
	         SheetNames = c(substr(paste0("Overexp_in_", gp),0,30), substr(paste0("Underexp_in_", gp),0,30) ), perl = "perl", verbose = FALSE, row.names = FALSE, col.names = TRUE, AdjWidth = TRUE, AutoFilter = TRUE, BoldHeaderRow = TRUE, na = "", FreezeRow = 1, FreezeCol = 1)
	save(Overexpressed,file=file.path(RDataSupdir,paste0("Overexpressed_",gp,"_vs_",ref,".RData")))
	} # end of groups loop

