
##############################################################################################
###############			LIBRARIES AND FUNCTIONS					##############################
##############################################################################################
options(stringsAsFactors=FALSE, width=180)
organism <- "hg38"
GencodeVersion <- ""
chRangeHM=T
distHC <- c("distPearson","distCosine","euclidean","maximum","manhattan","canberra","binary","minkowski")[1]
methHC <- c("ward","ward.D","ward.D2","single","complete","average")[2]

##ConsClust
repsCC <- 1000
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

library(geco.supervised);
library(geco.utils)
library(geco.visu)
library(edgeR); library(ggplot2); library(rgl); library(RColorBrewer); library(genefilter); library(xtable); library(geco.RNAseq); library(WriteXLS); library(data.table); library(stringr); library(limma); library(edgeR);library(dplyr) ####################################################################
library(scTools)
library(Matrix)
library(dplyr)
library(corrplot)
library(geco.utils)
library(geco.visu)
library(geco.unsupervised)
library(scatterplot3d)
library(scater)
library(Rtsne)
library(ccRemover)
library(colorRamps)
library(geco.supervised)
library(viridis)
library(colorRamps)
library(geco.supervised);  library(edgeR); library(ggplot2); library(rgl); library(RColorBrewer); library(genefilter); library(xtable); library(geco.RNAseq); library(WriteXLS); library(data.table); library(stringr); library(limma); library(edgeR)

library(scTools)
library(monocle3)
library(dplyr)
library(WriteXLS)
library(clValid)
library(ape)
library(ConsensusClusterPlus)
devtools::load_all("/media/pprompsy/LaCie/InstitutCurie/Documents/GitLab/R_packages/GeCo.Rpackages/geco.visu/")

MSigDBFile1 <- "/media/pprompsy/LaCie/InstitutCurie/Documents/GitLab/ChromSCape_devel/data/hg38.MSIG.gs.rda" 
MSigDBFile2 <- "/media/pprompsy/LaCie/InstitutCurie/Documents/GitLab/ChromSCape_devel/data/hg38.MSIG.ls.rda"
load(MSigDBFile1)
load(MSigDBFile2)
MSIG.ls = hg38.MSIG.ls
MSIG.gs = hg38.MSIG.gs
### PARAMETERS	
####################################################################


useClusterInfoFromUnsupp <- TRUE ## TRUE , FALSE
####################################################################
### DIRECTORIES and FILES 
####################################################################

######VARIABLES
options(width=180)
maindir <- "/media/pprompsy/LaCie/InstitutCurie/Z_server/Manuscripts/2020_ChemoTolerance/Raw_analysis/scRNAseq/"
resdir <- "/media/pprompsy/LaCie/InstitutCurie/Z_server/Manuscripts/2020_ChemoTolerance/Raw_analysis/scRNAseq/Results_PDX/Supervised/Co_occurence";if(!file.exists(resdir)){dir.create(resdir)}
RDatadir <- "/media/pprompsy/LaCie/InstitutCurie/Z_server/Manuscripts/2020_ChemoTolerance/Raw_analysis/scRNAseq/Results_PDX/Supervised/RData/"
RDatadir_unsup <- "/media/pprompsy/LaCie/InstitutCurie/Z_server/Manuscripts/2020_ChemoTolerance/Raw_analysis/scRNAseq/Results_PDX/Unsupervised//RData/"
RDataSupdir <-  file.path(resdir,"RData");if(!file.exists(RDataSupdir)){dir.create(RDataSupdir)}

FC_threshold <- 2
Signif_threshold <- 0.01

## Import annotation file
load(file.path(RDatadir,"Supervised_res_object_edgeR.Rdata"))
load(file.path(RDatadir_unsup,"annot_anocol_final.RData"))
overexpressed_pers_genes = my.res$gene_short_name[which(my.res$log2FC.persister_6_vs_UNT > FC_threshold & my.res$qval.persister_6_vs_UNT < Signif_threshold)]
# save(overexpressed_pers_genes, file=file.path(RDatadir,"overexpressed_pers6_genes_vs_UNT.RData"))
metadata <- annot_sel2
anocol = anocol2
####################################################################
### COMPARISON SETUP  
####################################################################

 mygps <- list(
'persister'=metadata[which(metadata$sample_id %in% c("HBCx95_persister_6souris") ),"cell_id"]
)
 myrefs <- list(
  'UNT'=metadata[which(metadata$sample_id %in% c("HBCx95_v3_UNT") ),"cell_id"]
 )
refs <- names(myrefs) 
groups <- names(mygps) 


##############################################################################################
###############			LOAD FILES			        			##############################
##############################################################################################
load(file.path(RDatadir_unsup,"LogCounts.RData"))

# Select only Untreated & Persister (6x mouse) cells
set.seed(2047); 

set.seed(2047);
pers_151 = sample(mygps$persister,151,replace = F)
unt_151 = sample(myrefs$UNT,151,replace = F)
unt_pers_cells_151 = c(pers_151,unt_151)
rownames(anocol) = metadata$cell_id

unt_pers_all_cells =  c(mygps$persister, myrefs$UNT)

cell_list = c("unt_pers_cells_151","unt_pers_all_cells")

# Different set of genes
setwd(maindir)
load("common_over_genes_pers_vs_unt.RData")
predictor_genes_MM468 = c("KRT16","KRT17","AKR1B1","TPM2","AL161431.1","TPM2","IGFL2-AS1")
predictor_genes_MM468 = intersect(predictor_genes_MM468,rownames(LogCounts))
diff_genes <-  c("NNMT","TAGLN","INHBA","KRT6A","KRT6B","KRT16","KRT14","KLK10","KLK5","KLK7")
diff_genes_figures <-  c("BMP6","CDKN2B","TGFBR2","TGFBR3","TGFB2","INHBA","INHBB","SMAD2","TGFB1","FOSL1","NNMT","KRT14","KLK10","KLK5","TAGLN","NNMT","ELN","TGFBR3")
gene_lists <- c("common_overexpressed_genes","overexpressed_pers_genes","predictor_genes_MM468","diff_genes","diff_genes_figures")


# TEST MULTIPLE COMBINATION : PERS GENES / PATHWAYS - UNTvsPers6 / UNTvsPersAll
for(gene_sel in gene_lists){
  
  genes = get(gene_sel)
  anocol. = anocol[unt_pers_cells_151,]
  
  mat <- LogCounts[rownames(LogCounts) %in% genes,unt_pers_cells_151]
  tmat <- t(mat)
  dim(mat)
  
  distHC <- c("distPearson","distCosine","euclidean","maximum","manhattan","canberra","binary","minkowski")[3]
  methHC <- c("ward","ward.D","ward.D2","single","complete","average")[3]
  if(distHC=="distPearson")	d <- distPearson(tmat)
  if(distHC=="distCosine")	d <- distCosine(tmat)
  if(distHC %in% c("euclidean","maximum","manhattan","canberra","binary","minkowski"))	
    d <- dist(tmat,distHC)
  hc <- hclust(d,method=methHC)
  
  mat.so <- as.matrix(mat[,hc$order])
  
  if(chRangeHM){
    for(i in 1:nrow(mat.so)){
      mat.so[i,] <- geco.changeRange(mat.so[i,],newmin=0,newmax=1)
    }
  }
  gc()
  png(file.path(resdir,paste0("Cooccurrence_pers_genes_",gene_sel,"_",distHC,"_",methHC,".png")), height=4000,width=4000,res=300)
  
  geco.hclustAnnotHeatmapPlot(x=(mat.so),
                              hc=hc,
                              hmColors=hmColors,
                              anocol=anocol.[hc$order,c(1,2),drop=F],
                              xpos=c(0.15,0.9,0.164,0.885),
                              ypos=c(0.1,0.5,0.5,0.6,0.62,0.95),
                              dendro.cex=0.01,
                              xlab.cex=1,
                              hmRowNames=TRUE,
                              hmRowNames.cex=0.5
  )
  dev.off()
}

########## ########## ########## ########## 
########## Most co-occuring in DMSO ########## 

anocol. = anocol[unt_pers_all_cells,]
genes = overexpressed_pers_genes
mat <- LogCounts[rownames(LogCounts) %in% genes,unt_pers_all_cells]
hist(colSums(mat),breaks=150)
bin_mat = mat
bin_mat[(bin_mat>0)]=1

coocurrence_persister_genes_score = as.data.frame(colSums(bin_mat) / nrow(bin_mat))
colnames(coocurrence_persister_genes_score) = "coocurrence_score"
coocurrence_persister_genes_score$sample = ""
coocurrence_persister_genes_score[mygps$persister,"sample"] = "Persister"
coocurrence_persister_genes_score[myrefs$UNT,"sample"] = "Untreated"

hist(coocurrence_persister_genes_score$coocurrence_score,breaks=150)
coocurrence_persister_genes_score$sample = factor(coocurrence_persister_genes_score$sample,levels=
                                                    c("Untreated","Persister"))

png(file.path(resdir,paste0("Coocurrence_persiter_genes_score_violin.png")),height=1000,width=1151,res=300)
ggplot(coocurrence_persister_genes_score, aes(x=sample,y=coocurrence_score,
                                              fill=sample)) + 
  geom_violin(alpha=0.8) + theme_classic() + scale_fill_manual(values=as.character(unique(anocol.[,"sample_id"])[2:1])) +
  stat_summary(fun=median, geom="point", size=2, color="black") + geom_jitter(size=0.2, alpha=0.2) + 
  ggpubr::stat_compare_means(method = "t.test",ref.group = "Untreated")
dev.off()

png(file.path(resdir,paste0("Coocurrence_persiter_genes_score_density.png")),height=1000,width=1151,res=300)
ggplot(coocurrence_persister_genes_score, aes(coocurrence_score,fill=sample,color=sample)) +
  geom_density(alpha=0.2) + theme_classic() + 
  scale_color_manual(values=as.character(unique(anocol.[,"sample_id"])[2:1])) +
  scale_fill_manual(values=as.character(unique(anocol.[,"sample_id"])[2:1])) +
  geom_vline(xintercept = quantile(coocurrence_persister_genes_score$coocurrence_score[which(coocurrence_persister_genes_score$sample=="Persister_6")],0.25), col=unique(anocol.[,"sample_id"])[1],lty=3)
dev.off()

############### ########## ########## 
########## ########## ########## ####

##########  Gene to Gene correlation: untreated vs persister 6########## 
for(i in 1:2){
  cells = list(mygps$persister, myrefs$UNT)[[i]]
  name = c("persister","untreated")[i]
  genes = overexpressed_pers_genes
  mat <- LogCounts[rownames(LogCounts) %in% genes,cells]
  if(length(which(rowSums(mat)==0))>0) mat = mat[-which(rowSums(mat)==0),]
  annot_genes = data.frame(gene= rownames(mat), total_counts = rowSums(mat))
  rownames(annot_genes) = rownames(mat)
  anocol_genes = geco.annotToCol4(annot_genes)
  cor_genes = cor(t(as.matrix(mat)),method = "pearson")
  hc_genes = hclust(as.dist(1 - cor_genes), method="ward.D")
  
  distHC <- c("distPearson","distCosine","euclidean","maximum","manhattan","canberra","binary","minkowski")[3]
  methHC <- c("ward","ward.D","ward.D2","single","complete","average")[3]
  if(distHC=="distPearson")	d <- distPearson(tmat)
  if(distHC=="distCosine")	d <- distCosine(tmat)
  if(distHC %in% c("euclidean","maximum","manhattan","canberra","binary","minkowski"))	
    d <- dist(tmat,distHC)
  hc <- hclust(d,method=methHC)
  
  mat.so <- as.matrix(cor_genes[hc_genes$order,hc_genes$order])
  if(chRangeHM){
    for(i in 1:nrow(mat.so)){
      mat.so[i,] <- geco.changeRange(mat.so[i,],newmin=0,newmax=1)
    }
  }
  png(file.path(resdir,paste0("Cooccurrence_pers_genes_heatmap_",name,".png")), height=8000,width=8000,res=600)
  geco.hclustAnnotHeatmapPlot(x=(mat.so),
                              hc=hc_genes,
                              hmColors=hmColors,
                              anocol=anocol_genes[hc_genes$order,],#[,ncol(cc.col):1]
                              xpos=c(0.15,0.9,0.164,0.885),
                              ypos=c(0.1,0.5,0.5,0.6,0.62,0.95),
                              dendro.cex=0.01,
                              xlab.cex=0.3,
                              hmRowNames=TRUE,
                              hmRowNames.cex=0.1
  )
  dev.off()
}

###########################################################
## Minimum set of features to find the sample of origin ###
###########################################################

library(doParallel)  # for parallel backend to foreach
library(foreach)     # for parallel processing with for loops
library(caret)       # for general model fitting
library(rpart)       # for fitting decision trees
library(ipred)  
library(mlbench)
library(doMC)

mat = as.data.frame(as.matrix(t(LogCounts[rownames(LogCounts) %in% 
                                                     overexpressed_pers_genes,
                                                   unt_pers_cells_151])))
dim(mat)
class = rep(1,length(unt_pers_cells_151))
class[grep("UNT",unt_pers_cells_151)] = 0


# load the data
# define the control using a random forest selection function
control <- rfeControl(functions=rfFuncs, method="cv", number=5)
# run the RFE algorithm
registerDoMC(cores = 7)
set.seed(2047)
results <- rfe(mat, class, sizes=c(2:25,seq(30,length(overexpressed_pers_genes)-1,5)), rfeControl=control)

save(results, file = file.path(RDatadir,"results_rfe_PDX.RData"))
gc()
# OR load results :
# load(file.path(RDatadir,"results_rfe_MM468.RData"))

#summarize the results
print(results)

# plot the results
png(file.path(resdir,paste0("curve_predictor_genes.png")), height=1000,1000)
plot(results, type=c("g", "o"))
dev.off()

# Set predictor genes as TOP 7 genes (as in PDX) (top 7 of subsampling = 7 variables)
predictor_genes = unique(results$variables$var[which(results$variables$Variables==10)])[1:10]

mat = LogCounts[rownames(LogCounts) %in% predictor_genes,unt_pers_cells_151]
tmat <- t(mat)
distHC <- c("distPearson","distCosine","euclidean","maximum","manhattan","canberra","binary","minkowski")[3]
methHC <- c("ward","ward.D","ward.D2","single","complete","average")[3]
if(distHC=="distPearson")	d <- distPearson(tmat)
if(distHC=="distCosine")	d <- distCosine(tmat)
if(distHC %in% c("euclidean","maximum","manhattan","canberra","binary","minkowski"))	
  d <- dist(tmat,distHC)
hc <- hclust(d,method=methHC)
mat.so <- as.matrix(mat)
cor_genes = cor(t(as.matrix(mat)),method = "pearson")
hc_genes = hclust(as.dist(1 - cor_genes), method="ward.D")
hc$labels = rep("",length(hc$labels))

if(chRangeHM){
  for(i in 1:nrow(mat.so)){
    mat.so[i,] <- geco.changeRange(mat.so[i,],newmin=0,newmax=1)
  }
}

png(file.path(resdir,paste0("Cooccurrence_pers_cells_predictor_genes.png")), height=4000,width=4000,res=600)
corColors <- colorRampPalette(c("royalblue","white","indianred1"))(256)
geco.hclustAnnotHeatmapPlot.withColumn(x=mat.so[hc_genes$order,hc$order],
                                       hc=hc,
                                       hc_row=hc_genes,
                                       hmColors=corColors,
                                       anocol=anocol.[hc$order,c(1,2),drop=F],
                                       anorow=anocol_genes[hc_genes$order,c(1,2)],
                                       xpos=c(0.375, 0.9, 0.3745, 0.885,0.05,0.325),
                                       ypos=c(0.1, 0.5, 0.5, 0.6, 0.62, 0.95),
                                       dendro.cex=0.5,
                                       xlab.cex=0.8,
                                       hmRowNames=FALSE,
                                       hmRowNames.cex=0.5
)
dev.off()

png(file.path(resdir,paste0("Cooccurrence_pers_cells_predictor_genes_hc.png")), height=800,width=2000)
plot(hc)
dev.off()

######## GINI coefficients ############
load("common_over_genes_pers_vs_unt.RData")
# house_keeping_genes = c("GAPDH","VGF","CCNA2","LMNA","GAPDH","UBB","UBC")
# From  https://www.genomics-online.com/resources/16/5049/housekeeping-genes/
# Eisenberg, Levanon: “Human housekeeping genes, revisited.

house_keeping_genes = c("GAPDH","ACTB","RRN18S","PGK1","PPIA","RPL13A","RPLP0",
                        "ARBP","B2M","YWHAZ","SDHA","TFRC","GUSB","HMBS","HPRT1","TBP")
house_keeping_genes = (intersect(house_keeping_genes,rownames(LogCounts)))

load(file.path(RDatadir_unsup,"RawCounts_features.RData"))
rownames(RawCounts) = gsub("hg19_","",rownames(RawCounts))
mat = RawCounts[which(rownames(RawCounts) %in% c(overexpressed_pers_genes,
                                                 house_keeping_genes)), intersect(myrefs$UNT,colnames(RawCounts))]

bin_mat =  (mat)
bin_mat[bin_mat>0] = 1

Gin_unt = edgeR::gini(t(bin_mat[]))

Gini_focus = Gin_unt[which(names(Gin_unt) %in% c(house_keeping_genes,common_overexpressed_genes))]
Gini_focus = data.frame(Gini_focus,0)
colorlist = ifelse(rownames(Gini_focus) %in% house_keeping_genes, "grey", "red")
pdf(file.path(resdir,"GiniScores_common_pers_genes_hk_inInitial.pdf"),width = 12)
plot(0,0,type="n",xlim=c(-0.06,1.06), ylim=c(0,2), yaxt = 'n', main ="Gini Scores")
for (i in 1:length(Gini_focus$Gini_focus)) { stripchart(Gini_focus$Gini_focus[i],
                                                        add = T, bg = colorlist[i],
                                                        vertical = F, pch =21, method = "jitter") }
text(Gini_focus$Gini_focus,
     (runif(length(Gini_focus$Gini_focus))/2) + rep(c(0,1.5),length(Gini_focus$Gini_focus)/2), 
     cex = 0.75,
     labels = rownames(Gini_focus), srt = 75, col = colorlist)
dev.off()

mat_pers = LogCounts[rownames(LogCounts) %in% c(overexpressed_pers_genes,
                                                                  house_keeping_genes), mygps$persister]

bin_mat_pers =  as.matrix(mat_pers)
bin_mat_pers[bin_mat_pers>0] = 1
Gin_pers = edgeR::gini(t(bin_mat_pers))

Gin_pers_focus = Gin_pers[which(names(Gin_pers) %in% c(house_keeping_genes,common_overexpressed_genes))]
Gin_pers_focus = data.frame(Gin_pers_focus,0)
colorlist = ifelse(rownames(Gin_pers_focus) %in% house_keeping_genes, "grey", "red")
pdf(file.path(resdir,"GiniScores_common_pers_genes_hk_inPersister.pdf"),width = 12)
plot(0,0,type="n",xlim=c(-0.06,1.06), ylim=c(0,2), yaxt = 'n', main ="Gini Scores")
for (i in 1:length(Gin_pers_focus$Gin_pers_focus)) { stripchart(Gin_pers_focus$Gin_pers_focus[i],
                                                                add = T, bg = colorlist[i],
                                                                vertical = F, pch =21, method = "jitter") }
text(Gin_pers_focus$Gin_pers_focus,
     (runif(length(Gin_pers_focus$Gin_pers_focus))/2) + rep(c(0,1.5),
                                                            length(Gin_pers_focus$Gin_pers_focus)/2), 
     cex = 0.75,
     labels = rownames(Gin_pers_focus), srt = 75, col = colorlist)
dev.off()

# Which genes "gained the most homogeneity" between initial and persister state ?
genes_expressed_both = intersect(rownames(Gini_focus),rownames(Gin_pers_focus))
difference_in_homogeneity = cbind(Gini_focus[which(rownames(Gini_focus) %in% genes_expressed_both),1,drop=F],
  Gin_pers_focus[which( rownames(Gin_pers_focus) %in% genes_expressed_both),1,drop=F])
difference_in_homogeneity$gain_homogeneity = difference_in_homogeneity$Gini_focus - difference_in_homogeneity$Gin_pers_focus
difference_in_homogeneity$Gene = rownames(difference_in_homogeneity)
png(file.path(resdir,paste0("Top6_Gain_in_homogeneity.png")), height=1500,width=1500,res=300)
difference_in_homogeneity %>% dplyr::arrange(dplyr::desc(gain_homogeneity)) %>%  head %>%
  dplyr::mutate(Gene = factor(Gene,levels=Gene)) %>% ggplot() +
  geom_bar(aes(x=Gene,y=gain_homogeneity),stat="identity") +
  theme_classic() + xlab("") + ylab("Difference in Gini score")
dev.off()

bin_mat =  mat
bin_mat[bin_mat>0] = 1

# high_gini_genes = names(Gin_unt)[which(Gin_unt>0.75)]
odd_ratio_mat = matrix(0,nrow = nrow(bin_mat), ncol=nrow(bin_mat)
                       ,dimnames = list(rownames(bin_mat),rownames(bin_mat)))
fisher_pvalue_mat = matrix(0,nrow = nrow(bin_mat), ncol=nrow(bin_mat)
                           ,dimnames = list(rownames(bin_mat),rownames(bin_mat)))

for(gene_k in rownames(bin_mat)){
  for(gene_i in rownames(bin_mat)){
    if(gene_k != gene_i){
      cells_i = which(bin_mat[gene_i,] > 0)
      cells_k = which(bin_mat[gene_k,] > 0)
      
      odd_ratio_mat[gene_i,gene_k] =  
        (length(intersect(cells_i,cells_k))/length(cells_i)) * (ncol(bin_mat)/length(cells_k))
      contingency_table = data.frame(a= c(length(intersect(cells_i,cells_k)),
                                          length(cells_i) - length(intersect(cells_i,cells_k))),
                                     b = c(length(cells_k),ncol(bin_mat)-length(cells_k))
      )
      fisher_pvalue_mat[gene_i,gene_k] = fisher.test(contingency_table,alternative = "greater")$p.value
    }
  }
}

odd_ratio_mat = log2(odd_ratio_mat)
odd_ratio_mat[is.infinite(odd_ratio_mat)] = 0
odd_ratio_mat[is.nan(odd_ratio_mat)] = 0

jackpot_genes = intersect(common_overexpressed_genes, rownames(odd_ratio_mat))

fisher_pvalue_mat_all = fisher_pvalue_mat
fisher_pvalue_mat = fisher_pvalue_mat[jackpot_genes,jackpot_genes]

pvalues = fisher_pvalue_mat[upper.tri(fisher_pvalue_mat)]
adjusted_pvalues = p.adjust(pvalues,"BH")
fisher_qvalue_mat = fisher_pvalue_mat
fisher_qvalue_mat[upper.tri(fisher_qvalue_mat)] = adjusted_pvalues

fisher_pvalue_mat = -log10(fisher_pvalue_mat)
fisher_pvalue_mat[fisher_pvalue_mat<0]=0
fisher_pvalue_mat[lower.tri(fisher_pvalue_mat)] = 0
fisher_pvalue_mat[is.infinite(fisher_pvalue_mat)] = 0

fisher_qvalue_mat = -log10(fisher_qvalue_mat)
fisher_qvalue_mat[fisher_qvalue_mat<0]=0
fisher_qvalue_mat[lower.tri(fisher_qvalue_mat)] = 0
fisher_qvalue_mat[is.infinite(fisher_qvalue_mat)] = 0

odd_ratio_mat_all = odd_ratio_mat
odd_ratio_mat = odd_ratio_mat[jackpot_genes,jackpot_genes]
odd_ratio_mat[odd_ratio_mat<0]=0
odd_ratio_mat[lower.tri(odd_ratio_mat)] = 0

png(file.path(resdir,paste0("LogOddRatio_rarity.png")), height=4000,width=4000,res=400)
gplots::heatmap.2((odd_ratio_mat),Rowv = F,Colv = F
                  , trace="none", density = "none", scale ="none",
                  col =rev(inferno(20)))
dev.off()

png(file.path(resdir,paste0("Pvalue_rarity.png")), height=4000,width=4000,res=400)
gplots::heatmap.2((fisher_pvalue_mat),Rowv = F,Colv = F
                  , trace="none", density = "none", scale ="none",
                  col =rev(inferno(20)))
dev.off()

png(file.path(resdir,paste0("Qvalue_rarity.png")), height=4000,width=4000,res=400)
gplots::heatmap.2((fisher_qvalue_mat),Rowv = F,Colv = F
                  , trace="none", density = "none", scale ="none",
                  col =rev(inferno(20)))
dev.off()

fisher_pvalue_mat <- geco.changeRange(fisher_pvalue_mat,newmin=0,newmax=1)
odd_ratio_mat <- geco.changeRange(odd_ratio_mat,newmin=0,newmax=1)

combined_mat = fisher_pvalue_mat + odd_ratio_mat

png(file.path(resdir,paste0("Combined_pvalue_logOddRatio_rarity.png")), height=4000,width=4000,res=400)
gplots::heatmap.2((combined_mat),Rowv = F,Colv = F
                  , trace="none", density = "none", scale ="none",
                  col =rev(inferno(20)))
dev.off()


library(RColorBrewer)
p = list()
for(gene in jackpot_genes){
  tab = as.data.frame(mat[gene,])
  colnames(tab)="Gene"
  q = quantile(tab$Gene,0.98)
  tab$Expression = ifelse(tab$Gene>=q,"high","low")
  p[[gene]] = ggplot(tab,aes(Gene)) + geom_histogram() + geom_rug(aes(color=Expression),sides="b") +
    theme_classic() + ggtitle(gene) + theme(legend.position = "None")
}

pdf(file.path(resdir,"histogram_transcript_abundance_marker_genes_in_untreated.pdf"))
gridExtra::grid.arrange(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]],p[[6]],p[[7]],p[[8]],p[[9]],nrow = 3)
gridExtra::grid.arrange(p[[10]],p[[11]],p[[12]],p[[13]],p[[14]],p[[15]],p[[16]],nrow = 3)
dev.off()

top =sort(fisher_pvalue_mat[],decreasing = T)[1:5]
for(gene_k in rownames(fisher_pvalue_mat)){
  for(gene_i in rownames(fisher_pvalue_mat)){
    if(fisher_pvalue_mat[gene_i,gene_k] %in% top){
      png(file.path(resdir,paste0(gene_i,"_",gene_k,"_dotplot.png")), height=1500,width=1500,res=300)
      tab = as.data.frame(cbind(mat[gene_i,],mat[gene_k,]))
      tab$cell_id = rownames(tab)
      colnames(tab)[1:2] = c("A","B")
      tab$total_counts = annot_sel2$total_counts[which(annot_sel2$cell_id %in% tab$cell_id)]
      print(ggplot(tab) + geom_jitter(aes(x=A,y=B,color=total_counts),alpha=0.35) + 
              theme_classic() + scale_color_viridis() +
              geom_vline(xintercept = 1, col = "red",lty=3) +
              geom_hline(yintercept = 1, col = "red",lty=3) +
              xlab(gene_i) + ylab(gene_k)
      )
      dev.off()
    }
  }
}

library(eulerr)
if(!dir.exists(file.path(resdir,"Venns"))) dir.create(file.path(resdir,"Venns"))
for(gene_k in jackpot_genes){
  p=list()
  for(gene_i in jackpot_genes){
    if(gene_k != gene_i){
      cells_i = which(bin_mat[gene_i,] > 0)
      cells_k = which(bin_mat[gene_k,] > 0)
      A = 0
      B = 0
      A_B = length(intersect(cells_i,cells_k))
      C = ncol(bin_mat)
      B_C = length(cells_k)
      A_B_C = A_B
      A_C = length(cells_i)
      venn_c  = c(A=A, B=B, C=C, "A&B"=A_B, "A&C"=A_C, "B&C"=B_C, "A&B&C"=A_B_C)
      names(venn_c) = c(gene_i,gene_k,"All",
                        paste0(gene_i,"&",gene_k),
                        paste0(gene_i,"&All"),
                        paste0(gene_k,"&All"),
                        paste0(gene_i,"&",gene_k,"&All")
      )
      e = euler(venn_c)
      p[[paste0(gene_i,"_",gene_k)]] = e
    }
  }
  pdf(file.path(resdir,"Venns",paste0("Venns_",gene_k,".pdf")))
  gridExtra::grid.arrange(plot(p[[1]],fills = c("#D6BA1C", "#C20606", "#DBDBDB"), labels = list(cex = c(0.5,0.5,1.5))),plot(p[[2]],fills = c("#D6BA1C", "#C20606", "#DBDBDB"), labels = list(cex = c(0.5,0.5,1.5))),plot(p[[3]],fills = c("#D6BA1C", "#C20606", "#DBDBDB"), labels = list(cex = c(0.5,0.5,1.5))),
                          plot(p[[4]],fills = c("#D6BA1C", "#C20606", "#DBDBDB"), labels = list(cex = c(0.5,0.5,1.5))),plot(p[[5]],fills = c("#D6BA1C", "#C20606", "#DBDBDB"), labels = list(cex = c(0.5,0.5,1.5))),plot(p[[6]],fills = c("#D6BA1C", "#C20606", "#DBDBDB"), labels = list(cex = c(0.5,0.5,1.5))),
                          plot(p[[7]],fills = c("#D6BA1C", "#C20606", "#DBDBDB"), labels = list(cex = c(0.5,0.5,1.5))),plot(p[[8]],fills = c("#D6BA1C", "#C20606", "#DBDBDB"), labels = list(cex = c(0.5,0.5,1.5))),plot(p[[9]],fills = c("#D6BA1C", "#C20606", "#DBDBDB"), labels = list(cex = c(0.5,0.5,1.5))),nrow = 3)
  gridExtra::grid.arrange(plot(p[[10]],fills = c("#D6BA1C", "#C20606", "#DBDBDB"), labels = list(cex = c(0.5,0.5,1.5))),plot(p[[11]],fills = c("#D6BA1C", "#C20606", "#DBDBDB"), labels = list(cex = c(0.5,0.5,1.5))),plot(p[[12]],fills = c("#D6BA1C", "#C20606", "#DBDBDB"), labels = list(cex = c(0.5,0.5,1.5))),
                          plot(p[[13]],fills = c("#D6BA1C", "#C20606", "#DBDBDB"), labels = list(cex = c(0.5,0.5,1.5))),plot(p[[14]],fills = c("#D6BA1C", "#C20606", "#DBDBDB"), labels = list(cex = c(0.5,0.5,1.5))),plot(p[[15]],fills = c("#D6BA1C", "#C20606", "#DBDBDB"), labels = list(cex = c(0.5,0.5,1.5))),
                          nrow = 3)
  
  dev.off()
}

png(file.path(resdir,"Venns",paste0("Venn_KLK10_KLK5_KRT14.png")), height=1500,width=1500,res=300)
v = data.frame("KLK10" =(bin_mat["KLK10",] > 0),
               "KLK5"= (bin_mat["KLK5",] > 0),
               "KRT14" = (bin_mat["KRT14",] > 0),
               "All" = T)
e = euler(v,shape = "ellipse")
plot(e, quantities = TRUE,labels=list(cex=),fills = c("#4b164bff", "#fbefafff","#f66f5eec","#DBDBDB"))
dev.off()


fisher_pvalue_mat = fisher_pvalue_mat_all[house_keeping_genes,house_keeping_genes]

pvalues = fisher_pvalue_mat[upper.tri(fisher_pvalue_mat)]
adjusted_pvalues = p.adjust(pvalues,"BH")
fisher_qvalue_mat = fisher_pvalue_mat
fisher_qvalue_mat[upper.tri(fisher_qvalue_mat)] = adjusted_pvalues

fisher_pvalue_mat = -log10(fisher_pvalue_mat)
fisher_pvalue_mat[fisher_pvalue_mat<0]=0
fisher_pvalue_mat[lower.tri(fisher_pvalue_mat)] = 0
fisher_pvalue_mat[is.infinite(fisher_pvalue_mat)] = 0

fisher_qvalue_mat = -log10(fisher_qvalue_mat)
fisher_qvalue_mat[fisher_qvalue_mat<0]=0
fisher_qvalue_mat[lower.tri(fisher_qvalue_mat)] = 0
fisher_qvalue_mat[is.infinite(fisher_qvalue_mat)] = 0

odd_ratio_mat = odd_ratio_mat_all[house_keeping_genes,house_keeping_genes]
odd_ratio_mat[odd_ratio_mat<0]=0
odd_ratio_mat[lower.tri(odd_ratio_mat)] = 0

png(file.path(resdir,paste0("LogOddRatio_rarity_HK.png")), height=4000,width=4000,res=400)
gplots::heatmap.2((odd_ratio_mat),Rowv = F,Colv = F
                  , trace="none", density = "none", scale ="none",
                  col =rev(inferno(20)))
dev.off()

png(file.path(resdir,paste0("Pvalue_rarity_HK.png")), height=4000,width=4000,res=400)
gplots::heatmap.2((fisher_pvalue_mat),Rowv = F,Colv = F
                  , trace="none", density = "none", scale ="none",
                  col =rev(inferno(20)))
dev.off()

png(file.path(resdir,paste0("Qvalue_rarity_HK.png")), height=4000,width=4000,res=400)
gplots::heatmap.2((fisher_qvalue_mat),Rowv = F,Colv = F
                  , trace="none", density = "none", scale ="none",
                  col =rev(inferno(20)))
dev.off()

library(RColorBrewer)
p = list()
for(gene in house_keeping_genes){
  tab = as.data.frame(mat[gene,])
  colnames(tab)="Gene"
  q = quantile(tab$Gene,0.98)
  tab$Expression = ifelse(tab$Gene>=q,"high","low")
  p[[gene]] = ggplot(tab,aes(Gene)) + geom_histogram() + geom_rug(aes(color=Expression),sides="b") +
    theme_classic() + ggtitle(gene) + theme(legend.position = "None")
}

pdf(file.path(resdir,"histogram_transcript_abundance_hk.pdf"))
gridExtra::grid.arrange(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]],p[[6]],p[[7]],p[[8]],p[[9]],nrow = 3)
gridExtra::grid.arrange(p[[10]],p[[11]],p[[12]],p[[13]],nrow = 3)
dev.off()



######## CO-OCCURRENCE Persister Genes in DMSO / Persisters -UMAP #############
load(file.path(RDatadir_unsup,"umap.RData"))
umap_res. = umap_res[rownames(coocurrence_persister_genes_score),]
anocol. = geco.annotToCol4(coocurrence_persister_genes_score)
png(file.path(resdir,paste0("UMAP_",colnames(anocol.)[1],".png")), height=1350,width=1200,res=300)
plot((umap_res.), col=alpha(anocol.[,1],0.3),pch=20,cex=0.6,
     main=paste0(colnames(anocol.)[1]," min=",round(min(coocurrence_persister_genes_score[,colnames(anocol.)[1]]),digits=3),
                 " max=",round(max(coocurrence_persister_genes_score[,colnames(anocol.)[1]]),digits=3)),
     xlab="component 1",ylab="component 2")
dev.off()


######## VARIABILITY OF PERSISTER GENES #############
mat = LogCounts[rownames(LogCounts) %in% overexpressed_pers_genes, unt_pers_cells_151]
rownames(mat) = overexpressed_pers_genes
bin_mat =  as.matrix(mat)
bin_mat[bin_mat>0] = 1
dim(bin_mat)
vars_p = apply(bin_mat[,1:151], MARGIN = 1, var)
vars_u = apply(bin_mat[,152:302], MARGIN = 1, var)

vars = data.frame("gene" = overexpressed_pers_genes ,"Untreated" = c(vars_u),
                  "Persister" = vars_p)

vars = vars %>% tidyr::gather("sample","Gene Variance",-gene)
vars = vars %>% dplyr::mutate(sample = factor(sample,levels=c("Untreated","Persister")))
vars[1:length(overexpressed_pers_genes),"Cells 'ON' (%)"] = rowSums(bin_mat[,152:302]) /5
vars[(length(overexpressed_pers_genes)+1):(2*length(overexpressed_pers_genes)),"Cells 'ON' (%)"] = rowSums(bin_mat[,1:151])/5
png(file.path(resdir,paste0("Variability_untreated_persister_rawcounts_min1.png")), height=1151,width=2000,res=300)
vars %>% ggplot(aes(x=sample,y=`Gene Variance`)) + geom_violin(alpha=0.75,aes(fill=sample)) +
  geom_jitter(alpha = 0.5,aes(color=`Cells 'ON' (%)`),width = 0.15) + theme_classic() +
  scale_fill_manual(values=as.character(unique(anocol[unt_pers_cells_151,"sample_id"])[2:1])) +
  scale_color_continuous(type="viridis")
dev.off()


vars %>% dplyr::filter((sample == "Untreated" & `Gene Variance` > 0.2 & `Cells 'ON' (%)` > 50) | 
                         (sample == "Persister" & `Gene Variance` < 0.05  & `Cells 'ON' (%)` > 50) )

# Cell to Cell variability of persister genes
mat = LogCounts[rownames(LogCounts) %in% overexpressed_pers_genes, c(mygps$persister,myrefs$UNT)]
bin_mat =  as.matrix(mat)
bin_mat[bin_mat>0] = 1
vars = apply(bin_mat, MARGIN = 2, var)
coocurrence_persister_genes_score$homogeneity = 1-vars
anocol. = geco.annotToCol4(coocurrence_persister_genes_score,plotLegend=T,
                           plotLegendFile=file.path(resdir,"Annotation_legends.pdf"))
png(file.path(resdir,paste0("UMAP_homogeneity.png")), height=1350,width=1200,res=300)
plot((umap_res.), col=alpha(anocol.[,"homogeneity"],0.3),pch=20,cex=0.6,
     main=paste0("homogeneity"," min=",round(min(coocurrence_persister_genes_score[,"homogeneity"]),digits=3),
                 " max=",round(max(coocurrence_persister_genes_score[,"homogeneity"]),digits=3)),
     xlab="component 1",ylab="component 2")
dev.off()
