######## GINI coefficients PDX - scChIP ############

library(here)
library(Matrix)
library(readxl)
library(ggplot2)
library(rtracklayer)
library(GenomicRanges)
library(geco.visu)
library(geco.supervised)
library(geco.unsupervised)
library(ape)
library(viridis)
source(file.path(here(),"Scripts","functions.R"))
source(file.path(here(),"Scripts","global_var.R"))

# Parameters
percent_coocurrence_to_define_top_cells = 0.3

maindir = here()
outdir = file.path(maindir,"output","scChIPseq","PDX","Coocurrence_K4")
if(!dir.exists(outdir)) dir.create(outdir)

# Find association Gene to Peak by using PDX K4 peaks "consensus annotation" obtained from bulk ChIPseq 
annot_PDX_K4 <- rtracklayer::import(file.path(maindir,"annotation","HBCx95_peaks_K4.bed.gz"))
load(file.path(maindir,"annotation","hg38.GeneTSS.rda"))
hg38.GeneTSS = as(hg38.GeneTSS,"GRanges")
annot10k <- read.table(unz(file.path(maindir,"annotation","gencode.v34.annotation.transcriptTSS_10k.bed.zip"),
                           "gencode.v34.annotation.transcriptTSS_10k.bed"))
colnames(annot10k) <- c("chr","start","end","transcripts","gene","sense")
annot10k = as(annot10k[,c(1:3,5)],"GRanges")
annot10k$ID = paste0(seqnames(annot10k),":",start(annot10k),"-",end(annot10k))
  
differential_results_scRNA <- readxl::read_xlsx(file.path(
  maindir,"output","scRNAseq","PDX","Supervised","Tables",
  "Differential_analysis_Limma_logFC_2.xlsx"),
  sheet=1)

overexpressed_PDX_persister_genes = differential_results_scRNA$Symbol
transcripts = annot10k[which(annot10k$gene %in% overexpressed_PDX_persister_genes)]
regions = hg38.GeneTSS[which(hg38.GeneTSS$gene %in% overexpressed_PDX_persister_genes)]
regions = transcripts[subjectHits(findOverlaps(regions,transcripts)),]
overexpressed_PDX_persister_regions = regions$ID

# house_keeping_genes = c("GAPDH","VGF","CCNA2","LMNA","GAPDH","UBB","UBC")
# From  https://www.genomics-online.com/resources/16/5049/housekeeping-genes/
# Eisenberg, Levanon: “Human housekeeping genes, revisited.

house_keeping_genes = c("TFRC", "SDHA", "TBP", "ACTB", "PPIA", "GUSB", "YWHAZ",
                        "HMBS", "GAPDH", "RPLP0", "B2M", "RPL13A", "PGK1","HPRT1")
house_keeping_genes_row <- c(11477, 13853, 19080, 19230, 19841, 20044, 23094,
                             30210, 30751, 32985, 36858, 46622, 51202, 51630) #manually annotated with the right TSS

annot10k_housekeeping <- annot10k[house_keeping_genes_row,]
house_keeping_genes_region = annot10k_housekeeping$ID

bivalent_genes <- c("BMP6","LAMB1", "VIM", "FOSL1", "ABCC4",
                    "COL4A2", "KRT14", "TGFB1", "KLK10" )
bivalent_genes_row <- c(16677, 20736, 25943, 29195, 34181, 34369, 41923, 46164, 46727)
annot10k_bivalent <- annot10k[bivalent_genes_row,]

path_matrices = file.path(  maindir,"input","scChIPseq",
                            "PDX","Count_Matrices")

HBCx95_UNT = scater::readSparseCounts(
  unzip(file.path(path_matrices,"HBCx95_m43_UNT_K4.tsv.zip"),
        files="HBCx95_m43_UNT_K4.tsv"))

#Determine the cells with the top cells with the most efficient IP (the most presence of H3K4me3 fragments for housekeeping genes)
mat <- HBCx95_UNT
mat_housekeeping <- mat[house_keeping_genes_row,]
bin_mat = mat_housekeeping
bin_mat[(bin_mat>0)]=1
coocurrence_K4me3_score_housekeeping = as.data.frame(Matrix::colSums(bin_mat) / nrow(bin_mat))
colnames(coocurrence_K4me3_score_housekeeping) = "Coocurrence_score"
top_cells_init <- rownames(coocurrence_K4me3_score_housekeeping)[which(
  coocurrence_K4me3_score_housekeeping$Coocurrence_score>percent_coocurrence_to_define_top_cells)]

mat_top  <- mat[,top_cells_init]
Gin_unt = edgeR::gini(t(mat_top))
Gini_focus = Gin_unt[which(names(Gin_unt) %in% c(annot10k_housekeeping$ID,annot10k_bivalent$ID))]
Gini_focus = data.frame(Gini_focus,0)
colorlist = ifelse(rownames(Gini_focus) %in% annot10k_housekeeping$ID, "grey", "red")
Gini_focus_DMSOi <- Gini_focus

mat_bivalent <- mat_top[bivalent_genes_row,]
mat_bivalent_DMSOi <- mat_bivalent

bin_mat = mat_bivalent
bin_mat[(bin_mat>0)]=1
coocurrence_K4me3_score_bivalent_DMSOi = as.data.frame(Matrix::colSums(bin_mat) / nrow(bin_mat))
colnames(coocurrence_K4me3_score_bivalent_DMSOi) = "Coocurrence_score"

pdf(file.path(outdir,"GiniScores_H3K4me3_bivalent_genes_inK4_DMSOi.pdf"),width = 12)
plot(0,0,type="n",xlim=c(-0.06,1.06), ylim=c(0,2), yaxt = 'n', main ="Gini Scores")
for (i in 1:length(Gini_focus$Gini_focus)) { stripchart(Gini_focus$Gini_focus[i],
                                                        add = T, bg = colorlist[i],
                                                        vertical = F, pch =21, method = "jitter") }
text(Gini_focus$Gini_focus,
     (runif(length(Gini_focus$Gini_focus))/2) + rep(c(0,1.5),length(Gini_focus$Gini_focus)/2), 
     cex = 0.75,
     labels = annot10k$gene[annot10k$ID %in% rownames(Gini_focus)], srt = 75, col = colorlist)
dev.off()

# Violin plot co-occurence score
mat_init = HBCx95_UNT

png(file.path(outdir,paste0("Coocurrence_score_overpers_K4_violin.png")),height=1000,width=1500,res=300)
violin_plot_cooccurence(mat_init,NULL,
                        sub_regions =  overexpressed_PDX_persister_regions,
                        sub_cells = list(top_cells_init),levels = "Untreated",colors = "#dfdfdfff")
dev.off()

png(file.path(outdir,paste0("Coocurrence_score_hk_K4_violin.png")),height=1000,width=1500,res=300)
violin_plot_cooccurence(mat_init,NULL,
                        sub_regions =  house_keeping_genes_region,
                        sub_cells = list(top_cells_init),levels = "Untreated",colors = "#dfdfdfff")
dev.off()

#Heatmap on top cells of DMSOi and 5FU6
tmat = Matrix::cBind(mat_pers,mat_init)
tmat = tmat[match(c(overexpressed_PDX_persister_regions,house_keeping_genes_region),rownames(tmat)),
            match(c(top_cells_init),colnames(tmat))]
tmat = tmat[-which(duplicated(rownames(tmat))),]
bin_mat = tmat
bin_mat[bin_mat>0] = 1

methHC <- c("ward","ward.D","ward.D2","single","complete","average")[3]
hc <- hclust(as.dist(1-cor(as.matrix(bin_mat))),method=methHC)
# chRangeHM = TRUE
mat.so <- as.matrix(bin_mat)
# if(chRangeHM){
#   for(i in 1:nrow(mat.so)){
#     mat.so[i,] <- geco.changeRange(mat.so[i,],newmin=0,newmax=1)
#   }
# }

#Create cell annotation data.frame
annot_cell = data.frame(cell_id = c(top_cells_init),
                        sample_id = c(rep("Initital",length(top_cells_init))))
rownames(annot_cell) = annot_cell$cell_id
anocol = geco.annotToCol4(annot_cell)

hc_gene = hclust(dist(as.matrix(bin_mat)),method=methHC)
hc_gene$labels = rep("",length(hc_gene$labels))
#Create gene/region annotation data.frame
annot_gene = data.frame(region = c(overexpressed_PDX_persister_regions,house_keeping_genes_region),
                        genes = c(as.character(regions$gene),as.character(annot10k_housekeeping$gene)),
                        type = c(rep("Persister",length(overexpressed_PDX_persister_regions)),
                                   rep("Housekeeping",length(house_keeping_genes_region))))
annot_gene = annot_gene[-which(duplicated(annot_gene$region)),]
rownames(annot_gene) = annot_gene$region
anocol_gene = geco.annotToCol4(annot_gene)

png(file.path(outdir,paste0("Cooccurrence_pers_hk_topcells_K4_heatmap.png")),
    height=8000,width=8000,res=600)
geco.hclustAnnotHeatmapPlot.withColumn(
  x=mat.so[hc_gene$order,hc$order],
  hc=hc,
  hmColors=hmColors,
  anocol=anocol[hc$order,],
  hc_row = hc_gene,
  anorow = anocol_gene[hc_gene$order,],
  dendro.cex=0.01,
  xlab.cex=1,
  hmRowNames=FALSE,
  hmRowNames.cex=0.5)
dev.off()


###############################################################
#K27me3 study
###############################################################

#Use peak annot to select top cells and then try to switch to 10k
# samples = c("PDX_5FU6_day33","PDX_DMSO6_day60")
# for(file in samples){
#   
# }
# M5FU6_day33 <- read.table("PDX_K27_peaks_1000/PDX_5FU6_day33.tsv")
# DMSO6_day60 <- read.table("PDX_K27_peaks_1000/PDX_DMSO6_day60.tsv")
# 
# 
# loci_control <- c("chr7_127610001_127640000","chr7_127610001_127640000","chr19_34686001_34730000","chr19_13276001_13725000","chr19_15609001_16005000","chr19_50514001_50788000","chr20_16764001_17491000")
# gene_control <- c("CACNA1A","PAX4","ZNF302","CYP4F2","SHANK1","PCSK2")
# 
# mat <- DMSO6_day60
# mat_loci_control <- mat[loci_control,]
# bin_mat = mat_loci_control
# bin_mat[(bin_mat<2)]=0
# bin_mat[(bin_mat>1)]=1
# 
# coocurrence_K4me3_score_housekeeping = as.data.frame(colSums(bin_mat) / nrow(bin_mat))
# top_cells <- rownames(coocurrence_K4me3_score_housekeeping)[which(coocurrence_K4me3_score_housekeeping$`colSums(bin_mat)/nrow(bin_mat)`>percent_coocurrence_to_define_top_cells)]
# 
# 
# annotZerone$Gini_K4_DMSO6_day60 <- Gin_unt
# Gini_focus = Gin_unt[which(names(Gin_unt) %in% c(house_keeping_genes,overexpressed_pers_genes))]
# Gini_focus = data.frame(Gini_focus,0)
# 
# colorlist = ifelse(rownames(Gini_focus) %in% house_keeping_genes, "grey", "red")
# pdf(file.path(outdir,"GiniScores_common_pers_genes_hk_inInitial.pdf"),width = 12)
# plot(0,0,type="n",xlim=c(-0.06,1.06), ylim=c(0,2), yaxt = 'n', main ="Gini Scores")
# for (i in 1:length(Gini_focus$Gini_focus)) { stripchart(Gini_focus$Gini_focus[i],
#                                                         add = T, bg = colorlist[i],
#                                                         vertical = F, pch =21, method = "jitter") }
# text(Gini_focus$Gini_focus,
#      (runif(length(Gini_focus$Gini_focus))/2) + rep(c(0,1.5),length(Gini_focus$Gini_focus)/2), 
#      cex = 0.75,
#      labels = rownames(Gini_focus), srt = 75, col = colorlist)
# dev.off()
