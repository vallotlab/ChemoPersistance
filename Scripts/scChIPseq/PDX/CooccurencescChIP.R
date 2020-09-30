######## GINI coefficients MM468 - scChIP ############

library(here)
library(Matrix)

maindir = here()
outdir = file.path(maindir,"output","scChIPseq","MM468","Coocurrence_K4")
if(!dir.exists(outdir)) dir.create(outdir)

load(file.path(maindir,"output","scRNAseq","common_over_genes_pers_vs_unt.RData"))

annot10k <- read.table(unz(file.path(maindir,"annotation","gencode.v34.annotation.transcriptTSS_10k.bed.zip"),
                           "gencode.v34.annotation.transcriptTSS_10k.bed"))

colnames(annot10k) <- c("chr","start","end","transcripts","gene","sense")
annot10k$ID <- paste(paste(annot10k$chr,annot10k$start,sep=":"),annot10k$end,sep="-")

# house_keeping_genes = c("GAPDH","VGF","CCNA2","LMNA","GAPDH","UBB","UBC")
# From  https://www.genomics-online.com/resources/16/5049/housekeeping-genes/
# Eisenberg, Levanon: “Human housekeeping genes, revisited.

house_keeping_genes = c("GAPDH","ACTB","RRN18S","PGK1","PPIA","RPL13A","RPLP0",
                        "ARBP","B2M","YWHAZ","SDHA","TFRC","GUSB","HMBS","HPRT1","TBP")
house_keeping_genes_row <- c(11745,52102,60829,52851,30219,14438,56781,45738,42778,53080,11123,51916) #manually annotated with the right TSS
annot10k_housekeeping <- annot10k[house_keeping_genes_row,]

bivalent_genes <- c("TGFB1","BMP6","KLK10","ABCC4","FOSL1","VIM","COL4A2","LAMB1","KRT14")
bivalent_genes_row <- c(29719,49126,30317,15839,9916,6113,16055,53950,24851)
annot10k_bivalent <- annot10k[bivalent_genes_row,]

file_names = unzip(file.path(
  maindir,"input","scChIPseq","MM468","MM468_K4_transcripts_10k_1000.zip"),
  list = T)$Name
list_files = vector("list", length(file_names))
names(list_files) = gsub(".tsv","",file_names)
path_matrices = file.path(  maindir,"input","scChIPseq",
                            "MM468","MM468_K4_transcripts_10k_1000.zip")

#10k counts
for(file in names(list_files)) {
  list_files[[file]] = scater::readSparseCounts(
    unzip(zipfile = path_matrices,files = paste0(file,".tsv")))
}

#Determine the cells with the top cells with the most efficient IP (the most presence of H3K4me3 fragments for housekeeping genes)
mat <- list_files$MM468_Rpop_DMSO1_K4me3
mat_housekeeping <- mat[house_keeping_genes_row,]
bin_mat = mat_housekeeping
bin_mat[(bin_mat>0)]=1
coocurrence_K4me3_score_housekeeping = as.data.frame(Matrix::colSums(bin_mat) / nrow(bin_mat))
colnames(coocurrence_K4me3_score_housekeeping) = "Coocurrence_score"
top_cells <- rownames(coocurrence_K4me3_score_housekeeping)[which(coocurrence_K4me3_score_housekeeping$Coocurrence_score>0.4)]

mat_top  <- mat[,top_cells]
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

#Determine the cells with the top cells with the most efficient IP (the most presence of H3K4me3 fragments for housekeeping genes)
mat <- list_files$MM468_5FU6_D60_K4
mat_housekeeping <- mat[house_keeping_genes_row,]
bin_mat = mat_housekeeping
bin_mat[(bin_mat>0)]=1
coocurrence_K4me3_score_housekeeping = as.data.frame(colSums(bin_mat) / nrow(bin_mat))
colnames(coocurrence_K4me3_score_housekeeping) = "Coocurrence_score"

top_cells <- rownames(coocurrence_K4me3_score_housekeeping)[which(coocurrence_K4me3_score_housekeeping$Coocurrence_score>0.4)]

mat_top  <- mat[,top_cells]
#bin_mat[bin_mat>0] = 1
Gin_unt = edgeR::gini(t(mat_top))
Gini_focus = Gin_unt[which(names(Gin_unt) %in% c(annot10k_housekeeping$ID,annot10k_bivalent$ID))]
Gini_focus = data.frame(Gini_focus,0)
Gini_focus_5FU6 <- Gini_focus
colorlist = ifelse(rownames(Gini_focus) %in% annot10k_housekeeping$ID, "grey", "red")

mat_bivalent <- mat_top[bivalent_genes_row,]
mat_bivalent_5FU6 <- mat_bivalent
bin_mat = mat_bivalent
bin_mat[(bin_mat>0)]=1
coocurrence_K4me3_score_bivalent_5FU6 = as.data.frame(colSums(bin_mat) / nrow(bin_mat))
colnames(coocurrence_K4me3_score_bivalent_5FU6) = "Coocurrence_score"


pdf(file.path(outdir,"GiniScores_H3K4me3_bivalent_genes_in5FU6.pdf"),width = 12)
plot(0,0,type="n",xlim=c(-0.06,1.06), ylim=c(0,2), yaxt = 'n', main ="Gini Scores")
for (i in 1:length(Gini_focus$Gini_focus)) { stripchart(Gini_focus$Gini_focus[i],
                                                        add = T, bg = colorlist[i],
                                                        vertical = F, pch =21, method = "jitter") }
text(Gini_focus$Gini_focus,
     (runif(length(Gini_focus$Gini_focus))/2) + rep(c(0,1.5),length(Gini_focus$Gini_focus)/2), 
     cex = 0.75,
     labels = annot10k$gene[annot10k$ID %in% rownames(Gini_focus)], srt = 75, col = colorlist)
dev.off()


pdf(file.path(outdir,"GiniScores_H3K4me3_5FU6_vs_DMSO_dotplot.pdf"),width = 12)
plot(Gini_focus_5FU6$Gini_focus~Gini_focus_DMSOi$Gini_focus,col=colorlist,pch=19,ylim=c(0.3,1),xlim=c(0.3,1))
abline(a=0,b=1)
text(Gini_focus_5FU6$Gini_focus~Gini_focus_DMSOi$Gini_focus,
     
     cex = 0.75,
     labels = annot10k$gene[annot10k$ID %in% rownames(Gini_focus)], srt = 75, col = colorlist)
dev.off()

#Heatmap on top cells of DMSOi and 5FU6
tmat <- cbind(mat_bivalent_5FU6,mat_bivalent_DMSOi)
distHC <- c("distPearson","distCosine","euclidean","maximum","manhattan","canberra","binary","minkowski")[3]
methHC <- c("ward","ward.D","ward.D2","single","complete","average")[3]
if(distHC=="distPearson")	d <- distPearson(tmat)
if(distHC=="distCosine")	d <- distCosine(tmat)
if(distHC %in% c("euclidean","maximum","manhattan","canberra","binary","minkowski"))	
  d <- dist(tmat,distHC)
hc <- hclust(d,method=methHC)
chRangeHM = TRUE
mat.so <- as.matrix(cor_genes[hc_genes$order,hc_genes$order])
if(chRangeHM){
  for(i in 1:nrow(mat.so)){
    mat.so[i,] <- geco.changeRange(mat.so[i,],newmin=0,newmax=1)
  }
}
png(file.path(outdir,paste0("Cooccurrence_pers_genes_heatmap_",name,".png")), height=8000,width=8000,res=600)
geco.hclustAnnotHeatmapPlot(x=(mat.so),
                            hc=hc,
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


###############################################################
#K27me3 study
###############################################################

#Use peak annot to select top cells and then try to switch to 10k
annotZerone
M5FU6_day33 <- read.table("MM468_K27_peaks_1000/MM468_5FU6_day33.tsv")
DMSO6_day60 <- read.table("MM468_K27_peaks_1000/MM468_DMSO6_day60.tsv")


loci_control <- c("chr7_127610001_127640000","chr7_127610001_127640000","chr19_34686001_34730000","chr19_13276001_13725000","chr19_15609001_16005000","chr19_50514001_50788000","chr20_16764001_17491000")
gene_control <- c("CACNA1A","PAX4","ZNF302","CYP4F2","SHANK1","PCSK2")

mat <- DMSO6_day60
mat_loci_control <- mat[loci_control,]
bin_mat = mat_loci_control
bin_mat[(bin_mat<2)]=0
bin_mat[(bin_mat>1)]=1

coocurrence_K4me3_score_housekeeping = as.data.frame(colSums(bin_mat) / nrow(bin_mat))
top_cells <- rownames(coocurrence_K4me3_score_housekeeping)[which(coocurrence_K4me3_score_housekeeping$`colSums(bin_mat)/nrow(bin_mat)`>0.4)]


annotZerone$Gini_K4_DMSO6_day60 <- Gin_unt
Gini_focus = Gin_unt[which(names(Gin_unt) %in% c(house_keeping_genes,overexpressed_pers_genes))]
Gini_focus = data.frame(Gini_focus,0)

colorlist = ifelse(rownames(Gini_focus) %in% house_keeping_genes, "grey", "red")
pdf(file.path(outdir,"GiniScores_common_pers_genes_hk_inInitial.pdf"),width = 12)
plot(0,0,type="n",xlim=c(-0.06,1.06), ylim=c(0,2), yaxt = 'n', main ="Gini Scores")
for (i in 1:length(Gini_focus$Gini_focus)) { stripchart(Gini_focus$Gini_focus[i],
                                                        add = T, bg = colorlist[i],
                                                        vertical = F, pch =21, method = "jitter") }
text(Gini_focus$Gini_focus,
     (runif(length(Gini_focus$Gini_focus))/2) + rep(c(0,1.5),length(Gini_focus$Gini_focus)/2), 
     cex = 0.75,
     labels = rownames(Gini_focus), srt = 75, col = colorlist)
dev.off()
