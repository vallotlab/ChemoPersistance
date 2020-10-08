######## GINI coefficients MM468 - scChIP ############

library(here)
source(file.path(here(),"Scripts","functions.R"))
source(file.path(here(),"Scripts","global_var.R"))

# Parameters
percent_coocurrence_to_define_top_cells = 0.3

maindir = here()
outdir = file.path(maindir,"output","scChIPseq","MM468","Coocurrence_K4")
if(!dir.exists(outdir)) dir.create(outdir)

# Find association Gene to Peak by using MM468 K4 peaks "consensus annotation" obtained from bulk ChIPseq 
annot_MM468_K4 <- rtracklayer::import(file.path(maindir,"annotation","MM468_peaks_K4_29714_peaks.bed.gz"))
load(file.path(maindir,"annotation","hg38.GeneTSS.rda"))
hg38.GeneTSS = as(hg38.GeneTSS,"GRanges")
annot10k <- read.table(unz(file.path(maindir,"annotation","gencode.v34.annotation.transcriptTSS_10k.bed.zip"),
                           "gencode.v34.annotation.transcriptTSS_10k.bed"))
colnames(annot10k) <- c("chr","start","end","transcripts","gene","sense")
annot10k = as(annot10k[,c(1:3,5)],"GRanges")
annot10k$ID = paste0(seqnames(annot10k),":",start(annot10k),"-",end(annot10k))
  
differential_results_scRNA <- readxl::read_xlsx(file.path(
  maindir,"output","scRNAseq","MM468","Persister","Supervised","Tables",
  "Differential_analysis_Limma.xlsx"),
  sheet=1)

overexpressed_MM468_persister_genes = differential_results_scRNA$Symbol
transcripts = annot10k[which(annot10k$gene %in% overexpressed_MM468_persister_genes)]
regions = hg38.GeneTSS[which(hg38.GeneTSS$gene %in% overexpressed_MM468_persister_genes)]
regions = transcripts[subjectHits(findOverlaps(regions,transcripts)),]
overexpressed_MM468_persister_regions = regions$ID

# house_keeping_genes = c("GAPDH","VGF","CCNA2","LMNA","GAPDH","UBB","UBC")
# From  https://www.genomics-online.com/resources/16/5049/housekeeping-genes/
# Eisenberg, Levanon: “Human housekeeping genes, revisited.

house_keeping_genes = c("TFRC", "SDHA", "TBP", "ACTB", "PPIA", "GUSB", "YWHAZ",
                        "HMBS", "GAPDH", "RPLP0", "B2M", "RPL13A", "PGK1","HPRT1")
house_keeping_genes_row <- c(11477, 13853, 19080, 19230, 19841, 20044, 23094,
                             30210, 30751, 32985, 36858, 46622, 51202, 51630) #manually annotated with the right TSS

annot10k_housekeeping <- annot10k[house_keeping_genes_row,]
house_keeping_genes_region = annot10k_housekeeping$ID

mark = list("K4" = c() , "K27" = c())
i=0
for(archive in c("MM468_K4_transcripts_10k_1000.zip","MM468_K27_transcripts_10k_1000.zip")){
  i = i +1
  file_names = basename(unzip(file.path(
    maindir,"input","scChIPseq","MM468","Count_Matrices",archive),
    list = T)$Name)
  list_files = vector("list", length(file_names[grep(".tsv",file_names)]))
  names(list_files) = gsub(".tsv","",file_names[grep(".tsv",file_names)])
  path_matrices = file.path(  maindir,"input","scChIPseq",
                              "MM468","Count_Matrices",
                              archive)
  
  #10k counts
  for(file in names(list_files)) {
    list_files[[file]] = scater::readSparseCounts(
      unzip(zipfile = path_matrices,files = paste0(file,".tsv")))
    colnames(list_files[[file]]) = paste0(file,"_",colnames(list_files[[file]]))
  }
  mark[[i]] = list_files
}

## Gini
bivalent_genes <- c("BMP6","LAMB1", "VIM", "FOSL1", "ABCC4",
                    "COL4A2", "KRT14", "TGFB1", "KLK10" )
bivalent_genes_row <- c(16677, 20736, 25943, 29195, 34181, 34369, 41923, 46164, 46727)
annot10k_bivalent <- annot10k[bivalent_genes_row,]

#Determine the cells with the top cells with the most efficient IP (the most presence of H3K4me3 fragments for housekeeping genes)
mat <- list_files$MM468_DMSO6_D0_K4
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

#Determine the cells with the top cells with the most efficient IP (the most presence of H3K4me3 fragments for housekeeping genes)
mat <- list_files$MM468_5FU6_D60_K4
mat_housekeeping <- mat[house_keeping_genes_row,]
bin_mat = mat_housekeeping
bin_mat[(bin_mat>0)]=1
coocurrence_K4me3_score_housekeeping = as.data.frame(colSums(bin_mat) / nrow(bin_mat))
colnames(coocurrence_K4me3_score_housekeeping) = "Coocurrence_score"

top_cells_pers <- rownames(coocurrence_K4me3_score_housekeeping)[which(
  coocurrence_K4me3_score_housekeeping$Coocurrence_score>percent_coocurrence_to_define_top_cells)]

mat_top  <- mat[,top_cells_pers]
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
text(Gini_focus_5FU6$Gini_focus~Gini_focus_DMSOi$Gini_focus,cex = 0.75,
     labels = annot10k$gene[annot10k$ID %in% rownames(Gini_focus)], srt = 75, col = colorlist)
dev.off()


# Violin plot co-occurence score
mat_pers = list_files$MM468_5FU6_D60_K4
mat_init = list_files$MM468_DMSO6_D0_K4

png(file.path(outdir,paste0("Coocurrence_score_overpers_K4_violin.png")),height=1000,width=1500,res=300)
violin_plot_cooccurence(mat_init,mat_pers,
                        sub_regions =  overexpressed_MM468_persister_regions,
                        sub_cells = list(top_cells_init,top_cells_pers))
dev.off()

png(file.path(outdir,paste0("Coocurrence_score_hk_K4_violin.png")),height=1000,width=1500,res=300)
violin_plot_cooccurence(mat_init,mat_pers,
                        sub_regions =  house_keeping_genes_region,
                        sub_cells = list(top_cells_init,top_cells_pers))
dev.off()

#Heatmap on top cells of DMSOi and 5FU6
mat_pers = list_files$MM468_5FU6_D60_K4
mat_init = list_files$MM468_DMSO6_D0_K4

mat_init_K27 = mark$K27$MM468_Ipop_DMSO6_K27me3
mat_pers_K27 = mark$K27$MM468_5FU6_D33_K27

tmat = Matrix::cBind(mat_pers,mat_init)
tmat = tmat[match(c(overexpressed_MM468_persister_regions,house_keeping_genes_region),rownames(tmat)),
            match(c(top_cells_pers,top_cells_init),colnames(tmat))]
bin_mat = tmat
bin_mat[bin_mat>0] = 1

methHC <- c("ward","ward.D","ward.D2","single","complete","average")[3]
hc <- hclust(as.dist(1-cor(as.matrix(bin_mat))),method=methHC)
mat.so <- as.matrix(bin_mat)

#Create cell annotation data.frame
annot_cell = data.frame(cell_id = c(top_cells_pers,top_cells_init),
                        sample_id = c(rep("Persister",length(top_cells_pers)),
                                      rep("Initital",length(top_cells_init))))
rownames(annot_cell) = annot_cell$cell_id
anocol = geco.annotToCol4(annot_cell)

hc_gene = hclust(dist(as.matrix(bin_mat)),method=methHC)
hc_gene$labels = rep("",length(hc_gene$labels))

#Create gene/region annotation data.frame
annot_gene = data.frame(region = c(overexpressed_MM468_persister_regions,house_keeping_genes_region),
                        genes = c(as.character(regions$gene),as.character(annot10k_housekeeping$gene)),
                        type = c(rep("Persister",length(overexpressed_MM468_persister_regions)),
                                 rep("Housekeeping",length(house_keeping_genes_region))),
                        total_cells = rowSums(bin_mat))
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


