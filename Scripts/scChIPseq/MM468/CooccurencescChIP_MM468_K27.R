######## GINI coefficients MM468 - scChIP ############

library(here)
source(file.path(here(),"Scripts","functions.R"))
source(file.path(here(),"Scripts","global_var.R"))

# Parameters
percent_coocurrence_to_define_top_cells = 0.3

maindir = here()
outdir = file.path(maindir,"output","scChIPseq","MM468","Coocurrence_K27")
if(!dir.exists(outdir)) dir.create(outdir)

# Find association Gene to Peak by using MM468 K4 peaks "consensus annotation" obtained from bulk ChIPseq 
annot_MM468_K27 <- rtracklayer::import(
  file.path(maindir,"annotation","MM468_peaks_K27.bed.gz"))
start(annot_MM468_K27) <- start(annot_MM468_K27) -1
annot_MM468_K27$ID = paste0(seqnames(annot_MM468_K27),":",
                            start(annot_MM468_K27),"-",
                            end(annot_MM468_K27))

load(file.path(maindir,"annotation","hg38.GeneTSS.rda"))

hg38.GeneTSS = as(hg38.GeneTSS,"GRanges")
annot10k <- read.table(
  unz(file.path(
    maindir,"annotation","gencode.v34.annotation.transcriptTSS_10k.bed.zip"),
    "gencode.v34.annotation.transcriptTSS_10k.bed"))
colnames(annot10k) <- c("chr","start","end","transcripts","gene","sense")
annot10k = as(annot10k[,c(1:3,5)],"GRanges")
annot10k$ID = paste0(seqnames(annot10k),":",start(annot10k),"-",end(annot10k))

differential_results_scRNA <- readxl::read_xlsx(file.path(
  maindir,"output","scRNAseq","MM468","Persister","Supervised","Tables",
  "Differential_analysis_Limma_logFC_1.58.xlsx"),
  sheet=1)

overexpressed_MM468_persister_genes = differential_results_scRNA$Symbol
transcripts = annot10k[which(annot10k$gene %in%
                               overexpressed_MM468_persister_genes)]
regions = hg38.GeneTSS[which(hg38.GeneTSS$gene %in%
                               overexpressed_MM468_persister_genes)]
regions = transcripts[subjectHits(findOverlaps(regions,transcripts)),]
overexpressed_MM468_persister_regions = regions$ID

annot_MM468_K27_overexpressed = annot_MM468_K27[queryHits(findOverlaps(annot_MM468_K27,transcripts)),]
annot_MM468_K27_overexpressed$gene = transcripts$gene[subjectHits(findOverlaps(annot_MM468_K27,transcripts))]
annot_MM468_K27_overexpressed = annot_MM468_K27_overexpressed[-which(duplicated(annot_MM468_K27_overexpressed)),]

house_keeping_genes = c("ZNF574","ZNF181","ZNF423","OR4F15","NLRP4",
                         "reg1","reg2","reg3","reg4","reg5")

house_keeping_genes_region = list(
  c("chr19",42064001,42071000),
  c("chr19",34686001,34730000),
  c("chr16",48633001,50023000),
  c("chr15",101765001,101888000),
  c("chr19",55707001,56072000),
  c("chr17",5497001,6433000),
  c("chr17",10247001,10659000),
  c("chr17",32885001,34230000),
  c("chr15",98086001,98438000),
  c("chr8",46102001,46284000)
)
house_keeping_genes_region <- as(
  setNames(as.data.frame(t(as.data.frame(house_keeping_genes_region))),
           c("chr","start","end")),"GRanges")
house_keeping_genes_region$gene = house_keeping_genes
annot_MM468_K27$housekeeping = ""
annot_MM468_K27$housekeeping[subjectHits(
  findOverlaps(house_keeping_genes_region, annot_MM468_K27))] = "housekeeping"

house_keeping_genes_region$ID=paste0(seqnames(house_keeping_genes_region),":",
                                       start(house_keeping_genes_region),"-",
                                       end(house_keeping_genes_region))

selected_regions = c(annot_MM468_K27_overexpressed$ID,
                     house_keeping_genes_region$ID)
combined = c(annot_MM468_K27_overexpressed, house_keeping_genes_region)
selected_genes = combined$gene[match(selected_regions,combined$ID)]

# Load bulk ChIP K27
load(file.path(maindir,"output","bulk_ChIPseq","MM468","K27_peaks_K27",
               "Supervised","RData","Supervised_res_Limma_object.Rdata"))
dim(res)
res_ids = paste0(res$chr,":",res$start,"-",res$end)

mark = list("K27" = c())
i=0
for(archive in c("MM468_K27_peaks_1000.zip")){
  i = i +1
  file_names = basename(unzip(file.path(
    maindir,"input","scChIPseq","MM468","Count_Matrices",archive),
    list = T, exdir = tempdir())$Name)
  list_files = vector("list", length(file_names[grep(".tsv",file_names)]))
  names(list_files) = gsub(".tsv","",file_names[grep(".tsv",file_names)])
  path_matrices = file.path(  maindir,"input","scChIPseq",
                              "MM468","Count_Matrices",
                              archive)
  
  #10k counts
  for(file in names(list_files)) {
    list_files[[file]] = scater::readSparseCounts(
      unzip(zipfile = path_matrices,files = paste0(file,".tsv"), 
            exdir = tempdir()))
    colnames(list_files[[file]]) = paste0(file,"_",
                                          colnames(list_files[[file]]))
  }
  mark[[i]] = list_files
}

## Gini

# Determine the cells with the top cells with the most efficient IP 
# (the most presence of H3K4me3 fragments for housekeeping genes)
mat <- mark$K27$MM468_DMSO6_day60
mat_housekeeping <- mat[which(rownames(mat) %in% house_keeping_genes_region$ID),]
bin_mat = mat_housekeeping
bin_mat[(bin_mat>0)]=1
coocurrence_K27me3_score_housekeeping = as.data.frame(
  Matrix::colSums(bin_mat) / nrow(bin_mat))
colnames(coocurrence_K27me3_score_housekeeping) = "Coocurrence_score"
top_cells_init <- rownames(coocurrence_K27me3_score_housekeeping)[which(
  coocurrence_K27me3_score_housekeeping$Coocurrence_score > 
    percent_coocurrence_to_define_top_cells)]

mat_top  <- mat[,top_cells_init]
Gin_unt = edgeR::gini(t(mat_top))
Gini_focus = Gin_unt[which(names(Gin_unt) %in% selected_regions)]
Gini_focus = data.frame(Gini_focus,0)
colorlist = ifelse(rownames(Gini_focus) %in%
                     house_keeping_genes_region$ID, "grey", "red")

mat_overexpressed <- mat_top[annot_MM468_K27_overexpressed$ID,]

bin_mat = mat_overexpressed
bin_mat[(bin_mat>0)]=1
coocurrence_K27me3_score_bivalent_DMSOi = as.data.frame(
  Matrix::colSums(bin_mat) / nrow(bin_mat))
colnames(coocurrence_K27me3_score_bivalent_DMSOi) = "Coocurrence_score"

pdf(file.path(outdir,"GiniScores_H3K27me3_bivalent_genes_inK27_DMSOi.pdf"),
    width = 12)
plot(0,0,type="n",xlim=c(-0.06,1.06), ylim=c(0,2),
     yaxt = 'n', main ="Gini Scores")
for (i in 1:length(Gini_focus$Gini_focus)) { 
  stripchart(Gini_focus$Gini_focus[i],
             add = T, bg = colorlist[i],
             vertical = F, pch =21, method = "jitter") }
text(Gini_focus$Gini_focus,
     (runif(length(Gini_focus$Gini_focus))/2) +
       rep(c(0,1.5),length(Gini_focus$Gini_focus)/2), 
     cex = 0.75,
     labels = combined$gene[match(rownames(Gini_focus),combined$ID)],
     srt = 75, col = colorlist)
dev.off()
Gini_focus_unt <- Gini_focus

#Determine the cells with the top cells with the most efficient IP (the most presence of H3K4me3 fragments for housekeeping genes)
mat <- mark$K27$MM468_5FU6_day33
mat_housekeeping <- mat[which(rownames(mat) %in% house_keeping_genes_region$ID),]
bin_mat = mat_housekeeping
bin_mat[(bin_mat>0)]=1
coocurrence_K27me3_score_housekeeping = as.data.frame(
  Matrix::colSums(bin_mat) / nrow(bin_mat))
colnames(coocurrence_K27me3_score_housekeeping) = "Coocurrence_score"
top_cells_pers <- rownames(coocurrence_K27me3_score_housekeeping)[which(
  coocurrence_K27me3_score_housekeeping$Coocurrence_score > 
    percent_coocurrence_to_define_top_cells)]

mat_top  <- mat[,top_cells_pers]
Gin_unt = edgeR::gini(t(mat_top))
Gini_focus = Gin_unt[which(names(Gin_unt) %in% selected_regions)]
Gini_focus = data.frame(Gini_focus,0)
colorlist = ifelse(rownames(Gini_focus) %in%
                     house_keeping_genes_region$ID, "grey", "red")
Gini_focus_DMSOi <- Gini_focus

mat_overexpressed <- mat_top[annot_MM468_K27_overexpressed$ID,]

bin_mat = mat_overexpressed
bin_mat[(bin_mat>0)]=1
coocurrence_K27me3_score_bivalent_pers = as.data.frame(
  Matrix::colSums(bin_mat) / nrow(bin_mat))
colnames(coocurrence_K27me3_score_bivalent_pers) = "Coocurrence_score"

pdf(file.path(outdir,"GiniScores_H3K27me3_bivalent_genes_inK27_persister.pdf"),
    width = 12)
plot(0,0,type="n",xlim=c(-0.06,1.06), ylim=c(0,2),
     yaxt = 'n', main ="Gini Scores")
for (i in 1:length(Gini_focus$Gini_focus)) { 
  stripchart(Gini_focus$Gini_focus[i],
             add = T, bg = colorlist[i],
             vertical = F, pch =21, method = "jitter") }
text(Gini_focus$Gini_focus,
     (runif(length(Gini_focus$Gini_focus))/2) +
       rep(c(0,1.5),length(Gini_focus$Gini_focus)/2), 
     cex = 0.75,
     labels = combined$gene[match(rownames(Gini_focus),combined$ID)],
     srt = 75, col = colorlist)
dev.off()
Gini_focus_pers <- Gini_focus


pdf(file.path(outdir,"GiniScores_H3K4me3_5FU6_vs_DMSO_dotplot.pdf"),
    width = 12)
plot(Gini_focus_pers$Gini_focus~Gini_focus_unt$Gini_focus,
     col=colorlist,pch=19,ylim=c(0.3,1),xlim=c(0.3,1))
abline(a=0,b=1)
text(Gini_focus_pers$Gini_focus~Gini_focus_unt$Gini_focus,cex = 0.75,
     labels = combined$gene[combined$ID %in% rownames(Gini_focus)],
     srt = 75, col = colorlist)
dev.off()

#Heatmap on top cells of DMSOi and 5FU6
# selected_regions = c(overexpressed_MM468_persister_regions,house_keeping_genes_region$ID)
selected_cells = c(top_cells_pers,top_cells_init)

mat_pers_K27 = mark$K27$MM468_5FU6_day33[match(selected_regions,rownames(mark$K27$MM468_5FU6_day33)),]
mat_init_K27 = mark$K27$MM468_DMSO6_day60[match(selected_regions,rownames(mark$K27$MM468_DMSO6_day60)),]

#remove 0 cells
# mat_init_K27 = mat_init_K27[,-which(colSums(mat_init_K27) == 0)]
# mat_pers_K27 = mat_pers_K27[,-which(colSums(mat_pers_K27) == 0)]

mat_K27 = cbind(mat_init_K27,mat_pers_K27)

tmat = mat_K27
tmat = tmat[,match(selected_cells,colnames(tmat))]
bin_mat = tmat
bin_mat[bin_mat>0] = 1

methHC <- c("ward","ward.D","ward.D2","single","complete","average")[3]
hc <- hclust(as.dist(1-cor(as.matrix(bin_mat))),method=methHC)
mat.so <- as.matrix(bin_mat)

#Create cell annotation data.frame
annot_cell = data.frame(cell_id = selected_cells,
                        sample_id = c(rep("Persister",length(top_cells_pers)),
                                      rep("Initital",length(top_cells_init))))
rownames(annot_cell) = annot_cell$cell_id
anocol = geco.annotToCol4(annot_cell)

hc_gene = hclust(dist(as.matrix(bin_mat)),method=methHC)
hc_gene$labels = rep("",length(hc_gene$labels))

#Create gene/region annotation data.frame
annot_gene = data.frame(
  region = selected_regions,
  genes =  selected_genes,
  type = c(rep("Persister",length(annot_MM468_K27_overexpressed)),
           rep("Housekeeping",length(house_keeping_genes_region))),
  total_cells = rowSums(bin_mat))

genes = annot_gene$genes
i=1
counter =1
while(i<length(genes)){
  if(annot_gene$genes[i] == annot_gene$genes[i+1]){
    counter = counter+1
    genes[i+1] = paste0(annot_gene$genes[i],".",counter)
  } else{
    counter=1
  }
  i=i+1
}

rownames(annot_gene) = genes
rownames(mat.so) <- genes
anocol_gene = geco.annotToCol4(annot_gene)

png(file.path(outdir,paste0("Cooccurrence_pers_hk_topcells_pers_init_K27_heatmap.png")),
    height=8000, width=8000, res=600)
geco.hclustAnnotHeatmapPlot.withColumn(
  x=mat.so[hc_gene$order,hc$order],
  hc=hc,
  hmColors=hmColors,
  anocol=anocol[hc$order,],
  hc_row = hc_gene,
  anorow = anocol_gene[hc_gene$order,-c(1,2)],
  dendro.cex=0.01,
  xlab.cex=1,
  hmRowNames=TRUE,
  hmRowNames.cex=0.25,
  hmCategNamesRows = TRUE,
  hmCategNamesRows.cex = 0.35)
dev.off()

#### K27 init only #####
differential_results_scRNA_all <- readxl::read_xlsx(file.path(
  maindir,"output","scRNAseq","MM468","Persister","Supervised","Tables",
  "Differential_analysis_Limma_logFC_1.58.xlsx"),
  sheet=3)

# Load scRNA counts
load(file.path( maindir,"output","scRNAseq","MM468","Persister","Unsupervised",
                "RData","persister_RawCounts.RData"))
load(file.path( maindir,"output","scRNAseq","MM468","Persister","Unsupervised",
                "RData","persister_LogCounts.RData"))
non_treated_cell <- colnames(persister_LogCounts)[grep("initial",colnames(persister_LogCounts))]
persister_LogCounts <- persister_LogCounts[match(selected_genes[1:60],rownames(persister_LogCounts)),
                                           non_treated_cell]
persister_LogCounts = rbind(persister_LogCounts,matrix(0,nrow = 5,ncol=ncol(persister_LogCounts)))
rownames(persister_LogCounts)[61:65] = paste0("reg",1:5)
persister_LogCounts[61:65,1] = 1 

RawCounts =  RawCounts[match(selected_genes[1:60],rownames(RawCounts)),
                                 non_treated_cell]
RawCounts = rbind(RawCounts,matrix(0,nrow = 5,ncol=ncol(RawCounts)))
rownames(RawCounts)[61:65] = paste0("reg",1:5)
RawCounts[61:65,1] = 1 
  
scRNA_log2FC = differential_results_scRNA_all[match(selected_genes,
                                                    differential_results_scRNA_all$Symbol),"log2FC.C2_pers"]
scRNA_log2FC$log2FC.C2_pers[is.na(scRNA_log2FC)] <- 0

binary_scRNA = RawCounts
binary_scRNA[binary_scRNA>0] = 1

meanLevel_nonZeroCells <- sapply(rownames(persister_LogCounts), function(x){
  vec <- persister_LogCounts[x,]
  sel <- which(binary_scRNA[x,]>0)
  mean(vec[sel])
})
meanLevel_nonZeroCells[is.nan(meanLevel_nonZeroCells)] = 0

bulk_log2FC = res[match(selected_regions,res_ids),"log2FC.X5FU2_3_5"]

for(samp in c("MM468_DMSO3_day77","MM468_DMSO6_day60")){
  mat <- mark$K27[[samp]]
  rownames(mat) = annot_MM468_K27$ID
  mat_housekeeping <- mat[which(rownames(mat) %in% house_keeping_genes_region$ID),]
  bin_mat = mat_housekeeping
  bin_mat[(bin_mat>0)]=1
  coocurrence_K27me3_score_housekeeping = as.data.frame(
    Matrix::colSums(bin_mat) / nrow(bin_mat))
  colnames(coocurrence_K27me3_score_housekeeping) = "Coocurrence_score"
  top_cells_init <- rownames(coocurrence_K27me3_score_housekeeping)[which(
    coocurrence_K27me3_score_housekeeping$Coocurrence_score > 
      percent_coocurrence_to_define_top_cells)]
  
  selected_cells = top_cells_init
  tmat = mat
  tmat = tmat[selected_regions,match(selected_cells,colnames(tmat))]
  bin_mat = tmat
  bin_mat[bin_mat>0] = 1
  
  methHC <- c("ward","ward.D","ward.D2","single","complete","average")[3]
  hc <- hclust(as.dist(1-cor(as.matrix(bin_mat))),method=methHC)
  mat.so <- as.matrix(bin_mat)
  
  #Create cell annotation data.frame
  annot_cell = data.frame(cell_id = selected_cells,
                          sample_id = c(rep("Initital",length(top_cells_init))))
  rownames(annot_cell) = annot_cell$cell_id
  anocol = geco.annotToCol4(annot_cell)
  
  hc_gene = hclust(dist(as.matrix(bin_mat)),method=methHC)
  hc_gene$labels = rep("",length(hc_gene$labels))
  
  #Create gene/region annotation data.frame
  annot_gene = data.frame(
    region = selected_regions,
    genes = selected_genes,
    type = c(rep("Persister",length(annot_MM468_K27_overexpressed)),
             rep("Housekeeping",length(house_keeping_genes_region))),
    total_cells = rowSums(bin_mat),
    bulk_log2FC = bulk_log2FC,
    scRNA_log2FC = scRNA_log2FC$log2FC.C2_pers,
    scRNA_logCounts = rowSums(persister_LogCounts),
    scRNA_percentCells = rowSums(binary_scRNA)/ncol(binary_scRNA),
    scRNA_meanLevel = meanLevel_nonZeroCells
    )
  
  genes = annot_gene$genes
  i=1
  counter =1
  while(i<length(genes)){
    if(annot_gene$genes[i] == annot_gene$genes[i+1]){
      counter = counter+1
      genes[i+1] = paste0(annot_gene$genes[i],".",counter)
    } else{
      counter=1
    }
    i=i+1
  }
  
  rownames(annot_gene) = genes
  rownames(mat.so) <- genes
  anocol_gene = geco.annotToCol4(annot_gene)
  
  order_log2FC_bulk_K27 = annot_gene$bulk_log2FC[1:55] 
  order_log2FC_bulk_K27 = order(order_log2FC_bulk_K27)
  order_log2FC_bulk_K27 = c(order_log2FC_bulk_K27,56:65)
  mat.so = mat.so[order_log2FC_bulk_K27,]
  annot_gene = annot_gene[order_log2FC_bulk_K27,]
  
  anocol_gene = geco.annotToCol4(annot_gene)
  # anocol_gene. = geco.annotToCol4(annot_gene,scale_q = "viridis")
  # anocol_gene[,"bulk_log2FC"] = anocol_gene.[,"bulk_log2FC"]
  # anocol_gene[,"scRNA_log2FC"] = anocol_gene.[,"scRNA_log2FC"]
  
  for(ord in c("force_order", "hierarchical")){
    if(ord != "force_order") order = hc_gene$order else order =1:nrow(annot_gene)
  png(file.path(outdir,paste0("Cooccurrence_pers_hk_topcells_",samp,"_",ord,"K27_heatmap.png")),
      height=8000, width=8000, res=600)
  geco.hclustAnnotHeatmapPlot.withColumn(
    x=mat.so[order,hc$order],
    hc=hc,
    hmColors=hmColors,
    anocol=anocol[hc$order,],
    hc_row = hc_gene,
    anorow = anocol_gene[order,-c(1,2)],
    dendro.cex=0.01,
    xlab.cex=1,
    hmRowNames=TRUE,
    hmRowNames.cex=0.25,
    hmCategNamesRows = TRUE,
    hmCategNamesRows.cex = 0.35)
  dev.off()
  }
}

