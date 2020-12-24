######## GINI coefficients MM468 - scChIP ############

library(here)
source(file.path(here(),"Scripts","functions.R"))
source(file.path(here(),"Scripts","global_var.R"))

# Parameters
percent_coocurrence_to_define_top_cells = 0.5

maindir = here()
outdir = file.path(maindir,"output","scChIPseq","MM468","Coocurrence_K4")
if(!dir.exists(outdir)) dir.create(outdir)

metadata = data.frame(sample_id = c(
  "MM468_DMSO1_day0", "MM468_5FU1_day60"),
  sample_id_color = c("#afafafff", "#1AB8AD"))# "#009688ff"))

# Find association Gene to Peak by using MM468 K4 peaks "consensus annotation" obtained from bulk ChIPseq

# Load Transcripts of hg38, merged if closer to 10k:
annot10k <- read.table(
  unz(file.path(
    maindir,"annotation","gencode.v34.transcripts10k_K27.zip"),
                           "gencode.v34.transcripts10k_K27.bed"), sep="\t")[,-c(6,9)]
colnames(annot10k) <- c("chr","start","end","transcripts","gene","log2FC_bulk","qvalue_bulk","K27_status")
annot10k = as(annot10k,"GRanges")
annot10k$ID = paste0(seqnames(annot10k),":",start(annot10k),"-",end(annot10k))
annot10k$ID_bis = paste0(seqnames(annot10k),"_",start(annot10k),"_",end(annot10k))

# Load Genes overexpressed between Persister and Non treated  in the scRNAseq data
differential_results_K4 <- readxl::read_xlsx(file.path(
  maindir,"output","bulk_ChIPseq","MM468","K4_transcripts_10k","Supervised","Tables",
  "Differential_analysis_Limma.xlsx"))

differential_results_scRNA <- readxl::read_xlsx(file.path(
  maindir,"output","scRNAseq","MM468","Persister","Supervised","Tables",
  "Differential_analysis_Limma_logFC_1.58.xlsx"),sheet = 1)

persister_K27_regulated = readxl::read_xls(file.path(
  maindir,"output","scRNAseq","MM468","Persister","ComparisonChIPseq",
  "K27","persister_K27_status.xls"),sheet = 1)
persister_K27_regulated_genes = persister_K27_regulated$Symbol[which(persister_K27_regulated$K27_status %in% c("Depleted FC < -3","Depleted FC < -2"))]
persister_K27_regulated_regions = annot10k_byGene$ID[which(annot10k_byGene$gene %in% persister_K27_regulated_genes)]

# differential_results_K4$ID = sub("-",":",gsub("_","-",differential_results_K4$ID))
overexpressed_MM468_persister_genes = differential_results_scRNA$Symbol  

# Overexpressed & K27 regulated genes
annot_transcript_10k_over = annot10k[which(annot10k$gene %in% overexpressed_MM468_persister_genes),]
annot_transcript_10k_over = as.data.frame(annot_transcript_10k_over[which(annot_transcript_10k_over$K27_status!= "Not Overlapping"),])
annot10k_K27_byGene = annot_transcript_10k_over %>% group_by(gene) %>%  slice_max(order_by = abs(log2FC_bulk), n = 1) # Select only the transcript associated with top differential K27 log2FC 

table(annot10k_K27_byGene$K27_status)

overexpressed_MM468_persister_region_K27 = annot10k_K27_byGene$ID
overexpressed_MM468_persister_genes_K27 = annot10k_K27_byGene$gene

# Retrieve K27 genes that are deregulated but not overexpressed
load(file.path(maindir,"output","scRNAseq","MM468","Persister","Supervised","RData","Supervised_res_object_edgeR.Rdata"))
res.scRNA = my.res

annot10k_K27_ovlp = as.data.frame(annot10k)[which(annot10k$K27_status !="Not Overlapping"),]
annot10k_K27_byGene_all = annot10k_K27_ovlp %>% group_by(gene) %>%  slice_max(order_by = abs(log2FC_bulk), n = 1) # Select only the transcript associated with top differential K27 log2FC 

subset_int = annot10k_K27_byGene_all[annot10k_K27_byGene_all$log2FC_bulk < -log2(3) & annot10k_K27_byGene_all$qvalue_bulk < 0.1,]
subset_int = subset_int[which(subset_int$gene %in% my.res$Symbol),]
subset_int$scRNA_status = "Not Overexpressed"
subset_int$log2FC_scRNA = res.scRNA$log2FC.C2_pers[match(subset_int$gene, res.scRNA$Symbol)]
subset_int$qvalue_scRNA = res.scRNA$qval.C2_pers[match(subset_int$gene, res.scRNA$Symbol)]

subset_int$scRNA_status[which(subset_int$log2FC_scRNA > log2(1.5))] = "Overexpressed - FC > 1.5"
subset_int$scRNA_status[which(subset_int$log2FC_scRNA > log2(2))] = "Overexpressed - FC > 2"
subset_int$scRNA_status[which(subset_int$log2FC_scRNA > log2(3))] = "Overexpressed - FC > 3"

subset_int$scRNA_status[which(subset_int$log2FC_scRNA < -log2(2) )] = "Underexpressed - FC < -2"

table(subset_int$scRNA_status)
not_overexpressed_but_deregulated_K27_regions = subset_int$ID[subset_int$scRNA_status=="Not Overexpressed"]
not_overexpressed_but_deregulated_K27_genes = subset_int$gene[subset_int$scRNA_status=="Not Overexpressed"]

# Load house keeping genes
# From  https://www.genomics-online.com/resources/16/5049/housekeeping-genes/
# Eisenberg, Levanon: “Human housekeeping genes, revisited.

# house_keeping_genes = c("TFRC", "SDHA", "TBP", "ACTB", "PPIA", "GUSB", "YWHAZ",
#                         "HMBS", "GAPDH", "RPLP0", "B2M", "RPL13A", "PGK1",
#                         "HPRT1")
house_keeping_genes=c("SDHA","ACTB","PPIA","YWHAZ","RPL13A")
#manually annotated with the right TSS
house_keeping_genes_row <- c(13853, 19230, 19841,23094, 46622) 
  
annot10k_housekeeping <- annot10k[house_keeping_genes_row,]
house_keeping_genes_region = annot10k_housekeeping$ID

selected_regions <- c(house_keeping_genes_region,overexpressed_MM468_persister_region_K27,not_overexpressed_but_deregulated_K27_regions)
selected_genes <- c(house_keeping_genes,overexpressed_MM468_persister_genes_K27,not_overexpressed_but_deregulated_K27_genes)

# Find transcripts corresponding to overexpressed "persister genes"
annot10k_selected = annot10k[match(selected_regions, annot10k$ID),]

mark = list()
i=0
dir = file.path(maindir,"input","scChIPseq","MM468","Count_Matrices")
for(files in list.files(dir,pattern = ".*H3K4me3_TSS.tsv.gz")){
  
  mat = read.csv(gzfile(file.path(dir,files)),sep="\t")
  name = gsub("_H3K4me3_TSS.tsv.gz","",files)
  colnames(mat) = paste0(name,"_", colnames(mat))
  mat = mat[,-1]
  #10k counts
  mark[[name]] = mat
}

# Determine the cells with the top cells with the most efficient IP 
# (the most presence of H3K4me3 fragments for housekeeping genes)
runGini = TRUE
if(runGini){

  mat <- mark$MM468_DMSO1_day0
  mat_housekeeping <- mat[which(annot10k$ID %in% house_keeping_genes_region) ,]
  bin_mat = mat_housekeeping
  bin_mat[(bin_mat>0)]=1
  coocurrence_K4me3_score_housekeeping = as.data.frame(
    Matrix::colSums(bin_mat) / nrow(bin_mat))
  colnames(coocurrence_K4me3_score_housekeeping) = "Coocurrence_score"
  top_cells_init <- rownames(coocurrence_K4me3_score_housekeeping)[which(
    coocurrence_K4me3_score_housekeeping$Coocurrence_score > 
      percent_coocurrence_to_define_top_cells)]
  
  mat_top  <- mat[,top_cells_init]
  Gin_unt = edgeR::gini(t(mat_top))
  Gini_focus = Gin_unt[match(selected_regions,annot10k$ID)]
  Gini_focus = data.frame(Gini_focus,0)
  colorlist = ifelse(selected_regions %in%
                       annot10k_housekeeping$ID, "grey", "red")
  Gini_focus_DMSOi <- Gini_focus
  
  mat_over <- mat_top[which((annot10k$ID %in% overexpressed_MM468_persister_region_K27)),]
  mat_over_DMSOi <- mat_over
  
  bin_mat = mat_over_DMSOi
  bin_mat[(bin_mat>0)]=1
  coocurrence_K4me3_score_bivalent_DMSOi = as.data.frame(
    Matrix::colSums(bin_mat) / nrow(bin_mat))
  colnames(coocurrence_K4me3_score_bivalent_DMSOi) = "Coocurrence_score"
  
  # pdf(file.path(outdir,"GiniScores_K4_DMSOi_ChIP_based.pdf"),
  #     width = 12)
  # plot(0,0,type="n",xlim=c(-0.06,1.06), ylim=c(0,2),
  #      yaxt = 'n', main ="Gini Scores")
  # for (i in 1:length(Gini_focus$Gini_focus)) { 
  #   stripchart(Gini_focus$Gini_focus[i],
  #              add = T, bg = colorlist[i],
  #              vertical = F, pch =21, method = "jitter") }
  # text(Gini_focus$Gini_focus,
  #      (runif(length(Gini_focus$Gini_focus))/2) +
  #        rep(c(0,1.5),length(Gini_focus$Gini_focus)/2), 
  #      cex = 0.75,
  #      labels = annot10k$gene[annot10k$ID %in% selected_regions],
  #      srt = 75, col = colorlist)
  # dev.off()
  
  #Determine the cells with the top cells with the most efficient IP (the most presence of H3K4me3 fragments for housekeeping genes)
  mat <- mark$MM468_5FU1_day60
  mat_housekeeping <- mat[which(annot10k$ID %in% house_keeping_genes_region),]
  bin_mat = mat_housekeeping
  bin_mat[(bin_mat>0)]=1
  coocurrence_K4me3_score_housekeeping = as.data.frame(
    Matrix::colSums(bin_mat) / nrow(bin_mat))
  colnames(coocurrence_K4me3_score_housekeeping) = "Coocurrence_score"
  top_cells_pers <- rownames(coocurrence_K4me3_score_housekeeping)[which(
    coocurrence_K4me3_score_housekeeping$Coocurrence_score > 
      percent_coocurrence_to_define_top_cells)]
  
  mat_top  <- mat[,top_cells_pers]
  Gin_unt = edgeR::gini(t(mat_top))
  Gini_focus = Gin_unt[match(selected_regions,annot10k$ID)]
  Gini_focus = data.frame(Gini_focus,0)
  colorlist = ifelse(selected_regions %in%
                       annot10k_housekeeping$ID, "grey", "red")
  Gini_focus_pers <- Gini_focus
  
  mat_over <- mat_top[match(overexpressed_MM468_persister_region_K27,annot10k$ID),]
  mat_over_pers <- mat_over
  
  bin_mat = mat_over_pers
  bin_mat[(bin_mat>0)]=1
  coocurrence_K4me3_score_pers = as.data.frame(
    Matrix::colSums(bin_mat) / nrow(bin_mat))
  colnames(coocurrence_K4me3_score_pers) = "Coocurrence_score"
  
  # pdf(file.path(outdir,"GiniScores_K4_DMSOi_ChIP_based.pdf"),
  #     width = 12)
  # plot(0,0,type="n",xlim=c(-0.06,1.06), ylim=c(0,2),
  #      yaxt = 'n', main ="Gini Scores")
  # for (i in 1:length(Gini_focus$Gini_focus)) { 
  #   stripchart(Gini_focus$Gini_focus[i],
  #              add = T, bg = colorlist[i],
  #              vertical = F, pch =21, method = "jitter") }
  # text(Gini_focus$Gini_focus,
  #      (runif(length(Gini_focus$Gini_focus))/2) +
  #        rep(c(0,1.5),length(Gini_focus$Gini_focus)/2), 
  #      cex = 0.75,
  #      labels = annot10k$gene[which(annot10k$ID %in% selected_regions)],
  #      srt = 75, col = colorlist)
  # dev.off()
  # 
  # 
  # pdf(file.path(outdir,"GiniScores_H3K4me3_5FU6_vs_DMSO_dotplot_ChIPbased.pdf"),
  #     width = 12)
  # plot(Gini_focus_pers$Gini_focus~Gini_focus_DMSOi$Gini_focus,
  #      col=colorlist,pch=19,ylim=c(0.3,1),xlim=c(0.3,1))
  # abline(a=0,b=1)
  # text(Gini_focus_pers$Gini_focus~Gini_focus_DMSOi$Gini_focus,cex = 0.75,
  #      labels = annot10k$gene[annot10k$ID %in% rownames(Gini_focus)],
  #      srt = 75, col = colorlist)
  # dev.off()
   
  
  # Violin plot co-occurence score
  mat_pers = mark$MM468_5FU1_day60[which(annot10k$ID %in% overexpressed_MM468_persister_region_K27),]
  mat_init = mark$MM468_DMSO1_day0[which(annot10k$ID %in% overexpressed_MM468_persister_region_K27),]
  
  png(file.path(outdir,paste0("Coocurrence_score_overpers_K4_violin.png")),
      height=1000,width=1500,res=300)
  print(violin_plot_cooccurence(mat_init,mat_pers,
                          sub_cells = list(top_cells_init,top_cells_pers)))
  dev.off()
}

################################################################################
##  Heatmaps Non Treated :
# Prepare gene/region annotation :
differential_results_K4_all <- readxl::read_xlsx(file.path(
  maindir,"output","bulk_ChIPseq","MM468","K4_transcripts_10k","Supervised","Tables",
  "Differential_analysis_Limma.xlsx"), sheet = 3)

differential_results_K4_all$ID = sub("-",":",gsub("_","-",differential_results_K4_all$ID))
bulk_log2FC_K4 = differential_results_K4_all$log2FC.X5FU6[match(selected_regions,differential_results_K4_all$ID)]
  
for(samp in c("MM468_DMSO1_day0")){
  mat <- mark[[samp]]
  mat = cbind(mat, mark$MM468_5FU1_day60)
  
  mat_housekeeping <- mat[which(annot10k$ID %in% house_keeping_genes_region),]
  bin_mat = mat_housekeeping
  bin_mat[(bin_mat>0)]=1
  coocurrence_K4me3_score_housekeeping = as.data.frame(
    colSums(bin_mat) / nrow(bin_mat))
  colnames(coocurrence_K4me3_score_housekeeping) = "Coocurrence_score"
  
  top_cells <- rownames(coocurrence_K4me3_score_housekeeping)[which(
    coocurrence_K4me3_score_housekeeping$Coocurrence_score > 
      percent_coocurrence_to_define_top_cells)]
  
  tmat = mat[match(c(house_keeping_genes_region,persister_K27_regulated_regions), annot10k$ID),top_cells]
  bin_mat = tmat
  bin_mat[bin_mat>0] = 1
  methHC <- c("ward","ward.D","ward.D2","single","complete","average")[3]
  hc <- hclust(as.dist(1-cor(as.matrix(bin_mat))),method=methHC)
  mat.so <- as.matrix(bin_mat)
  
  #Create cell annotation data.frame
  annot_cell = data.frame(cell_id = top_cells,
                          sample_id = c(rep("Initial",length(grep("DMSO",top_cells))),
                          rep("Persister",length(grep("5FU",top_cells)))),
                          total_K4 = colSums(bin_mat))
  rownames(annot_cell) = annot_cell$cell_id
  anocol = geco.annotToCol4(annot_cell)
  
  hc_gene = hclust(dist(as.matrix(bin_mat)),method=methHC)
  hc_gene$labels = rep("",length(hc_gene$labels))
  
  #Create gene/region annotation data.frame
  annot_gene = data.frame(
    region = c(house_keeping_genes_region,persister_K27_regulated_regions),
    genes = c(house_keeping_genes,persister_K27_regulated_genes),
    HouseKeeping = c(rep("HouseKeeping",length(house_keeping_genes_region)),
                     rep("Persister_K27_regulated",length(persister_K27_regulated_genes)))
    )

  rownames(mat.so) = annot_gene$genes
  anocol_gene = geco.annotToCol4(annot_gene)
  hmColors <- colorRampPalette(c("white","royalblue"))(256)
  anocol[,"sample_id"] = c(rep("#afafafff",length(grep("DMSO",top_cells))),
                              rep("#1AB8AD",length(grep("5FU",top_cells))))

  pdf(file.path(outdir,paste0(samp,"_cooccurrence_pers_hk_topcells_K4_heatmap_K27regulated_pers_genes_with_nonover.pdf")),
      height=12, width=12)
  geco.hclustAnnotHeatmapPlot.withColumn(
    x=mat.so[,],
    hc=hc,
    hmColors=hmColors,
    anocol=anocol[,-c(1)],
    hc_row = hc_gene,
    anorow = anocol_gene[,-c(1)],
    dendro.cex=0.01,
    xlab.cex=1,
    hmRowNames=TRUE,
    hmRowNames.cex=0.25,
    hmCategNamesRows = TRUE,
    hmCategNamesRows.cex = 0.5)

  dev.off()
    
}

