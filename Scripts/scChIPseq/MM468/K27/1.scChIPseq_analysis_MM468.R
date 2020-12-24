library(here)
source(file.path(here(),"Scripts","functions.R"))
source(file.path(here(),"Scripts","global_var.R"))

maindir = here()

dataset_name = "MM468_H3K27me3_peaks"
ref_genome ="hg38"

metadata = data.frame(sample_id = c(
    "MM468_DMSO1_day60", "MM468_DMSO3_day77", "MM468_DMSO5_day131",
    "MM468_5FU1_day33", "MM468_5FU2_day67",
    "MM468_5FU6_day131", "MM468_5FU3_day147", "MM468_5FU2_day171"),
    sample_id_color = c("#dfdfdfff", "#999999ff","#363636",
                        "#118675ff", "#8cc453ff",
                        "#ff5722ff", "#feb40fff", "#fd8508ff"))

# 
# sample_id_color = c("#dfdfdfff", "#999999ff","#363636",
#                     "#118675ff", "#55D42F",
#                     "#ffea2eff", "#ff9800ff", "#ff5722ff"))

datadir = file.path(maindir, "input", "scChIPseq", "MM468")
outdir = file.path(maindir, "output", "scChIPseq", "MM468",dataset_name)
if(!dir.exists(outdir)) dir.create(outdir)
plotdir = file.path(outdir, "Plots"); if(!dir.exists(plotdir)) dir.create(plotdir)
plotdir_unsup = file.path(plotdir, "Unsupervised"); if(!dir.exists(plotdir_unsup)) dir.create(plotdir_unsup)
dir.create(file.path(plotdir_unsup,"Heatmaps"),showWarnings = F)
dir.create(file.path(plotdir_unsup,"UMAPS"),showWarnings = F)
plotdir_sup = file.path(plotdir, "Supervised"); if(!dir.exists(plotdir_sup)) dir.create(plotdir_sup)
dir.create(file.path(plotdir_sup,"Heatmaps"),showWarnings = F)
dir.create(file.path(plotdir_sup,"UMAPS"),showWarnings = F)

# Creating a ChromSCape directory in order for interactive vizualisation with
# ChromSCape Shiny Application (see ChromSCape::launchApp())
ChromSCape_analyses = file.path(maindir, "output", "scChIPseq","ChromSCape_analyses")
if(!dir.exists(ChromSCape_analyses)) dir.create(ChromSCape_analyses)
ChromSCape_directory = file.path(ChromSCape_analyses,dataset_name)
if(!dir.exists(ChromSCape_directory)) dir.create(ChromSCape_directory)
filt_dir <- file.path(ChromSCape_directory,"Filtering_Normalize_Reduce"); dir.create(filt_dir, showWarnings = FALSE)
cor_dir <- file.path(ChromSCape_directory, "correlation_clustering"); dir.create(cor_dir, showWarnings = FALSE)
cor_plot_dir <- file.path(cor_dir,"Plots"); dir.create(cor_plot_dir, showWarnings = FALSE)
da_gsa_dir <- file.path(ChromSCape_directory, "Diff_Analysis_Gene_Sets"); dir.create(da_gsa_dir, showWarnings = FALSE)


# Run ChromSCape :
set.seed(47)

################
## Input Data ##
################

# Reading in matrices
out <- import_scExp_gz(file.path(datadir,"Count_Matrices"),
                       pattern = "*_H3K27me3_peaks.tsv.gz")

# Save raw
datamatrix = out$datamatrix
annot_raw = out$annot_raw
save(datamatrix,annot_raw, file = file.path(ChromSCape_directory,"scChIP_raw.RData"))

scExp <- create_scExp(datamatrix, annot_raw)

#################
## Filter Data ##
#################

# Quality filtering
min_reads_per_cell = 3000
max_quantile_read_per_cell = 95
min_percent_to_keep_feature = 0 # Keep all peaks as peaks were pre-defined (not genomic bins!)

upper_threshold = quantile(Matrix::colSums(counts(scExp)),probs = max_quantile_read_per_cell/100 )

# Plot distribution before filtering
plot_distribution_scExp(scExp, raw = TRUE, log10 = TRUE, pseudo_counts = 1,bins = 150) +
    geom_vline(xintercept = log10(min_reads_per_cell),color = "green") +
    geom_vline(xintercept = log10(upper_threshold),color = "red")

# Filter
scExp = filter_scExp(scExp, min_cov_cell = min_reads_per_cell,
                     quant_removal = max_quantile_read_per_cell, 
                     percentMin = min_percent_to_keep_feature) 

# Plot distribution after filtering
plot_distribution_scExp(scExp, raw = T, log10 = T, pseudo_counts = 1,bins = 150)  +
    geom_vline(xintercept = log10(min_reads_per_cell),color = "green") +
    geom_vline(xintercept = log10(upper_threshold),color = "red")

# Number of cells passing filtering
num_cell_after_QC_filt_scExp(scExp, annot_raw)

# Excluding known CNA loci
exclude_regions = rtracklayer::import(file.path(maindir,"annotation","MM468_identified_CNA.bed.gz"))
exclude_regions = as.data.frame(exclude_regions)
scExp = exclude_features_scExp(scExp, exclude_regions, by = "region")

subsample_n = 1000
if(length(subsample_n)>0){
    print("Doing subsampling")
    set.seed(47)
    scExp = subsample_scExp(scExp, n_cells = subsample_n)
    gc()
    
}

# Normalize 
scExp = normalize_scExp(scExp, type = "CPM")

# Annotate with genes
scExp = feature_annotation_scExp(scExp, ref = ref_genome)

# Reduce dimensions
set.seed(47)
scExp = reduce_dims_scExp(scExp, n = 50, dimension_reductions = c("PCA","UMAP"),
                          batch_correction = F, verbose = F)

# Color samples
color_df = metadata
scExp = colors_scExp(scExp, c("sample_id","total_counts","batch_id"),color_df = color_df) # add colouring

# Deselect GSK & GSK 5FU
prefix <- paste0(dataset_name,"_",min_reads_per_cell,"_",
                 min_percent_to_keep_feature,"_",
                 max_quantile_read_per_cell,"_uncorrected")

save(scExp, file = file.path(filt_dir, paste0(prefix,".RData")))

# Development version of plot_reduced_dim_scExp function
png(file.path(plotdir_unsup,"UMAPS","PCA_sample_HQ.png"), res=300,height=1200,width=2000)
print(plot_reduced_dim_scExp_devel(scExp,"sample_id","PCA",downsample = 5000, transparency = 0.35))
dev.off()

png(file.path(plotdir_unsup,"UMAPS","UMAP_sample_HQ.png"), res=300,height=1200,width=2000)
print(plot_reduced_dim_scExp_devel(scExp,"sample_id","UMAP",downsample = 5000, transparency = 0.35))
dev.off()

# Clear shots
png(file.path(plotdir_unsup,"UMAPS","PCA_sample_HQ_clear.png"), res=300,height=1200,width=1200)
print(plot_reduced_dim_scExp_devel(scExp,"sample_id","PCA",downsample = 5000, transparency = 0.35) +
          theme(legend.position = "None", text = element_blank()))
dev.off()

png(file.path(plotdir_unsup,"UMAPS","UMAP_sample_HQ_clear.png"), res=300,height=1200,width=1200)
print(plot_reduced_dim_scExp_devel(scExp,"sample_id","UMAP",downsample = 5000, transparency = 0.35) +
          theme(legend.position = "None", text = element_blank()))
dev.off()

png(file.path(plotdir_unsup,"UMAPS","PCA_counts_HQ.png"), res=300,height=1200,width=2000)
print(plot_reduced_dim_scExp_devel(scExp,"total_counts","PCA",downsample = 5000, transparency = 0.35))
dev.off()

png(file.path(plotdir_unsup,"UMAPS","UMAP_counts_HQ.png"), res=300,height=1200,width=2000)
print(plot_reduced_dim_scExp_devel(scExp,"total_counts","UMAP",downsample = 5000, transparency = 0.35))
dev.off()

# Calculate correlation matrix and hierarchical clustering
scExp = correlation_and_hierarchical_clust_scExp(scExp)

# Saving filtered data with correlation matrix
prefix <- paste0(dataset_name,"_",min_reads_per_cell,"_",
                 min_percent_to_keep_feature,"_",
                 max_quantile_read_per_cell,"_uncorrected")

save(scExp, file = file.path(filt_dir, paste0(prefix,".RData")))
gc()

############################
## Correlation clustering ##
############################
set.seed(47)
# Filter lowly correlated cells
scExp_cf = filter_correlated_cell_scExp(scExp,
                                        random_iter = 5,
                                        corr_threshold = 99,
                                        percent_correlation = 1)
gc()

# Run consensus clustering (RAM heavy) - optional
runConsensus <- TRUE
if(runConsensus){
    scExp_cf = consensus_clustering_scExp(scExp_cf, reps = 50,
                                          maxK = 6,
                                          clusterAlg = "hc",
                                          prefix = file.path(cor_plot_dir, prefix))
    
    # Plot item consensus score
    plot_cluster_consensus_scExp(scExp_cf)
    
}

# Choosing 4 clusters
nclust = 4
scExp_cf = choose_cluster_scExp(scExp_cf, nclust = nclust, consensus = runConsensus)

# Coloring clusters
color_df = data.frame(cell_cluster = unique(scExp_cf$cell_cluster))
color_df$cell_cluster_color = c("#ffd700ff", "#ffb5c5ff", "#3cb371ff" ,"#668b8bff") # , "#473c8bff"
# "#3cb371ff", "#94524aff", "#ba55d3ff")
scExp_cf = colors_scExp(scExp_cf, annotCol = "cell_cluster",
                        color_by = c("cell_cluster"),color_df = color_df) # add colouring


save(scExp_cf, file = file.path(ChromSCape_directory, "correlation_clustering",
                                paste0(prefix,".RData"))) #save the data 

# Plotting
png(file.path(plotdir_unsup,"Heatmaps","ClusteringHeatMap.png"), res=300,height=4000,width=4000)
print(plot_heatmap_scExp(scExp_cf))
dev.off()

# Plotting
png(file.path(plotdir_unsup,"Heatmaps","Consensus_score.png"), res=300,height=4000,width=4000)
print( plot_cluster_consensus_scExp(scExp_cf))
dev.off()

#plot UMAP cluster
png(file.path(plotdir_unsup,"UMAPS","UMAP_cluster_HQ.png"), res=300,height=1200,width=2000)
print(plot_reduced_dim_scExp_devel(scExp_cf,"cell_cluster","UMAP",downsample = 5000, transparency = 0.35, size = 2))
dev.off()

png(file.path(plotdir_unsup,"UMAPS","UMAP_cluster_HQ_clear.png"), res=300,height=1200,width=1200)
print(plot_reduced_dim_scExp_devel(scExp_cf,"cell_cluster","UMAP",downsample = 5000, transparency = 0.35, size = 2) +
          theme(legend.position = "None", text = element_blank()))
dev.off()

#plot UMAP cluster
png(file.path(plotdir_unsup,"UMAPS","UMAP_sample_cf_HQ.png"), res=300,height=1200,width=2000)
print(plot_reduced_dim_scExp_devel(scExp_cf,"sample_id","UMAP",downsample = 5000, transparency = 0.35, size = 2))
dev.off()

png(file.path(plotdir_unsup,"UMAPS","UMAP_sample_cf_HQ_clear.png"), res=300,height=1200,width=1200)
print(plot_reduced_dim_scExp_devel(scExp_cf,"sample_id","UMAP",downsample = 5000, transparency = 0.35, size = 2) +
          theme(legend.position = "None", text = element_blank()))
dev.off()

############################
## Intra Correlation      ##
############################

load(file = file.path(ChromSCape_directory, "correlation_clustering",
                      paste0(prefix,".RData")))

#Table association
num_cell_in_cluster_scExp(scExp_cf) 
tab = table(scExp_cf$sample_id, scExp_cf$cell_cluster)
write.csv(tab, file.path(plotdir_unsup,"Heatmaps", "Contingency_sample_cluster_bargraph.csv"))

# Plot contingency bar graph - scChIP seq
annot = SingleCellExperiment::colData(scExp_cf)

annot = as_tibble(annot) %>% group_by(sample_id,sample_id_color,cell_cluster,cell_cluster_color) %>%
    summarise(count=n())

annot = left_join(annot,
                  annot %>% ungroup %>% group_by(sample_id) %>% summarise(sum = sum(count)),
                  by ="sample_id")
annot = annot %>% mutate(freq = count/sum)

annot = annot %>% mutate(sample_id=factor(sample_id,levels=metadata$sample_id)) 

png(file.path(plotdir_unsup,"Heatmaps", "Contingency_sample_cluster_bargraph.png"),width=800,height=1500,res=300)
print(ggplot(annot) + geom_bar(aes(y=freq,x=sample_id,fill=cell_cluster), stat="identity") +
          scale_fill_manual(values=color_df$cell_cluster_color) +
          theme_classic() + 
          theme(axis.text.x = element_text(angle = 90)) )
dev.off()

png(file.path(plotdir_unsup,"Heatmaps", "Contingency_sample_cluster_bargraph_wolabel.png"),width=1500,height=1500,res=300)
print(ggplot(annot) + geom_bar(aes(y=freq,x=sample_id,fill=cell_cluster), stat="identity") +
          scale_fill_manual(values=color_df$cell_cluster_color) +
          theme_classic() + 
          theme(text = element_blank()) )
dev.off()


## Intra-correlation score 
mati = t(reducedDim(scExp_cf,"PCA"))
cor_mat <- cor(mati)
hc_PCA <- hclust(as.dist(1 - cor_mat), method = "ward.D"); gc()
annot = scExp_cf@colData

intra_corr=data.frame()
for(i in unique(annot$cell_cluster)){
    cells = as.character(annot$cell_id[which(annot$cell_cluster==i)])
    tmp = cor_mat[cells,cells]
    tab = data.frame("cluster" = rep(i,ncol(tmp)), "intra_corr" = colMeans(tmp))
    intra_corr=rbind(intra_corr,tab)
}

mat.so.cor <- mati[,hc_PCA$order]
size_group <- as.numeric(table(annot$cell_cluster))
IntraCorr_test <- data.frame(expressionGroup=color_df$cell_cluster,pvalue=0)

samp = data.frame()
for(i in paste0("C",1:nclust)) {
    samp = rbind(samp,intra_corr[sample(which(intra_corr$cluster==i),500,replace = T),])
}

#ref sample for wilcox testing, C5 cluster
MAT <- samp[samp$cluster %in% c("C3","C4"),]
MAT_ref <- as.matrix(MAT[,-1])

for(i in 1:length(unique(annot$cell_cluster))){
    
    Group <- as.character(unique(annot$cell_cluster))[i]
    Col_clust <- as.character(color_df$cell_cluster_color[i])
    MAT <- samp[samp$cluster==Group,]
    test <- wilcox.test(as.matrix(MAT[,-1]),MAT_ref)
    IntraCorr_test$pvalue[i] <- test$p.value
    
    png(file.path(plotdir, paste0("Pearson_Correlation_Inferno_scores_",Group,".png")), height=1150,width=1000,res=300)
    MAT <- mat.so.cor[,colnames(mat.so.cor) %in% annot$cell_id[annot$cell_cluster==Group]]
    debut <- round(100-(1-min(cor(MAT)))/0.02,0)
    image(cor(MAT),col=(viridis::inferno(100)[debut:100]))
    dev.off()
}

gc()

WriteXLS(IntraCorr_test,file.path(plotdir_unsup,"Heatmaps","IntraCorr_test.xls"),row.names = T)


png(file.path(plotdir_unsup,"Heatmaps","ExpressionGroup_intracorrelation_violin.png"),height=1000,width=1000,res=300)
ggplot(samp,aes(x=cluster,y=intra_corr, fill=cluster)) + 
    geom_violin(alpha=0.8) + theme_classic() + scale_fill_manual(values=c("slateblue4","gold","paleturquoise4",
                                                                          "pink1","mediumseagreen","#94524A","mediumorchid")) +
    stat_summary(fun=median, geom="point", size=2, color="black")
dev.off()

# By sample
annot = SingleCellExperiment::colData(scExp_cf)
intra_corr=data.frame()
for(i in unique(annot$sample_id)){
    cells = as.character(annot$cell_id[which(annot$sample_id==i)])
    tmp = cor_mat[cells,cells]
    tab = data.frame("sample_id" = rep(i,ncol(tmp)), "intra_corr" = colMeans(tmp))
    intra_corr=rbind(intra_corr,tab)
}

mat.so.cor <- mati[,hc_PCA$order]

samp = data.frame()
for(i in unique(annot$sample_id)) {
    samp = rbind(samp,intra_corr[sample(which(intra_corr$sample_id==i),151,replace = T),])
}

png(file.path(plotdir_unsup,"Heatmaps", "Sample_id_intracorrelation_violin.png"),
    height=1000,width=1000,res=300)
ggplot(samp,aes(x=sample_id, y=intra_corr, fill=sample_id)) + 
    geom_violin(alpha=0.8) + theme_classic() +
    scale_fill_manual(values= unique(annot$sample_id_color)) +
    stat_summary(fun=median, geom="point", size=2, color="black") + theme(axis.text.x = element_text(angle=90))
dev.off()

## Inter-correlation score 
inter_corr = data.frame()
for(i in unique(annot$cell_cluster)){
    cells_i = as.character(annot$cell_id[which(annot$cell_cluster==i)])
    for(j in unique(annot$cell_cluster)){
            cells_j = as.character(annot$cell_id[which(annot$cell_cluster==j)])
            tmp = cor_mat[cells_i,cells_j]
            tab = data.frame("cluster_i" = rep(i,nrow(tmp)),
                             "cluster_j" = rep(j,nrow(tmp)),
                             "inter_corr" = rowMeans(tmp))
            inter_corr=rbind(inter_corr,tab)
        
    }
}

inter_corr = inter_corr %>% filter( (cluster_i == "C1" & cluster_j == "C1") |
                                        (cluster_i == "C2" & cluster_j == "C1") |
                                        (cluster_i == "C4" & cluster_j == "C1")) 
inter_corr = inter_corr %>% dplyr::mutate(association = factor(paste0(cluster_i,"_",cluster_j),
                                          levels=c("C4_C1","C2_C1","C1_C1")))

table(inter_corr$association)
inter_corr %>% group_by(association) %>% filter(inter_corr > 0.2) %>% summarise(n())

my_comparisons <- list(c("C4_C1", "C2_C1"),c("C2_C1", "C1_C1") )
png(file.path(plotdir_unsup,"Heatmaps", "Clusters_intercorrelation_violin.png"),
    height=1000,width=1000,res=300)
ggplot(inter_corr,aes(x=association, y=inter_corr, fill=association)) + 
    geom_violin(alpha=0.8) + theme_classic() +
    scale_fill_manual(values= unique(annot$cell_cluster_color)[c(4,2,1)]) +
    stat_summary(fun=median, geom="point", size=2, color="black") +
    theme(axis.text.x = element_text(angle=90)) +
    ggpubr::stat_compare_means(comparisons = my_comparisons,method = "t.test")
dev.off()

# Generate BAM for all clusters
inputBam = "/media/pacome/LaCie/InstitutCurie/Documents/Data/results/BAMs/MM468_K27/"
peakCall = FALSE
if(peakCall == TRUE) subset_bam_call_peaks(scExp_cf,odir = file.path(ChromSCape_directory, "Diff_Analysis_Gene_Sets"),
inputBam = list.files(file.path(inputBam),full.names =T ))

############################
## Supervised Analysis    ##
############################
load(file.path(ChromSCape_directory, "correlation_clustering",
               paste0(prefix,".RData")))
annotK27 <- rtracklayer::import(file.path(
    maindir,"annotation","MM468_peaks_K27.bed.gz"))
start(annotK27)  = start(annotK27) -1
annotK27$ID = paste0(seqnames(annotK27),"_",start(annotK27),"_",end(annotK27))
annot_10k = read.table(gzfile(file.path(maindir,"annotation",
                                          "gencode.v34.annotation.transcriptTSS_10k.bed.gz")),
                         sep="\t", header = F)
colnames(annot_10k) = c("chr","start","end","transcripts","gene","strand")
annot_10k = as(annot_10k,"GRanges")
annot_10k$ID = paste0(seqnames(annot_10k),"_",start(annot_10k),"_", end(annot_10k))
start(annot_10k) = start(annot_10k) + 4500
end(annot_10k) = start(annot_10k)  + 1000
width(annot_10k)

annotK27$gene = ""
hits <- findOverlaps(annotK27, annot_10k)
agg <- aggregate(annot_10k, hits, gene=paste(gene, collapse = ","))
annotK27$gene[match(subsetByOverlaps(annotK27,annot_10k)$name,annotK27$ID)] = agg$gene

#################################################################################
###  Differential Analysis
#################################################################################

#################################################################################
qval.th = 0.1
logFC.th =1
enrichment_qval = 0.1
#################################################################################

#################################################################################
########################## DA & GSA cluster sc DMSO C4 vs C3 ####################
#################################################################################

scExp_C4_vs_C2 = scExp_cf

scExp_C4_vs_C2 = scExp_C4_vs_C2[,which(scExp_C4_vs_C2$cell_cluster %in% c("C2","C4"))]
scExp_C4_vs_C2$cell_cluster[which(scExp_C4_vs_C2$cell_cluster == "C4")] = "C1"
table(scExp_C4_vs_C2$cell_cluster,scExp_C4_vs_C2$sample_id)
mat = NULL
for(cluster in c("C1","C2")){
    samps = unique(scExp_C4_vs_C2$sample_id[which(scExp_C4_vs_C2$cell_cluster == cluster)])
    tot = table(scExp_C4_vs_C2$sample_id[which(scExp_C4_vs_C2$cell_cluster == cluster)])
    for(sample in names(tot)){
        if(tot[sample]>50){
            col = as.matrix(rowSums(counts(scExp_C4_vs_C2)[,which(scExp_C4_vs_C2$cell_cluster == cluster &
                                                          scExp_C4_vs_C2$sample_id == sample)]))
            colnames(col) = paste0(sample,"_",cluster)
            if(is.null(mat)) mat = col else mat = cbind(mat,col)
        } 
    } 
}

myrefs <- list(
    DMSO_C4 = colnames(mat)[grep("C1",colnames(mat))]
)
mygps <- list(
    DMSO_C2 = colnames(mat)[grep("C2",colnames(mat))]
)

refs <- names(myrefs)
groups <- names(mygps)

#selection of peaks with at least a log2 RPKM of 1 in one sample
RZ_sel = mat
feature = as.data.frame(scExp_C4_vs_C2@rowRanges)
feature = feature[,c(1,2,3)]
annot = data.frame(sample_id = colnames(mat), cluster = c("C1","C1","C1","C2","C2"))
res <- geco.ChIPseqCompareLimma(mat=RZ_sel,
                                metadata =annot,
                                ref=myrefs,
                                groups=mygps,
                                featureTab=feature
)
res$Gene = annotK27$gene[match(res$id,annotK27$ID)]
under_res = res %>% dplyr::filter(log2FC.DMSO_C2 < -logFC.th & qval.DMSO_C2 < qval.th) %>% 
    dplyr::arrange(qval.DMSO_C2) %>% dplyr::select(id,seqnames,start,end,log2FC.DMSO_C2,qval.DMSO_C2,Gene)
over_res = res %>% dplyr::filter(log2FC.DMSO_C2 > logFC.th & qval.DMSO_C2 < qval.th) %>% 
    dplyr::arrange(qval.DMSO_C2) %>% dplyr::select(id,seqnames,start,end,log2FC.DMSO_C2,qval.DMSO_C2,Gene)

WriteXLS(c("over_res","under_res","res"), 
         ExcelFileName = file.path(ChromSCape_directory, "Diff_Analysis_Gene_Sets",
                                     paste0("DA_grouped_DMSO_C2_vs_DMSO_C4_peaks.xlsx")),
         SheetNames = c(paste0("Over_",logFC.th,"_",qval.th,"_n",nrow(over_res)),
                        paste0("Under_-",logFC.th,"_",qval.th,"_n",nrow(under_res)),"All"),
         perl = "perl", verbose = FALSE, row.names = FALSE, col.names = TRUE, AdjWidth = T,
         AutoFilter = TRUE, BoldHeaderRow = TRUE, na = "", FreezeRow = 1, FreezeCol = 1)

# Volcano plot
res = readxl::read_xlsx(file.path(ChromSCape_directory, "Diff_Analysis_Gene_Sets",paste0("DA_grouped_DMSO_C2_vs_DMSO_C4_peaks.xlsx")),sheet = 3)
res = as.data.frame(res)
# colnames(res) = c("seqnames", "start", "end", "id", "cdiff.C2", "pval.C2", "qval.C2","Gene")
res$qval.DMSO_C2 = as.numeric(res$qval.DMSO_C2)
res$pval.DMSO_C2 = as.numeric(res$pval.DMSO_C2)
summaryTab <- geco.summaryCompareedgeR(restab=res,
                                       ref=myrefs,
                                       groups=mygps,
                                       qval.th=qval.th,
                                       fc.th=logFC.th,
                                       plotdir=file.path(plotdir_sup))


# Gene Set Analysis GROUPED
run_GSA = TRUE
if(run_GSA){
    annotbase <- "MSigDB" 
    database <- MSIG.ls ##MSigDB
    Overexpressed  <- Underexpressed <- data.frame()
    data("hg38.GeneTSS")
    GencodeGene = hg38.GeneTSS 
    # possibleGenes <- unique(unlist(strsplit(annotK27$gene,split=",")))
    possibleGenes <- union( unique(unlist(strsplit(res$Gene,split=","))),
                            c(as.character(unique(GencodeGene$gene)),"ELFN2"))
    gp <- groups
    ref <- refs
    print(paste0("Processing ",gp, " vs ", ref, " _ ",annotbase ))
    
    signific <- which(res[,paste("qval",gp,sep=".")] <= qval.th & abs(res[,paste("log2FC",gp,sep=".")]) > logFC.th)
    over <- which(res[,paste("qval",gp,sep=".")] <= qval.th & res[,paste("log2FC",gp,sep=".")] > logFC.th)
    under <- which(res[,paste("qval",gp,sep=".")] <= qval.th & res[,paste("log2FC",gp,sep=".")] < -logFC.th)
    print(paste0("significant = ", length(signific))) ; print(paste0("over = ", length(over))) ; print(paste0("under = ", length(under)))
    
    if(length(over)){
        enrich.test <- geco.enrichmentTest(gene.sets=database,
                                           mylist = unique(unlist(strsplit(res$Gene[over],split=","))),
                                           possibleIds=possibleGenes)
        enrich.test <- data.frame(Gene_set_name=rownames(enrich.test), enrich.test, check.names=FALSE)
        enrich.test <- merge( subset(MSIG.gs, select=-Genes), enrich.test, by.x="Gene.Set", by.y="Gene_set_name", all.y=TRUE, sort=FALSE ) ## Get class of gene set
        enrich.test <- enrich.test[order(enrich.test$`p-value`),]
        Overexpressed  <- enrich.test
    }
    if(length(under)){
        enrich.test <- geco.enrichmentTest(gene.sets=database,
                                           mylist= unique(unlist(strsplit(res$Gene[under],split=","))),
                                           possibleIds=possibleGenes)
        enrich.test <- data.frame(Gene_set_name=rownames(enrich.test), enrich.test, check.names=FALSE)
        enrich.test <- merge( subset(MSIG.gs, select=-Genes), enrich.test, by.x="Gene.Set", by.y="Gene_set_name", all.y=TRUE, sort=FALSE ) ## Get class of gene set
        enrich.test <- enrich.test[order(enrich.test$`p-value`),]
        ind <- which(enrich.test$`q-value`<= 0.1);if(!length(ind)){ind <- 1:20}
        Underexpressed <- enrich.test[ind,]		
        }
    
    WriteXLS(
        c("Overexpressed", "Underexpressed"),
        ExcelFileName = file.path(ChromSCape_directory, "Diff_Analysis_Gene_Sets",
                                  paste0("Enrichment_test_",gp,"_vs_",ref,
                                         "_logFC",round(logFC.th,2),"_peaks.xlsx")), 
        SheetNames = c( paste0("Overexp_in_", gp), paste0("Underexp_in_", gp) ),
        perl = "perl", verbose = FALSE, row.names = FALSE, col.names = TRUE,
        AdjWidth = TRUE, AutoFilter = TRUE, BoldHeaderRow = TRUE, na = "",
        FreezeRow = 1, FreezeCol = 1)
}

#################################################################################
### DA & GSA grouped Persister (C1) vs DMSO (C2 + C4) with BULK - CONSENSUS PEAKS #######
#################################################################################

scExp_Pers_vs_DMSO = scExp_cf
scExp_Pers_vs_DMSO = scExp_Pers_vs_DMSO[,which(scExp_Pers_vs_DMSO$cell_cluster %in% c("C1","C2","C4"))]
scExp_Pers_vs_DMSO$cell_cluster[which(scExp_Pers_vs_DMSO$cell_cluster %in% c("C1"))] = "bis"
scExp_Pers_vs_DMSO$cell_cluster[which(scExp_Pers_vs_DMSO$cell_cluster %in% c("C2","C4"))] = "C1"
scExp_Pers_vs_DMSO$cell_cluster[which(scExp_Pers_vs_DMSO$cell_cluster %in% c("bis"))] = "C2"

# Grouped DA
mat = NULL
for(cluster in c("C1","C2")){
    samps = unique(scExp_Pers_vs_DMSO$sample_id[which(scExp_Pers_vs_DMSO$cell_cluster == cluster)])
    tot = table(scExp_Pers_vs_DMSO$sample_id[which(scExp_Pers_vs_DMSO$cell_cluster == cluster)])
    for(sample in names(tot)){
        if(tot[sample]>50){
            col = as.matrix(rowSums(counts(scExp_Pers_vs_DMSO)[,which(scExp_Pers_vs_DMSO$cell_cluster == cluster &
                                                                        scExp_Pers_vs_DMSO$sample_id == sample)]))
            colnames(col) = paste0(sample,"_",cluster)
            if(is.null(mat)) mat = col else mat = cbind(mat,col)
        } 
    } 
}

# Read bulk Count matrices
resultZ <- read.csv(gzfile(file.path(maindir,"input","bulk_ChIPseq","MM468","CountTable_bulk_MM468_H3K27me3_peaks.csv.gz")))
resultZ$ID  <- paste(resultZ$Chromosome,resultZ$Begin,resultZ$End,sep="_")
resultZ = resultZ[,-c(1,2,3)]
bulk_RZ = resultZ
colnames(bulk_RZ) [grep("DMSO",colnames(bulk_RZ))] =
    paste0(colnames(bulk_RZ) [grep("DMSO",colnames(bulk_RZ))],"_C1")
colnames(bulk_RZ) [grep("5FU",colnames(bulk_RZ))] =
    paste0(colnames(bulk_RZ) [grep("5FU",colnames(bulk_RZ))],"_C2")

bulk_RZ = bulk_RZ[match(rownames(mat),bulk_RZ$ID),]
rownames(bulk_RZ) = bulk_RZ$ID
bulk_RZ = bulk_RZ[,-ncol(bulk_RZ)]

mat = cbind(mat, bulk_RZ)
myrefs <- list(
    DMSO = colnames(mat)[grep("C1",colnames(mat))]
)
mygps <- list(
    Persister = colnames(mat)[grep("C2",colnames(mat))]
)

refs <- names(myrefs)
groups <- names(mygps)

#selection of peaks with at least a log2 RPKM of 1 in one sample
RZ_sel = mat
feature = as.data.frame(scExp_Pers_vs_DMSO@rowRanges)
feature = feature[,c(1,2,3)]
annot = data.frame(sample_id = colnames(mat), cluster = gsub(".*_","",colnames(mat)))
res <- geco.ChIPseqCompareLimma(mat=RZ_sel,
                                metadata =annot,
                                ref=myrefs,
                                groups=mygps,
                                featureTab=feature
)

res$Gene = annotK27$gene[match(res$id,gsub(":|-","_",annotK27$ID))]
under_res = res %>% dplyr::filter(log2FC.Persister < -logFC.th & qval.Persister < qval.th) %>% 
    dplyr::arrange(qval.Persister) %>% dplyr::select(id,seqnames,start,end,
                                                     log2FC.Persister,qval.Persister, Gene)
over_res = res %>% dplyr::filter(log2FC.Persister > logFC.th & qval.Persister < qval.th) %>% 
    dplyr::arrange(qval.Persister) %>% dplyr::select(id,seqnames,start,end,
                                                     log2FC.Persister,qval.Persister, Gene)

WriteXLS(c("over_res","under_res","res"), 
         ExcelFileName = file.path(ChromSCape_directory, "Diff_Analysis_Gene_Sets",
                                   paste0("DA_grouped_Persister_vs_DSMO_with_bulk_peaks.xlsx")),
         SheetNames = c(paste0("Over_",round(logFC.th,2),"_",qval.th,"_n",nrow(over_res)),
                        paste0("Under_-",round(logFC.th,2),"_",qval.th,"_n",nrow(under_res)),"All"), perl = "perl", verbose = FALSE, row.names = FALSE, col.names = TRUE, AdjWidth = T, AutoFilter = TRUE, BoldHeaderRow = TRUE, na = "", FreezeRow = 1, FreezeCol = 1)

# Volcano plot
res = readxl::read_xlsx(file.path(ChromSCape_directory, "Diff_Analysis_Gene_Sets",paste0("DA_grouped_Persister_vs_DSMO_with_bulk_peaks.xlsx")),sheet = 3)
res = as.data.frame(res)

res$qval.Persister = as.numeric(res$qval.Persister)
res$pval.Persister = as.numeric(res$pval.Persister)
summaryTab <- geco.summaryCompareedgeR(restab=res,
                                       ref=myrefs,
                                       groups=mygps,
                                       qval.th=qval.th,
                                       fc.th=logFC.th,
                                       plotdir=file.path(plotdir_sup))

# Anticorrelation with scRNA
scRNA = readxl::read_xlsx(file.path(maindir, "output","scRNAseq","MM468","Persister",
                                    "Supervised","Tables","Differential_analysis_Limma_logFC_1.58.xlsx"),sheet= 3)
res_unique_gene = res %>% tidyr::separate_rows(Gene) %>% unique %>% dplyr::select(id, log2FC.Persister,
                                                                                  qval.Persister,Gene)
scRNA_scChIP = left_join(scRNA, res_unique_gene,by =c("Symbol"="Gene"))
top_demethylated <- head(scRNA_scChIP$Symbol[order(scRNA_scChIP$log2FC.Persister,decreasing=F)],n=15)

pdf(file.path(plotdir_sup,"scRNA_scChIP_log2FC_Persister_sc_grouped_vs_DMSO_peaks.pdf"))
sp <- ggscatter(scRNA_scChIP, x = "log2FC.Persister", y = "log2FC.C2_pers",
                add = "reg.line", ylab = "scRNA log2FC", xlab="scChIP K27 log2FC",  # Add regressin line
                add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE,
                label="Symbol",repel=TRUE,
                label.select= top_demethylated) #add size of dot function of initial expression in DMSO
# Add correlation coefficient
print(sp + stat_cor(method = "pearson",label.x=1,label.sep = "\n"))

# With only diff chip seq
scRNA_scChIP_filt = scRNA_scChIP %>% dplyr::filter(abs(log2FC.Persister) > logFC.th,
                                                   qval.Persister < qval.th)
sp <- ggscatter(scRNA_scChIP_filt, x = "log2FC.Persister", y = "log2FC.C2_pers",
                add = "reg.line", ylab = "scRNA log2FC", xlab="scChIP K27 log2FC",  # Add regressin line
                add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE,
                label="Symbol",repel=TRUE,
                label.select= top_demethylated) #add size of dot function of initial expression in DMSO
# Add correlation coefficient
print(sp + stat_cor(method = "pearson",label.x=1,label.sep = "\n"))
dev.off()

# Gene Set Analysis GROUPED
run_GSA = TRUE
if(run_GSA){
    annotbase <- "MSigDB" 
    database <- MSIG.ls ##MSigDB
    Overexpressed  <- Underexpressed <- data.frame()
    data("hg38.GeneTSS")
    GencodeGene = hg38.GeneTSS 
    # possibleGenes <- unique(unlist(strsplit(annotK27$gene,split=",")))
    possibleGenes <- union( unique(unlist(strsplit(res$Gene,split=","))),
                            c(as.character(unique(GencodeGene$gene)),"ELFN2"))
    gp <- groups
    ref <- refs
    print(paste0("Processing ",gp, " vs ", ref, " _ ",annotbase ))
    
    signific <- which(res[,paste("qval",gp,sep=".")] <= qval.th & abs(res[,paste("log2FC",gp,sep=".")]) > logFC.th)
    over <- which(res[,paste("qval",gp,sep=".")] <= qval.th & res[,paste("log2FC",gp,sep=".")] > logFC.th)
    under <- which(res[,paste("qval",gp,sep=".")] <= qval.th & res[,paste("log2FC",gp,sep=".")] < -logFC.th)
    print(paste0("significant = ", length(signific))) ; print(paste0("over = ", length(over))) ; print(paste0("under = ", length(under)))
    
    if(length(over)){
        enrich.test <- geco.enrichmentTest(gene.sets=database,
                                           mylist = unique(unlist(strsplit(res$Gene[over],split=","))),
                                           possibleIds=possibleGenes)
        enrich.test <- data.frame(Gene_set_name=rownames(enrich.test), enrich.test, check.names=FALSE)
        enrich.test <- merge( subset(MSIG.gs, select=-Genes), enrich.test, by.x="Gene.Set", by.y="Gene_set_name", all.y=TRUE, sort=FALSE ) ## Get class of gene set
        enrich.test <- enrich.test[order(enrich.test$`p-value`),]
        Overexpressed  <- enrich.test
    }
    if(length(under)){
        enrich.test <- geco.enrichmentTest(gene.sets=database,
                                           mylist= unique(unlist(strsplit(res$Gene[under],split=","))),
                                           possibleIds=possibleGenes)
        enrich.test <- data.frame(Gene_set_name=rownames(enrich.test), enrich.test, check.names=FALSE)
        enrich.test <- merge( subset(MSIG.gs, select=-Genes), enrich.test, by.x="Gene.Set", by.y="Gene_set_name", all.y=TRUE, sort=FALSE ) ## Get class of gene set
        enrich.test <- enrich.test[order(enrich.test$`p-value`),]
        ind <- which(enrich.test$`q-value`<= 0.1);if(!length(ind)){ind <- 1:20}
        Underexpressed <- enrich.test[ind,]		
    }
    
    WriteXLS(
        c("Overexpressed", "Underexpressed"),
        ExcelFileName = file.path(ChromSCape_directory, "Diff_Analysis_Gene_Sets",
                                  paste0("Enrichment_test_",gp,"_vs_",ref,
                                         "_logFC",round(logFC.th,2),"_peaks.xlsx")), 
        SheetNames = c( paste0("Overexp_in_", gp), paste0("Underexp_in_", gp) ),
        perl = "perl", verbose = FALSE, row.names = FALSE, col.names = TRUE,
        AdjWidth = TRUE, AutoFilter = TRUE, BoldHeaderRow = TRUE, na = "",
        FreezeRow = 1, FreezeCol = 1)
}

#################################################################################
### DA & GSA grouped Persister (C1) vs DMSO (C2 + C4) with BULK - TSS 10K #######
#################################################################################

run10k = TRUE
if(run10k){
    # 10k mat 
    out <- import_scExp_gz(file.path(datadir,"Count_Matrices"),
                           pattern = "*_H3K27me3_TSS.tsv.gz")
    
    # Save raw
    datamatrix_10k = out$datamatrix
    annot_raw_10k = out$annot_raw
    scExp_10k = create_scExp(datamatrix_10k, annot_raw_10k)
    scExp_10k = exclude_features_scExp(scExp_10k, exclude_regions)
    
    datamatrix_10k = counts(scExp_10k)
    mat = NULL
    for(cluster in c("C1","C2")){
        samps = unique(scExp_Pers_vs_DMSO$sample_id[which(scExp_Pers_vs_DMSO$cell_cluster == cluster)])
        tot = table(scExp_Pers_vs_DMSO$sample_id[which(scExp_Pers_vs_DMSO$cell_cluster == cluster)])
        for(sample in names(tot)){
            cells = scExp_Pers_vs_DMSO$cell_id[which(scExp_Pers_vs_DMSO$cell_cluster == cluster &
                                                         scExp_Pers_vs_DMSO$sample_id == sample)]
            cells = intersect(cells,colnames(datamatrix_10k))
            
            if(length(cells)>50){
                col = as.matrix(rowSums(datamatrix_10k[,match(cells, colnames(datamatrix_10k))]))
                print(paste0(sample, " - ", length(cells)))
                colnames(col) = paste0(sample,"_",cluster)
                if(is.null(mat)) mat = col else mat = cbind(mat,col)
            } 
        } 
    }
    
    # Read bulk Count matrices
    resultZ <- read.csv(gzfile(file.path(maindir,"input","bulk_ChIPseq","MM468","CountTable_bulk_MM468_H3K27me3_TSS.csv.gz")))
    row = paste0(resultZ$Chromosome, ":", resultZ$Begin, "-",resultZ$End)
    dups = duplicated(row)
    resultZ = resultZ[!dups,]
    rownames(resultZ) <- row[!dups]
    feature = as.data.frame(resultZ[,c(1,2,3)])
    rownames(feature) = gsub(":|-","_",rownames(feature))
    resultZ = resultZ[,-c(1,2,3)]
    
    bulk_RZ = resultZ
    rownames(bulk_RZ) = gsub(":|-","_",rownames(bulk_RZ))
    colnames(bulk_RZ) [grep("DMSO",colnames(bulk_RZ))] =
        paste0(colnames(bulk_RZ) [grep("DMSO",colnames(bulk_RZ))],"_bulk_C1")
    colnames(bulk_RZ) [grep("5FU",colnames(bulk_RZ))] =
        paste0(colnames(bulk_RZ) [grep("5FU",colnames(bulk_RZ))],"_bulk_C2")
    
    bulk_RZ = bulk_RZ[match(rownames(mat),rownames(bulk_RZ)),]
    all.equal(rownames(bulk_RZ), rownames(mat))
    mat = cbind(mat, bulk_RZ)
    
    myrefs <- list(
        DMSO = colnames(mat)[grep("C1",colnames(mat))]
    )
    mygps <- list(
        Persister = colnames(mat)[grep("C2",colnames(mat))]
    )
    
    refs <- names(myrefs)
    groups <- names(mygps)
    
    RZ_sel = mat
    annot = data.frame(sample_id = colnames(mat), cluster = gsub(".*_","",colnames(mat)))
    rownames(annot) = annot$sample_id
    res <- geco.ChIPseqCompareLimma(mat = RZ_sel,
                                    metadata = annot,
                                    ref = myrefs,
                                    groups = mygps,
                                    featureTab = feature
    )
    
    res$Gene = annot_10k$gene[match(res$id,annot_10k$ID)]
    under_res = res %>% dplyr::filter(log2FC.Persister < -logFC.th & qval.Persister < qval.th) %>% 
        dplyr::arrange(qval.Persister) %>% dplyr::select(id,Chromosome,Begin,End,
                                                         log2FC.Persister,qval.Persister, Gene)
    over_res = res %>% dplyr::filter(log2FC.Persister > logFC.th & qval.Persister < qval.th) %>% 
        dplyr::arrange(qval.Persister) %>% dplyr::select(id,Chromosome,Begin,End,
                                                         log2FC.Persister,qval.Persister, Gene)
    
    WriteXLS(c("over_res","under_res","res"), 
             ExcelFileName = file.path(ChromSCape_directory, "Diff_Analysis_Gene_Sets",
                                       paste0("DA_grouped_Persister_vs_DSMO_with_bulk_TSS.xlsx")),
             SheetNames = c(paste0("Over_",round(logFC.th,2),"_",qval.th,"_n",nrow(over_res)),
                            paste0("Under_-",round(logFC.th,2),"_",qval.th,"_n",nrow(under_res)),"All"), perl = "perl", verbose = FALSE, row.names = FALSE, col.names = TRUE, AdjWidth = T, AutoFilter = TRUE, BoldHeaderRow = TRUE, na = "", FreezeRow = 1, FreezeCol = 1)
    
    # Volcano plot
    res = readxl::read_xlsx(file.path(ChromSCape_directory, "Diff_Analysis_Gene_Sets",
                                      paste0("DA_grouped_Persister_vs_DSMO_with_bulk_TSS.xlsx")),sheet = 3)
    res = as.data.frame(res)
    res$qval.Persister = as.numeric(res$qval.Persister)
    res$pval.Persister = as.numeric(res$pval.Persister)
    summaryTab <- geco.summaryCompareedgeR(restab=res,
                                           ref=myrefs,
                                           groups=mygps,
                                           qval.th=qval.th,
                                           fc.th=logFC.th,
                                           plotdir=file.path(plotdir_sup))
    
    # Anticorrelation with scRNA
    scRNA = readxl::read_xlsx(file.path(maindir, "output","scRNAseq","MM468","Persister",
                                        "Supervised","Tables","Differential_analysis_Limma_logFC_1.58.xlsx"),sheet= 3)
    res_unique_gene = res %>% group_by(Gene) %>%
        slice_max(abs(log2FC.Persister)) %>% dplyr::select(id, log2FC.Persister,qval.Persister,Gene)
    scRNA_scChIP = left_join(scRNA, res_unique_gene,by =c("Symbol"="Gene"))
    
    top_demethylated <- scRNA_scChIP %>% filter(log2FC.C2_pers > log2(3), scRNA_scChIP$log2FC.Persister < -1) %>% 
        arrange(desc(log2FC.Persister)) %>% dplyr::select(Symbol) %>% unique %>%  head(n=35) 
    top_demethylated = top_demethylated$Symbol
    
    options(ggrepel.max.overlaps = 25)
    options(max.overlaps = 25)
    pdf(file.path(plotdir_sup,"scRNA_scChIP_log2FC_Persister_sc_grouped_vs_DMSO_TSS.pdf"))
    sp <- ggscatter(scRNA_scChIP, x = "log2FC.Persister", y = "log2FC.C2_pers",
                    add = "reg.line", ylab = "scRNA log2FC", xlab="scChIP K27 log2FC",  # Add regressin line
                    add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
                    conf.int = TRUE,
                    label="Symbol",repel=TRUE,
                    label.select= top_demethylated) #add size of dot function of initial expression in DMSO
    # Add correlation coefficient
    print(sp + stat_cor(method = "pearson",label.x=1,label.sep = "\n"))
    
    # With only diff chip seq
    scRNA_scChIP_filt = scRNA_scChIP %>% dplyr::filter(abs(log2FC.Persister) > logFC.th,
                                                       qval.Persister < qval.th)
    sp <- ggscatter(scRNA_scChIP_filt, x = "log2FC.Persister", y = "log2FC.C2_pers",
                    add = "reg.line", ylab = "scRNA log2FC", xlab="scChIP K27 log2FC",  # Add regressin line
                    add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
                    conf.int = TRUE,
                    label="Symbol",repel=TRUE,
                    label.select= top_demethylated, max.overlaps = 25) #add size of dot function of initial expression in DMSO
    # Add correlation coefficient
    print(sp + stat_cor(method = "pearson",label.x=1,label.sep = "\n"))
    dev.off()
    
    run_GSA = TRUE
    if(run_GSA){
        annotbase <- "MSigDB" 
        database <- MSIG.ls ##MSigDB
        Overexpressed  <- Underexpressed <- data.frame()
        data("hg38.GeneTSS")
        GencodeGene = hg38.GeneTSS 
        # possibleGenes <- unique(unlist(strsplit(annotK27$gene,split=",")))
        possibleGenes <- unique(union(res$Gene,c(as.character(unique(GencodeGene$gene)),"ELFN2")))
        gp <- groups
        ref <- refs
        print(paste0("Processing ",gp, " vs ", refs, " _ ",annotbase ))
        
        signific <- which(res[,paste("qval",gp,sep=".")] <= qval.th & abs(res[,paste("log2FC",gp,sep=".")]) > logFC.th)
        over <- which(res[,paste("qval",gp,sep=".")] <= qval.th & res[,paste("log2FC",gp,sep=".")] > logFC.th)
        under <- which(res[,paste("qval",gp,sep=".")] <= qval.th & res[,paste("log2FC",gp,sep=".")] < -logFC.th)
        print(paste0("significant = ", length(signific))) ; print(paste0("over = ", length(over))) ; print(paste0("under = ", length(under)))
        
        if(length(over)){
            enrich.test <- geco.enrichmentTest(gene.sets=database,
                                               mylist = unique(unlist(strsplit(res$Gene[over],split=","))),
                                               possibleIds=possibleGenes)
            enrich.test <- data.frame(Gene_set_name=rownames(enrich.test), enrich.test, check.names=FALSE)
            enrich.test <- merge( subset(MSIG.gs, select=-Genes), enrich.test, by.x="Gene.Set", by.y="Gene_set_name", all.y=TRUE, sort=FALSE ) ## Get class of gene set
            enrich.test <- enrich.test[order(enrich.test$`p-value`),]
            Overexpressed  <- enrich.test
        }
        if(length(under)){
            enrich.test <- geco.enrichmentTest(gene.sets=database,
                                               mylist= unique(unlist(strsplit(res$Gene[under],split=","))),
                                               possibleIds=possibleGenes)
            enrich.test <- data.frame(Gene_set_name=rownames(enrich.test), enrich.test, check.names=FALSE)
            enrich.test <- merge( subset(MSIG.gs, select=-Genes), enrich.test, by.x="Gene.Set", by.y="Gene_set_name", all.y=TRUE, sort=FALSE ) ## Get class of gene set
            enrich.test <- enrich.test[order(enrich.test$`p-value`),]
            ind <- which(enrich.test$`q-value`<= 0.1);if(!length(ind)){ind <- 1:20}
            Underexpressed <- enrich.test[ind,]		
        }
        
        WriteXLS(
            c("Overexpressed", "Underexpressed"),
            ExcelFileName = file.path(ChromSCape_directory, "Diff_Analysis_Gene_Sets",
                                      paste0("Enrichment_test_",gp,"_vs_",refs,
                                             "_logFC",round(logFC.th,2),"_TSS.xlsx")), 
            SheetNames = c( paste0("Overexp_in_", gp), paste0("Underexp_in_", gp) ),
            perl = "perl", verbose = FALSE, row.names = FALSE, col.names = TRUE,
            AdjWidth = TRUE, AutoFilter = TRUE, BoldHeaderRow = TRUE, na = "",
            FreezeRow = 1, FreezeCol = 1)
    }
    
}

################################################################################
################ Enrichment in TSS #############################################
################################################################################

# Load bulk + sc grouped peak DA 
res = readxl::read_xlsx(file.path(ChromSCape_directory, "Diff_Analysis_Gene_Sets",
                                  "DA_grouped_Persister_vs_DSMO_with_bulk_peaks.xlsx"), sheet = 3)
# Load annotation peak affectation 
annotK27_peak_affectation = read.table( file.path(maindir, "annotation", "annotK27_peak_affectation.tsv"), 
                                        sep = "\t", header = TRUE)
res$peak_affectation = annotK27_peak_affectation$peak_affectation[match(res$id,annotK27_peak_affectation$ID)]
col = data.frame("col"=c("#999999ff","#dededeff","#3b9ab2ff","#78b7c5ff","#ebcc2aff","#e1af00ff"))
sig_diff <- (abs(res$log2FC.Persister)> logFC.th & res$qval.Persister < qval.th)
sig_under <- (res$log2FC.Persister< -logFC.th & res$qval.Persister < qval.th)
sig_over <- (res$log2FC.Persister> logFC.th & res$qval.Persister < qval.th)
l = list("Diff" = sig_diff, "Under" = sig_under, "Over" = sig_over)

#Enrichment
fisher <- function(a,b,c,d){
    data <- matrix(c(a,b,c,d),ncol=2)
    c(p = fisher.test(data)$p.value)
}

sig = l[["Over"]]
sum_sig = sum(table(res$peak_affectation[sig]))
sum_all = sum(table(res$peak_affectation))
res$peak_affectation = factor(res$peak_affectation, levels =c(
    "intergenic", "intergenic_enhancer", "tss_pc", "genebody_pc", "tss_other", "genebody_other"
))
tab = res[sig,] %>% dplyr::count(peak_affectation,.drop=F) %>%
    left_join(res %>% dplyr::count(peak_affectation,.drop=F), by = "peak_affectation") %>%
    mutate(sum_sig =sum_sig, sum_all =sum_all, Freq = log2( (n.x/sum_sig) / (n.y/sum_all) ),
           color = col$col) %>% rowwise() %>% mutate(p=fisher(n.x,sum_sig-n.x,n.y,sum_all-n.y)) %>% 
    mutate(is_significative = ifelse(p <0.05, TRUE,FALSE),
           Freq = ifelse(is.infinite(Freq),0,Freq),
           fill = ifelse(is_significative,color,"white"))

png(file.path(plotdir,paste0("enrichment_HQ_enriched.png")),width = 1000,height = 1000,res=300)
p1 = tab %>% ggplot() + geom_bar(aes(x=peak_affectation, y=Freq),fill= tab$fill, color = tab$color,
                                 stat="identity") + ggtitle(paste0(i," enriched H3K27me3 peaks")) +
    theme_classic() + ggtitle(paste0("Enriched H3K27me3 peaks - ",length(which(sig))," peaks")) + 
    ylab("Log2 Enrichment") + xlab("")  + geom_hline(yintercept = 0, color ="grey")
print(p1 +  theme(axis.text.x = element_text(angle = 90, vjust =1), legend.position = "None"))
dev.off()

sig = l[["Under"]]
sum_sig = sum(table(res$peak_affectation[sig]))
sum_all = sum(table(res$peak_affectation))

tab = res[sig,] %>% dplyr::count(peak_affectation,.drop=F) %>%
    left_join(res %>% dplyr::count(peak_affectation,.drop=F), by = "peak_affectation") %>%
    mutate(sum_sig =sum_sig, sum_all =sum_all, Freq = log2( (n.x/sum_sig) / (n.y/sum_all) ),
           color = col$col) %>% rowwise() %>% mutate(p=fisher(n.x,sum_sig-n.x,n.y,sum_all-n.y)) %>% 
    mutate(is_significative = ifelse(p <0.05, TRUE,FALSE),
           Freq = ifelse(is.infinite(Freq),0,Freq),
           fill = ifelse(is_significative,color,"white"))

png(file.path(plotdir,paste0("enrichment_HQ_under.png")),width = 1000,height = 1000,res=300)
p2= tab %>% ggplot() + geom_bar(aes(x=peak_affectation, y=Freq),fill= tab$fill, color = tab$color,
                                stat="identity") + ggtitle(paste0(i," depleted H3K27me3 peaks")) +
    theme_classic() + ggtitle(paste0("Depleted H3K27me3 peaks - ",length(which(sig))," peaks")) +
    ylab("Log2 Enrichment") + xlab("")  + geom_hline(yintercept = 0, color ="grey")
print(p2 + theme(axis.text.x = element_text(angle = 90, vjust =1), legend.position = "None"))
dev.off()

png(file.path(plotdir,paste0("enrichment_HQ_over_wolabel.png")),width = 800,height = 1000,res=300)
print(p1 + theme(text=element_blank()))
dev.off()

png(file.path(plotdir,paste0("enrichment_HQ_under_wolabel.png")),width = 800,height = 1000,res=300)
print(p2 + theme(text=element_blank()))
dev.off()

