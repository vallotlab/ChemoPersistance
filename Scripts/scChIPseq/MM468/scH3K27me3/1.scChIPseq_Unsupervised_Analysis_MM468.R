library(here)
source(file.path(here(),"Scripts","functions.R"))
source(file.path(here(),"Scripts","global_var.R"))

# Directories -------------------------------------------------------------
maindir = here()
resdir <- file.path(maindir,"output","scChIPseq","MM468","scH3K27me3"); if(!file.exists(resdir)){dir.create(resdir)}
QCdir <- file.path(resdir,"QC"); if(!file.exists(QCdir)){dir.create(QCdir)}
input_dir  <- file.path(maindir,"input","scChIPseq","MM468")

# Unsupervised directories
resdir <- file.path(resdir, "Unsupervised") ; if(!file.exists(resdir)){dir.create(resdir)}
resdir_UMAP <- file.path(resdir,"UMAP") ; if(!file.exists(resdir_UMAP)){dir.create(resdir_UMAP)}
resdir_heatmaps = file.path(resdir,"Heatmaps"); if(!dir.exists(resdir_heatmaps)) dir.create(resdir_heatmaps)
RDatadir <- file.path(resdir,"RData") ; if(!file.exists(RDatadir)){dir.create(RDatadir)}

dataset_name = "MM468_H3K27me3_peaks"
ref_genome ="hg38"

metadata = data.frame(sample_id = c(
    "MM468_DMSO1_day60", "MM468_DMSO3_day77", "MM468_DMSO5_day131",
    "MM468_5FU1_day33", "MM468_5FU2_day67",
    "MM468_5FU6_day131", "MM468_5FU3_day147", "MM468_5FU2_day171"),
    sample_id_color = c("#dfdfdfff", "#999999ff","#363636",
                        "#118675ff", "#8cc453ff",
                        "#ff5722ff", "#feb40fff", "#fd8508ff"))

# Run ChromSCape :
set.seed(47)

################
## Input Data ##
################

# Reading in matrices
out <- import_scExp(
    grep("GSKJ4", invert = TRUE, value = TRUE,
         list.files(file.path(input_dir,"Count_Matrices"), full.names = TRUE, pattern = "H3K27me3_peaks")),
    remove_pattern = "_H3K27me3_peaks")

# Save raw
datamatrix = out$datamatrix
annot_raw = out$annot_raw
qs::qsave(datamatrix,file = file.path(RDatadir,"datamatrix.qs"))
qs::qsave(annot_raw, file = file.path(RDatadir,"annot_raw.qs"))

scExp <- create_scExp(datamatrix, annot_raw)

#################
## Filter Data ##
#################

# Quality filtering
min_reads_per_cell = 3000
max_quantile_read_per_cell = 95
n_features_to_keep = 10000 # Keep all peaks as peaks were pre-defined (not genomic bins!)

prefix <- paste0(dataset_name,"_",min_reads_per_cell,"_",
                 n_features_to_keep,"_",
                 max_quantile_read_per_cell,"_uncorrected")

upper_threshold = quantile(Matrix::colSums(counts(scExp)),probs = max_quantile_read_per_cell/100 )

# Plot distribution before filtering
plot_distribution_scExp(scExp, raw = TRUE, log10 = TRUE, pseudo_counts = 1,bins = 150) +
    geom_vline(xintercept = log10(min_reads_per_cell),color = "green") +
    geom_vline(xintercept = log10(upper_threshold),color = "red")

# Filter
scExp = filter_scExp(scExp, min_cov_cell = min_reads_per_cell,
                     quant_removal = max_quantile_read_per_cell) 
scExp = find_top_features(scExp, n = n_features_to_keep) 

# Plot distribution after filtering
plot_distribution_scExp(scExp, raw = T, log10 = T, pseudo_counts = 1,bins = 150)  +
    geom_vline(xintercept = log10(min_reads_per_cell),color = "green") +
    geom_vline(xintercept = log10(upper_threshold),color = "red")

# Number of cells passing filtering
num_cell_after_QC_filt_scExp(scExp, annot = annot_raw, datamatrix = datamatrix)

# Excluding known CNA loci
exclude_regions = rtracklayer::import(file.path(maindir,"annotation","MM468_identified_CNA.bed.gz"))
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


qs::qsave(scExp, file = file.path(RDatadir, paste0(prefix,".qs")))

# Plot UMAPS and PCA
png(file.path(resdir_UMAP,"PCA_sample_HQ.png"), res=300,height=1200,width=2000)
print(plot_reduced_dim_scExp(scExp, color_by = "sample_id", reduced_dim = "PCA",
                             downsample = 5000, transparency = 0.35))
dev.off()

png(file.path(resdir_UMAP,"UMAP_sample_HQ.png"), res=300,height=1200,width=2000)
print(plot_reduced_dim_scExp(scExp,"sample_id","UMAP",downsample = 5000, transparency = 0.35))
dev.off()

# Clear shots
png(file.path(resdir_UMAP,"PCA_sample_HQ_clear.png"), res=300,height=1200,width=1200)
print(plot_reduced_dim_scExp(scExp,"sample_id","PCA",downsample = 5000, transparency = 0.35) +
          theme(legend.position = "None", text = element_blank()))
dev.off()

png(file.path(resdir_UMAP,"UMAP_sample_HQ_clear.png"), res=300,height=1200,width=1200)
print(plot_reduced_dim_scExp(scExp,"sample_id","UMAP",downsample = 5000, transparency = 0.35) +
          theme(legend.position = "None", text = element_blank()))
dev.off()

png(file.path(resdir_UMAP,"PCA_counts_HQ.png"), res=300,height=1200,width=2000)
print(plot_reduced_dim_scExp(scExp,"total_counts","PCA",downsample = 5000, transparency = 0.35))
dev.off()

png(file.path(resdir_UMAP,"UMAP_counts_HQ.png"), res=300,height=1200,width=2000)
print(plot_reduced_dim_scExp(scExp,"total_counts","UMAP",downsample = 5000, transparency = 0.35))
dev.off()

# Calculate correlation matrix and hierarchical clustering
scExp = correlation_and_hierarchical_clust_scExp(scExp)


qs::qsave(scExp, file = file.path(RDatadir, paste0(prefix,"_correlation.qs")))
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
                                          plot_consclust = "pdf", 
                                          prefix = file.path(cor_plot_dir, prefix))
    
    # Plot item consensus score
    png(file.path(resdir_heatmaps,"cluster_consensus_score.png"), res=300,height=1200,width=2000)
    plot_cluster_consensus_scExp(scExp_cf)
    dev.off()
    
}

# Choosing 4 clusters
nclust = 4
scExp_cf = choose_cluster_scExp(scExp_cf, nclust = nclust, consensus = runConsensus)

# Coloring clusters
color_df = data.frame(cell_cluster = unique(scExp_cf$cell_cluster))
color_df$cell_cluster_color = c("#ffd700ff", "#ffb5c5ff", "#3cb371ff" ,"#668b8bff")
                                # "#473c8bff", "#3cb371ff", "#94524aff", "#ba55d3ff")
scExp_cf = colors_scExp(scExp_cf, annotCol = "cell_cluster",
                        color_by = c("cell_cluster"),color_df = color_df) # add colouring


qs::qsave(scExp_cf, file = file.path(RDatadir, paste0(prefix,"_correlation_filtered.qs"))) #save the data 

# Plotting
png(file.path(resdir_heatmaps,"ClusteringHeatMap.png"), res=300,height=4000,width=4000)
print(plot_heatmap_scExp(scExp_cf))
dev.off()

# Plotting
png(file.path(resdir_heatmaps,"Consensus_score.png"), res=300,height=4000,width=4000)
print( plot_cluster_consensus_scExp(scExp_cf))
dev.off()

#plot UMAP cluster
png(file.path(resdir_UMAP,"UMAP_cluster_HQ.png"), res=300,height=1200,width=2000)
print(plot_reduced_dim_scExp_devel(scExp_cf,"cell_cluster","UMAP",downsample = 5000, transparency = 0.35, size = 2))
dev.off()

png(file.path(resdir_UMAP,"UMAP_cluster_HQ_clear.png"), res=300,height=1200,width=1200)
print(plot_reduced_dim_scExp_devel(scExp_cf,"cell_cluster","UMAP",downsample = 5000, transparency = 0.35, size = 2) +
          theme(legend.position = "None", text = element_blank()))
dev.off()

#plot UMAP cluster
png(file.path(resdir_UMAP,"UMAP_sample_cf_HQ.png"), res=300,height=1200,width=2000)
print(plot_reduced_dim_scExp_devel(scExp_cf,"sample_id","UMAP",downsample = 5000, transparency = 0.35, size = 2))
dev.off()

png(file.path(resdir_UMAP,"UMAP_sample_cf_HQ_clear.png"), res=300,height=1200,width=1200)
print(plot_reduced_dim_scExp_devel(scExp_cf,"sample_id","UMAP",downsample = 5000, transparency = 0.35, size = 2) +
          theme(legend.position = "None", text = element_blank()))
dev.off()

############################
## Intra Correlation      ##
############################

#Table association
num_cell_in_cluster_scExp(scExp_cf) 
tab = table(scExp_cf$sample_id, scExp_cf$cell_cluster)
write.csv(tab, file.path(resdir_heatmaps, "Contingency_sample_cluster_bargraph.csv"))

# Plot contingency bar graph - scChIP seq
annot = SingleCellExperiment::colData(scExp_cf)

annot = as_tibble(annot) %>% group_by(sample_id,sample_id_color,cell_cluster,cell_cluster_color) %>%
    summarise(count=n())

annot = left_join(annot,
                  annot %>% ungroup %>% group_by(sample_id) %>% summarise(sum = sum(count)),
                  by ="sample_id")
annot = annot %>% mutate(freq = count/sum)

annot = annot %>% mutate(sample_id=factor(sample_id,levels=metadata$sample_id)) 

png(file.path(resdir_heatmaps, "Contingency_sample_cluster_bargraph.png"),width=800,height=1500,res=300)
print(ggplot(annot) + geom_bar(aes(y=freq,x=sample_id, fill=cell_cluster), stat="identity") +
          scale_fill_manual(values=color_df$cell_cluster_color) +
          theme_classic() + 
          theme(axis.text.x = element_text(angle = 90)) )
dev.off()

png(file.path(resdir_heatmaps, "Contingency_sample_cluster_bargraph_wolabel.png"),width=1500,height=1500,res=300)
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
    
    png(file.path(resdir_heatmaps, paste0("Pearson_Correlation_Inferno_scores_",Group,".png")), height=1150,width=1000,res=300)
    MAT <- mat.so.cor[,colnames(mat.so.cor) %in% annot$cell_id[annot$cell_cluster==Group]]
    debut <- round(100-(1-min(cor(MAT)))/0.02,0)
    image(cor(MAT),col=(viridis::inferno(100)[debut:100]))
    dev.off()
}
gc()

WriteXLS(IntraCorr_test,file.path(resdir_heatmaps,"IntraCorr_test.xls"),row.names = T)

png(file.path(resdir_heatmaps,"ExpressionGroup_intracorrelation_violin.png"),height=1000,width=1000,res=300)
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

png(file.path(resdir_heatmaps, "Sample_id_intracorrelation_violin.png"),
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
    samples_i = as.character(annot$sample_id[match(cells_i,annot$cell_id)])
    for(j in unique(annot$cell_cluster)){
            cells_j = as.character(annot$cell_id[which(annot$cell_cluster==j)])
            tmp = cor_mat[cells_i,cells_j]
            tab = data.frame("cells_i" = cells_i,
                             "sample_i" = samples_i,
                             "cluster_i" = rep(i,nrow(tmp)),
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
correlated_samp = (inter_corr %>% filter(association == "C2_C1" & inter_corr > 0.2))$sample_i
table(correlated_samp)

my_comparisons <- list(c("C4_C1", "C2_C1"),c("C2_C1", "C1_C1") )
png(file.path(resdir_heatmaps, "Clusters_intercorrelation_violin.png"),
    height=1000,width=1000,res=300)
ggplot(inter_corr,aes(x=association, y=inter_corr, fill=association)) + 
    geom_violin(alpha=0.8) + theme_classic() +
    scale_fill_manual(values= unique(annot$cell_cluster_color)[c(4,2,1)]) +
    stat_summary(fun=median, geom="point", size=2, color="black") +
    theme(axis.text.x = element_text(angle=90)) +
    ggpubr::stat_compare_means(comparisons = my_comparisons,method = "t.test")
dev.off()

