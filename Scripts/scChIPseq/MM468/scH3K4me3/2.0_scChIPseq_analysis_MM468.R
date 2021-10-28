library(here)
source(file.path(here(),"Scripts","functions.R"))
source(file.path(here(),"Scripts","global_var.R"))

# Directories -------------------------------------------------------------
maindir = here()
resdir <- file.path(maindir,"output","scChIPseq","MM468","scH3K4me3"); if(!file.exists(resdir)){dir.create(resdir)}
QCdir <- file.path(resdir,"QC"); if(!file.exists(QCdir)){dir.create(QCdir)}
input_dir  <- file.path(maindir,"input","scChIPseq","MM468")

# Unsupervised directories
resdir <- file.path(resdir, "Unsupervised") ; if(!file.exists(resdir)){dir.create(resdir)}
resdir_UMAP <- file.path(resdir,"UMAP") ; if(!file.exists(resdir_UMAP)){dir.create(resdir_UMAP)}
resdir_heatmaps = file.path(resdir,"Heatmaps"); if(!dir.exists(resdir_heatmaps)) dir.create(resdir_heatmaps)
RDatadir <- file.path(resdir,"RData") ; if(!file.exists(RDatadir)){dir.create(RDatadir)}

dataset_name = "MM468_H3K4me3_10k_TSS"
ref_genome ="hg38"

metadata = data.frame(sample_id = c(
  "MM468_DMSO1_day0", "MM468_5FU1_day60"),
  sample_id_color = c("#afafafff", "#118675ff"))

# Run ChromSCape :
set.seed(47)

################
## Input Data ##
################
out <- import_scExp(list.files(file.path(input_dir,"Count_Matrices"), full.names = TRUE, pattern = "H3K4me3_TSS"),
                       remove_pattern =  "_H3K4me3_TSS")

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
min_reads_per_cell = 1600
max_quantile_read_per_cell = 95
n_features_to_keep = 5000

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
num_cell_after_QC_filt_scExp(scExp = scExp, annot = annot_raw, datamatrix = datamatrix)


# Excluding known CNA loci
exclude_regions = rtracklayer::import(file.path(maindir,"annotation","MM468_identified_CNA.bed.gz"))
scExp = exclude_features_scExp(scExp, exclude_regions, by = "region")

# Normalize 
scExp = normalize_scExp(scExp, type = "CPM")

# Annotate with genes
scExp = feature_annotation_scExp(scExp, ref = ref_genome)

# Reduce dimensions
scExp = reduce_dims_scExp(scExp, n = 50, dimension_reductions = c("PCA","UMAP"),
    batch_correction = F, verbose = F)

# Color samples
color_df = metadata
scExp = colors_scExp(scExp, c("sample_id","total_counts","batch_id"),color_df = color_df) # add colouring

# Plot PCA + UMAP sample / total counts
# Development version of plot_reduced_dim_scExp function
set.seed(47)
png(file.path(resdir_UMAP,"PCA_sample_HQ.png"), res=300,height=1200,width=2000)
plot_reduced_dim_scExp(scExp,"sample_id","PCA",downsample = 5000, transparency = 0.55)
dev.off()

png(file.path(resdir_UMAP,"UMAP_sample_no_cf_TFIDF_HQ.png"), res=300,height=1200,width=2000)
plot_reduced_dim_scExp(scExp,"sample_id","UMAP",downsample = 5000, transparency = 0.45, size = 2)
dev.off()

# Clear shots
png(file.path(resdir_UMAP,"PCA_sample_HQ_clear.png"), res=300,height=1200,width=1200)
plot_reduced_dim_scExp(scExp,"sample_id","PCA",downsample = 5000, transparency = 0.55) +
    theme(legend.position = "None", text = element_blank())
dev.off()

png(file.path(resdir_UMAP,"UMAP_sample_no_cf_TFIDF_HQ_clear.png"), res=300,height=1200,width=1200)
plot_reduced_dim_scExp(scExp,"sample_id","UMAP",downsample = 5000, transparency = 0.45, size = 2) +
    theme(legend.position = "None", text = element_blank())
dev.off()

set.seed(47)
png(file.path(resdir_UMAP,"PCA_counts_HQ.png"), res=300,height=1200,width=2000)
plot_reduced_dim_scExp(scExp,"total_counts","PCA",downsample = 5000, transparency = 0.55)
dev.off()

png(file.path(resdir_UMAP,"UMAP_counts_HQ.png"), res=300,height=1200,width=2000)
plot_reduced_dim_scExp(scExp,"total_counts","UMAP",downsample = 5000, transparency = 0.55)
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
run_correlation_filtering = FALSE
if(run_correlation_filtering){
  scExp_cf = filter_correlated_cell_scExp(scExp,
                                          random_iter = 5,
                                          corr_threshold = 99,
                                          percent_correlation = 1)
  gc()
} else {
  scExp_cf = scExp
}

# Run consensus clustering (RAM heavy) - optional
runConsensus <- TRUE
if(runConsensus){
    scExp_cf = consensus_clustering_scExp(scExp_cf, reps = 10,
                                          maxK = 10,
                                          clusterAlg = "hc",
                                          prefix = file.path(resdir_heatmaps, prefix))
    
    # Plot item consensus score
    plot_cluster_consensus_scExp(scExp_cf)
    
}

# Choosing 2 clusters
nclust = 2
scExp_cf = choose_cluster_scExp(scExp_cf, nclust = nclust, consensus = runConsensus)

# Coloring clusters
scExp_cf = colors_scExp(scExp_cf, annotCol = "cell_cluster",
                     color_by = c("cell_cluster")) # add colouring


qs::qsave(scExp_cf, file = file.path(RDatadir, paste0(prefix,"_correlation.qs"))) #save the data 

# Plotting
png(file.path(resdir_heatmaps,"ClusteringHeatMap.png"), res=300,height=4000,width=4000)
plot_heatmap_scExp(scExp_cf)
dev.off()

#plot UMAP cluster
png(file.path(resdir_UMAP,"UMAP_cluster_HQ.png"), res=300,height=1200,width=2000)
plot_reduced_dim_scExp_devel(scExp_cf,"cell_cluster","UMAP",downsample = 5000, transparency = 0.35, size = 2)
dev.off()

png(file.path(resdir_UMAP,"UMAP_cluster_HQ_clear.png"), res=300,height=1200,width=1200)
plot_reduced_dim_scExp_devel(scExp_cf,"cell_cluster","UMAP",downsample = 5000, transparency = 0.35, size = 2) +
    theme(legend.position = "None", text = element_blank())
dev.off()

png(file.path(resdir_UMAP,"UMAP_sample_id_no_cf_HQ.png"), res=300,height=1200,width=2000)
plot_reduced_dim_scExp_devel(scExp_cf,"sample_id","UMAP",downsample = 5000, transparency = 0.35, size = 2)
dev.off()

png(file.path(resdir_UMAP,"UMAP_sample_id_no_cf_HQ_clear.png"), res=300,height=1200,width=1200)
plot_reduced_dim_scExp_devel(scExp_cf,"sample_id","UMAP",downsample = 5000, transparency = 0.35, size = 2) +
  theme(legend.position = "None", text = element_blank())
dev.off()


