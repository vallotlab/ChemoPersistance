library(here)
source(file.path(here(),"Scripts","functions.R"))
source(file.path(here(),"Scripts","global_var.R"))

maindir = here()

dataset_name = "MM468_H3K4me3_10k_TSS"
ref_genome ="hg38"

metadata = data.frame(sample_id = c(
  "MM468_DMSO1_day0", "MM468_5FU1_day60"),
  sample_id_color = c("#afafafff", "#1AB8AD"))


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

out <- import_scExp_gz(file.path(datadir,"Count_Matrices"),
                       pattern = "*_H3K4me3_TSS.tsv.gz")

# Save raw
datamatrix = out$datamatrix
annot_raw = out$annot_raw
save(datamatrix,annot_raw, file = file.path(ChromSCape_directory,"scChIP_raw.RData")) 

scExp <- create_scExp(datamatrix, annot_raw)

#################
## Filter Data ##
#################

# Quality filtering
min_reads_per_cell = 1600
max_quantile_read_per_cell = 95
min_percent_to_keep_feature = 1 # Keep all peaks as peaks were pre-defined (not genomic bins!)

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
num_cell_after_QC_filt_scExp(scExp = scExp, annot = annot_raw)


# Excluding known CNA loci
exclude_regions = rtracklayer::import(file.path(maindir,"annotation","MM468_identified_CNA.bed.gz"))
exclude_regions = as.data.frame(exclude_regions)
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
png(file.path(plotdir_unsup,"UMAPS","PCA_sample_HQ.png"), res=300,height=1200,width=2000)
plot_reduced_dim_scExp_devel(scExp,"sample_id","PCA",downsample = 5000, transparency = 0.55)
dev.off()

png(file.path(plotdir_unsup,"UMAPS","UMAP_sample_HQ.png"), res=300,height=1200,width=2000)
plot_reduced_dim_scExp_devel(scExp,"sample_id","UMAP",downsample = 5000, transparency = 0.55)
dev.off()

# Clear shots
png(file.path(plotdir_unsup,"UMAPS","PCA_sample_HQ_clear.png"), res=300,height=1200,width=1200)
plot_reduced_dim_scExp_devel(scExp,"sample_id","PCA",downsample = 5000, transparency = 0.55) +
    theme(legend.position = "None", text = element_blank())
dev.off()

png(file.path(plotdir_unsup,"UMAPS","UMAP_sample_HQ_clear.png"), res=300,height=1200,width=1200)
plot_reduced_dim_scExp_devel(scExp,"sample_id","UMAP",downsample = 5000, transparency = 0.55) +
    theme(legend.position = "None", text = element_blank())
dev.off()

set.seed(47)
png(file.path(plotdir_unsup,"UMAPS","PCA_counts_HQ.png"), res=300,height=1200,width=2000)
plot_reduced_dim_scExp_devel(scExp,"total_counts","PCA",downsample = 5000, transparency = 0.55)
dev.off()

png(file.path(plotdir_unsup,"UMAPS","UMAP_counts_HQ.png"), res=300,height=1200,width=2000)
plot_reduced_dim_scExp_devel(scExp,"total_counts","UMAP",downsample = 5000, transparency = 0.55)
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
    scExp_cf = consensus_clustering_scExp(scExp_cf, reps = 10,
                                          maxK = 10,
                                          clusterAlg = "hc",
                                          prefix = file.path(cor_plot_dir, prefix))
    
    # Plot item consensus score
    plot_cluster_consensus_scExp(scExp_cf)
    
}

# Choosing 7 clusters
nclust = 2
scExp_cf = choose_cluster_scExp(scExp_cf, nclust = nclust, consensus = runConsensus)

# Coloring clusters
scExp_cf = colors_scExp(scExp_cf, annotCol = "cell_cluster",
                     color_by = c("cell_cluster")) # add colouring


save(scExp_cf, file = file.path(ChromSCape_directory, "correlation_clustering",
                               paste0(prefix,".RData"))) #save the data 

# Plotting
png(file.path(plotdir_unsup,"Heatmaps","ClusteringHeatMap.png"), res=300,height=4000,width=4000)
plot_heatmap_scExp(scExp_cf)
dev.off()

#plot UMAP cluster
png(file.path(plotdir_unsup,"UMAPS","UMAP_cluster_HQ.png"), res=300,height=1200,width=2000)
plot_reduced_dim_scExp_devel(scExp_cf,"cell_cluster","UMAP",downsample = 5000, transparency = 0.35, size = 2)
dev.off()

png(file.path(plotdir_unsup,"UMAPS","UMAP_cluster_HQ_clear.png"), res=300,height=1200,width=1200)
plot_reduced_dim_scExp_devel(scExp_cf,"cell_cluster","UMAP",downsample = 5000, transparency = 0.35, size = 2) +
    theme(legend.position = "None", text = element_blank())
dev.off()

# Differential analysis
qval.th = 0.1
logFC.th = 1
enrichment_qval = 0.1

# Supervised
scExp_Pers_vs_DMSO = scExp_cf
scExp_Pers_vs_DMSO$cell_cluster[which(scExp_Pers_vs_DMSO$sample_id == "MM468_5FU1_day60")] = "C2"
scExp_Pers_vs_DMSO$cell_cluster[which(scExp_Pers_vs_DMSO$sample_id != "MM468_5FU1_day60")] = "C1"
de_type  = "one_vs_rest"
scExp_Pers_vs_DMSO = differential_analysis_scExp(scExp_Pers_vs_DMSO,de_type = de_type)

### Gene Set Analysis
scExp_Pers_vs_DMSO = gene_set_enrichment_analysis_scExp(scExp_Pers_vs_DMSO,enrichment_qval = 0.1,
                                                      qval.th = qval.th, cdiff.th = logFC.th)

diff_Pers_vs_DMSO = scExp_Pers_vs_DMSO@metadata$diff$res
diff_Pers_vs_DMSO = diff_Pers_vs_DMSO[,c(1:4,10:14)]
colnames(diff_Pers_vs_DMSO) = gsub("C2","Pers_vs_DMSO",colnames(diff_Pers_vs_DMSO))

GSA_Pers_vs_DMSO = scExp_Pers_vs_DMSO@metadata$enr

save(diff_Pers_vs_DMSO,GSA_Pers_vs_DMSO,
     file = file.path(ChromSCape_directory, "Diff_Analysis_Gene_Sets",
                      paste0("Pers_vs_DMSO_",qval.th,"_",logFC.th,"_",de_type,".RData"))) #save the data 

annot_10k = read.table(gzfile(file.path(maindir,"annotation",
                                        "gencode.v34.annotation.transcriptTSS_10k.bed.gz")),
                       sep="\t", header = F)
colnames(annot_10k) = c("chr","start","end","transcripts","gene","strand")

diff_Pers_vs_DMSO$gene = 
  annot_10k$gene[match(rownames(diff_Pers_vs_DMSO),paste0(annot_10k$chr,"_",annot_10k$start,"_",annot_10k$end))]

under_res = diff_Pers_vs_DMSO %>% dplyr::filter(cdiff.Pers_vs_DMSO < -logFC.th & qval.Pers_vs_DMSO < qval.th) %>% 
  dplyr::arrange(qval.Pers_vs_DMSO) %>% dplyr::select(chr,start,end,,cdiff.Pers_vs_DMSO,qval.Pers_vs_DMSO,gene)
over_res = diff_Pers_vs_DMSO %>% dplyr::filter(cdiff.Pers_vs_DMSO > logFC.th & qval.Pers_vs_DMSO < qval.th) %>% 
  dplyr::arrange(qval.Pers_vs_DMSO) %>% dplyr::select(chr,start,end,cdiff.Pers_vs_DMSO,qval.Pers_vs_DMSO,gene)

WriteXLS(c("over_res","under_res","diff_Pers_vs_DMSO"),
         ExcelFileName = file.path(
           ChromSCape_directory,"Diff_Analysis_Gene_Sets",paste0("DA_Pers_vs_DMSO_logFC_",
                                                                 round(logFC.th,2),".xlsx")),
         SheetNames = c(
           paste0("Over_",round(logFC.th,2),"_",qval.th,"_n",
                  nrow(over_res)),paste0(
                    "Under_-",round(logFC.th,2),"_",qval.th,
                    "_n",nrow(under_res)),"All"),
         perl = "perl", verbose = FALSE, row.names = FALSE,
         col.names = TRUE, AdjWidth = T, AutoFilter = TRUE,
         BoldHeaderRow = TRUE, na = "", FreezeRow = 1, FreezeCol = 1)

Overexpressed = GSA_Pers_vs_DMSO$Overexpressed[[2]]
Underexpressed = GSA_Pers_vs_DMSO$Underexpressed[[2]]

WriteXLS(
  c("Overexpressed", "Underexpressed"),
  ExcelFileName = file.path(ChromSCape_directory, "Diff_Analysis_Gene_Sets",
                            paste0("Enrichment_test_Pers_vs_DMSO_logFC",
                                   round(logFC.th,2),"_TSS.xlsx")), 
  SheetNames = c( paste0("Overexp_in_Pers"), paste0("Underexp_in_Pers") ),
  perl = "perl", verbose = FALSE, row.names = FALSE, col.names = TRUE,
  AdjWidth = TRUE, AutoFilter = TRUE, BoldHeaderRow = TRUE, na = "",
  FreezeRow = 1, FreezeCol = 1)
