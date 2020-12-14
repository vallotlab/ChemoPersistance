library(here)
source(file.path(here(),"Scripts","functions.R"))
source(file.path(here(),"Scripts","global_var.R"))

maindir = here()

dataset_name = "MM468_H3K27me3_peaks"
ref_genome ="hg38"

metadata = data.frame(sample_id = c(
    "MM468_DMSO6_day60", "MM468_DMSO3_day77", "MM468_DMSO1_day131",
    "MM468_5FU6_day33", "MM468_5FU5_day67", "MM468_GSKJ4_day91","MM468_GSKJ4_5FU3_day91",
    "MM468_5FU4_day131", "MM468_5FU3_day147", "MM468_5FU5_day171"),
    sample_id_color = c("#dfdfdfff", "#999999ff","#363636",
                        "#118675ff", "#55D42F","#8E3899", "#E676FA",
                        "#ffea2eff", "#ff9800ff", "#ff5722ff"))

datadir = file.path(maindir, "input", "scChIPseq", "MM468")
outdir = file.path(maindir, "output", "scChIPseq", "MM468",dataset_name)
if(!dir.exists(outdir)) dir.create(outdir)
plotdir = file.path(outdir, "Plots"); if(!dir.exists(plotdir)) dir.create(plotdir)
plotdir_unsup = file.path(plotdir, "Unsupervised"); if(!dir.exists(plotdir_unsup)) dir.create(plotdir_unsup)
dir.create(plotdir_unsup,"Heatmaps",showWarnings = F)
dir.create(plotdir_unsup,"UMAPS",showWarnings = F)
plotdir_sup = file.path(plotdir, "Supervised"); if(!dir.exists(plotdir_sup)) dir.create(plotdir_sup)
dir.create(plotdir_sup,"Heatmaps",showWarnings = F)
dir.create(plotdir_sup,"UMAPS",showWarnings = F)

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
out <- import_scExp(unzip(file.path(datadir,"Count_Matrices", "MM468_K27_peaks_1000.zip"),
                          exdir = tempdir()))

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

subsample_n = 500
if(length(subsample_n)>0){
    print("Doing subsampling")
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

sel = scExp$cell_id[which(!scExp$sample_id %in% c("MM468_GSKJ4_day91","MM468_GSKJ4_5FU3_day91"))]
scExp = scExp[,sel]

# Development version of plot_reduced_dim_scExp function
set.seed(47)
png(file.path(plotdir_unsup,"UMAPs","PCA_sample_HQ.png"), res=300,height=1200,width=2000)
plot_reduced_dim_scExp_devel(scExp,"sample_id","PCA",downsample = 5000, transparency = 0.35)
dev.off()

png(file.path(plotdir_unsup,"UMAPs","UMAP_sample_HQ.png"), res=300,height=1200,width=2000)
plot_reduced_dim_scExp_devel(scExp,"sample_id","UMAP",downsample = 5000, transparency = 0.35)
dev.off()

# Clear shots
png(file.path(plotdir_unsup,"UMAPs","PCA_sample_HQ_clear.png"), res=300,height=1200,width=1200)
plot_reduced_dim_scExp_devel(scExp,"sample_id","PCA",downsample = 5000, transparency = 0.35) +
    theme(legend.position = "None", text = element_blank())
dev.off()

png(file.path(plotdir_unsup,"UMAPs","UMAP_sample_HQ_clear.png"), res=300,height=1200,width=1200)
plot_reduced_dim_scExp_devel(scExp,"sample_id","UMAP",downsample = 5000, transparency = 0.35) +
    theme(legend.position = "None", text = element_blank())
dev.off()

set.seed(47)
png(file.path(plotdir_unsup,"UMAPs","PCA_counts_HQ.png"), res=300,height=1200,width=2000)
plot_reduced_dim_scExp_devel(scExp,"total_counts","PCA",downsample = 5000, transparency = 0.35)
dev.off()

png(file.path(plotdir_unsup,"UMAPs","UMAP_counts_HQ.png"), res=300,height=1200,width=2000)
plot_reduced_dim_scExp_devel(scExp,"total_counts","UMAP",downsample = 5000, transparency = 0.35)
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
nclust = 5
scExp_cf = choose_cluster_scExp(scExp_cf, nclust = nclust, consensus = runConsensus)

# Coloring clusters
color_df = data.frame(cell_cluster = unique(scExp_cf$cell_cluster))
color_df$cell_cluster_color = c("#473c8bff", "#ffd700ff", "#3cb371ff" ,"#668b8bff","#ffb5c5ff")
# "#3cb371ff", "#94524aff", "#ba55d3ff")
scExp_cf = colors_scExp(scExp_cf, annotCol = "cell_cluster",
                        color_by = c("cell_cluster"),color_df = color_df) # add colouring


save(scExp_cf, file = file.path(ChromSCape_directory, "correlation_clustering",
                                paste0(prefix,".RData"))) #save the data 

# Plotting
png(file.path(plotdir_unsup,"Heatmaps","ClusteringHeatMap.png"), res=300,height=4000,width=4000)
plot_heatmap_scExp(scExp_cf)
dev.off()

#plot UMAP cluster
png(file.path(plotdir_unsup,"UMAPs","UMAP_cluster_HQ.png"), res=300,height=1200,width=2000)
plot_reduced_dim_scExp_devel(scExp_cf,"cell_cluster","UMAP",downsample = 5000, transparency = 0.35, size = 2)
dev.off()

png(file.path(plotdir_unsup,"UMAPs","UMAP_cluster_HQ_clear.png"), res=300,height=1200,width=1200)
plot_reduced_dim_scExp_devel(scExp_cf,"cell_cluster","UMAP",downsample = 5000, transparency = 0.35, size = 2) +
    theme(legend.position = "None", text = element_blank())
dev.off()

#plot UMAP cluster
png(file.path(plotdir_unsup,"UMAPs","UMAP_sample_cf_HQ.png"), res=300,height=1200,width=2000)
plot_reduced_dim_scExp_devel(scExp_cf,"sample_id","UMAP",downsample = 5000, transparency = 0.35, size = 2)
dev.off()

png(file.path(plotdir_unsup,"UMAPs","UMAP_sample_cf_HQ_clear.png"), res=300,height=1200,width=1200)
plot_reduced_dim_scExp_devel(scExp_cf,"sample_id","UMAP",downsample = 5000, transparency = 0.35, size = 2) +
    theme(legend.position = "None", text = element_blank())
dev.off()

############################
## Intra Correlation      ##
############################

load(file = file.path(ChromSCape_directory, "correlation_clustering",
                                     paste0(prefix,".RData")))
#Table association
num_cell_in_cluster_scExp(scExp_cf) 
tab = table(scExp_cf$sample_id, scExp_cf$cell_cluster)

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
ggplot(annot) + geom_bar(aes(y=freq,x=sample_id,fill=cell_cluster), stat="identity") +
    scale_fill_manual(values=color_df$cell_cluster_color) +
    theme_classic() + 
    theme(axis.text.x = element_text(angle = 90)) 
dev.off()

png(file.path(plotdir_unsup,"Heatmaps", "Contingency_sample_cluster_bargraph_wolabel.png"),width=1500,height=1500,res=300)
ggplot(annot) + geom_bar(aes(y=freq,x=sample_id,fill=cell_cluster), stat="identity") +
    scale_fill_manual(values=color_df$cell_cluster_color) +
    theme_classic() + 
    theme(text = element_blank()) 
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

inter_corr = inter_corr %>% filter( (cluster_i == "C2" & cluster_j == "C2") |
                                        (cluster_i == "C4" & cluster_j == "C2") |
                                        (cluster_i == "C3" & cluster_j == "C2")) 
inter_corr = inter_corr %>% dplyr::mutate(association = factor(paste0(cluster_i,"_",cluster_j),
                                          levels=c("C3_C2","C4_C2","C2_C2")))

my_comparisons <- list(c("C3_C2", "C4_C2"),c("C4_C2", "C2_C2") )
png(file.path(plotdir_unsup,"Heatmaps", "Clusters_intercorrelation_violin.png"),
    height=1000,width=1000,res=300)
ggplot(inter_corr,aes(x=association, y=inter_corr, fill=association)) + 
    geom_violin(alpha=0.8) + theme_classic() +
    scale_fill_manual(values= unique(annot$cell_cluster_color)[c(3,4,2)]) +
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

###  Differential Analysis
qval.th = 0.1
logFC.th =1
enrichment_qval = 0.1

#################################################################################
########################## DA & GSA cluster sc DMSO C4 vs C3 ####################
#################################################################################

scExp_C4_vs_C3 = scExp_cf

scExp_C4_vs_C3 = scExp_C4_vs_C3[,which(scExp_C4_vs_C3$cell_cluster %in% c("C3","C4"))]
scExp_C4_vs_C3$cell_cluster[which(scExp_C4_vs_C3$cell_cluster == "C3")] = "C1"
scExp_C4_vs_C3$cell_cluster[which(scExp_C4_vs_C3$cell_cluster == "C4")] = "C2"
table(scExp_C4_vs_C3$cell_cluster,scExp_C4_vs_C3$sample_id)
mat = NULL
for(cluster in c("C1","C2")){
    samps = unique(scExp_C4_vs_C3$sample_id[which(scExp_C4_vs_C3$cell_cluster == cluster)])
    tot = table(scExp_C4_vs_C3$sample_id[which(scExp_C4_vs_C3$cell_cluster == cluster)])
    for(sample in names(tot)){
        if(tot[sample]>50){
            col = as.matrix(rowSums(counts(scExp_C4_vs_C3)[,which(scExp_C4_vs_C3$cell_cluster == cluster &
                                                          scExp_C4_vs_C3$sample_id == sample)]))
            colnames(col) = paste0(sample,"_",cluster)
            if(is.null(mat)) mat = col else mat = cbind(mat,col)
        } 
    } 
}


library(ChromSCape)
myrefs <- list(
    DMSO_C3 = colnames(mat)[grep("C1",colnames(mat))]
)
mygps <- list(
    DMSO_C4 = colnames(mat)[grep("C2",colnames(mat))]
)

refs <- names(myrefs)
groups <- names(mygps)

#selection of peaks with at least a log2 RPKM of 1 in one sample
RZ_sel = mat
feature = as.data.frame(scExp_C4_vs_C3@rowRanges)
feature = feature[,c(1,2,3)]
annot = data.frame(sample_id = colnames(mat), cluster = c("C1","C1","C1","C2","C2"))
res <- geco.ChIPseqCompareLimma(mat=RZ_sel,
                                metadata =annot,
                                ref=myrefs,
                                groups=mygps,
                                featureTab=feature
)
res$Gene = annotK27$gene[match(res$id,annotK27$ID)]
under_res = res %>% dplyr::filter(log2FC.DMSO_C4 < -logFC.th & qval.DMSO_C4 < qval.th) %>% 
    dplyr::arrange(qval.DMSO_C4) %>% dplyr::select(id,seqnames,start,end,log2FC.DMSO_C4,qval.DMSO_C4,Gene)
over_res = res %>% dplyr::filter(log2FC.DMSO_C4 > logFC.th & qval.DMSO_C4 < qval.th) %>% 
    dplyr::arrange(qval.DMSO_C4) %>% dplyr::select(id,seqnames,start,end,log2FC.DMSO_C4,qval.DMSO_C4,Gene)

WriteXLS(c("over_res","under_res","res"), 
         ExcelFileName = file.path(ChromSCape_directory, "Diff_Analysis_Gene_Sets",
                                     paste0("DA_grouped_C4_vs_C3_.xlsx")),
         SheetNames = c(paste0("Over_",logFC.th,"_",qval.th,"_n",nrow(over_res)),
                        paste0("Under_-",logFC.th,"_",qval.th,"_n",nrow(under_res)),"All"),
         perl = "perl", verbose = FALSE, row.names = FALSE, col.names = TRUE, AdjWidth = T,
         AutoFilter = TRUE, BoldHeaderRow = TRUE, na = "", FreezeRow = 1, FreezeCol = 1)

# Volcano plot
res = readxl::read_xlsx(file.path(ChromSCape_directory, "Diff_Analysis_Gene_Sets",paste0("DA_grouped_C4_vs_C3_.xlsx")),sheet = 3)
res = as.data.frame(res)
# colnames(res) = c("seqnames", "start", "end", "id", "cdiff.C2", "pval.C2", "qval.C2","Gene")
res$qval.DMSO_C4 = as.numeric(res$qval.DMSO_C4)
res$pval.DMSO_C4 = as.numeric(res$pval.DMSO_C4)
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
                                         "_logFC",round(logFC.th,2),".xlsx")), 
        SheetNames = c( paste0("Overexp_in_", gp), paste0("Underexp_in_", gp) ),
        perl = "perl", verbose = FALSE, row.names = FALSE, col.names = TRUE,
        AdjWidth = TRUE, AutoFilter = TRUE, BoldHeaderRow = TRUE, na = "",
        FreezeRow = 1, FreezeCol = 1)
}


run_DA_sc = FALSE
if(run_DA_sc) {
    de_type  = "one_vs_rest"
    scExp_C4_vs_C3 = differential_analysis_scExp(scExp_C4_vs_C3,de_type = de_type)
    
    ### Gene Set Analysis
    scExp_C4_vs_C3 = gene_set_enrichment_analysis_scExp(scExp_C4_vs_C3,enrichment_qval = 0.1,
                                                        qval.th = qval.th, logFC.th = logFC.th)
    
    diff_C4_vs_C3 = scExp_C4_vs_C3@metadata$diff$res
    diff_C4_vs_C3 = diff_C4_vs_C3[,c(1:4,10:14)]
    colnames(diff_C4_vs_C3) = gsub("C2","C4_vs_C3",colnames(diff_C4_vs_C3))
    
    GSA_C4_vs_C3 = scExp_C4_vs_C3@metadata$enr
    
    diff_C4_vs_C3$gene = annotK27$gene[match(diff_C4_vs_C3$ID,gsub(":|-","_",annotK27$ID))]
    
    save(diff_C4_vs_C3,GSA_C4_vs_C3,
         file = file.path(ChromSCape_directory, "Diff_Analysis_Gene_Sets",
                          paste0("C4_vs_C3_",qval.th,"_",logFC.th,"_",de_type,".RData"))) #save the data 
    
    under_res = diff_C4_vs_C3 %>% dplyr::filter(cdiff.C4_vs_C3 < -logFC.th & qval.C4_vs_C3 < qval.th) %>% 
        dplyr::arrange(qval.C4_vs_C3) %>% dplyr::select(chr,start,end,,cdiff.C4_vs_C3,qval.C4_vs_C3,gene)
    over_res = diff_C4_vs_C3 %>% dplyr::filter(cdiff.C4_vs_C3 > logFC.th & qval.C4_vs_C3 < qval.th) %>% 
        dplyr::arrange(qval.C4_vs_C3) %>% dplyr::select(chr,start,end,cdiff.C4_vs_C3,qval.C4_vs_C3,gene)
    
    WriteXLS(c("over_res","under_res","diff_C4_vs_C3"),
             ExcelFileName = file.path(
                 ChromSCape_directory,"Diff_Analysis_Gene_Sets",paste0("DA_C4_vs_C3_logFC_",
                                                                       round(logFC.th,2),".xlsx")),
             SheetNames = c(
                 paste0("Over_",round(logFC.th,2),"_",qval.th,"_n",
                        nrow(over_res)),paste0(
                            "Under_-",round(logFC.th,2),"_",qval.th,
                            "_n",nrow(under_res)),"All"),
             perl = "perl", verbose = FALSE, row.names = FALSE,
             col.names = TRUE, AdjWidth = T, AutoFilter = TRUE,
             BoldHeaderRow = TRUE, na = "", FreezeRow = 1, FreezeCol = 1)
    
    Overexpressed = GSA_C4_vs_C3$Overexpressed
    Underexpressed = GSA_C4_vs_C3$Underexpressed
    
    WriteXLS(
        c("", ""),
        ExcelFileName = file.path(ChromSCape_directory,
                                  "Diff_Analysis_Gene_Sets",paste0("GSA_C4_vs_C3_logFC_",
                                                                   round(logFC.th,2),".xlsx")), 
        SheetNames = c( paste0("Overexp_in_", gp), paste0("Underexp_in_", gp) ),
        perl = "perl", verbose = FALSE, row.names = FALSE, col.names = TRUE,
        AdjWidth = TRUE, AutoFilter = TRUE, BoldHeaderRow = TRUE, na = "",
        FreezeRow = 1, FreezeCol = 1)
}


#################################################################################
########################## DA & GSA cluster scPersister vs C3 ###################
#################################################################################

scExp_Pers_vs_C3 = scExp_cf
scExp_Pers_vs_C3 = scExp_Pers_vs_C3[,which(scExp_Pers_vs_C3$cell_cluster %in% c("C3","C2"))]
scExp_Pers_vs_C3$cell_cluster[which(scExp_Pers_vs_C3$cell_cluster == "C3")] = "C1"


# Grouped DA
mat = NULL
for(cluster in c("C1","C2")){
    samps = unique(scExp_Pers_vs_C3$sample_id[which(scExp_Pers_vs_C3$cell_cluster == cluster)])
    tot = table(scExp_Pers_vs_C3$sample_id[which(scExp_Pers_vs_C3$cell_cluster == cluster)])
    for(sample in names(tot)){
        if(tot[sample]>50){
            col = as.matrix(rowSums(counts(scExp_Pers_vs_C3)[,which(scExp_Pers_vs_C3$cell_cluster == cluster &
                                                                      scExp_Pers_vs_C3$sample_id == sample)]))
            colnames(col) = paste0(sample,"_",cluster)
            if(is.null(mat)) mat = col else mat = cbind(mat,col)
        } 
    } 
}

myrefs <- list(
    DMSO_C3 = colnames(mat)[grep("C1",colnames(mat))]
)
mygps <- list(
    Persister = colnames(mat)[grep("C2",colnames(mat))]
)

refs <- names(myrefs)
groups <- names(mygps)

#selection of peaks with at least a log2 RPKM of 1 in one sample
RZ_sel = mat
feature = as.data.frame(scExp_Pers_vs_C3@rowRanges)
feature = feature[,c(1,2,3)]
annot = data.frame(sample_id = colnames(mat), cluster = c("C1","C1","C1","C2","C2","C2","C2"))
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
                                   paste0("DA_grouped_Persister_vs_C3.xlsx")),
         SheetNames = c(paste0("Over_",round(logFC.th,2),"_",qval.th,"_n",nrow(over_res)),
                        paste0("Under_-",round(logFC.th,2),"_",qval.th,"_n",nrow(under_res)),"All"), perl = "perl", verbose = FALSE, row.names = FALSE, col.names = TRUE, AdjWidth = T, AutoFilter = TRUE, BoldHeaderRow = TRUE, na = "", FreezeRow = 1, FreezeCol = 1)

# Volcano plot
res = readxl::read_xlsx(file.path(ChromSCape_directory, "Diff_Analysis_Gene_Sets",paste0("DA_grouped_Persister_vs_C3.xlsx")),sheet = 3)
res = as.data.frame(res)
# colnames(res) = c("seqnames", "start", "end", "id", "cdiff.C2", "pval.C2", "qval.C2","Gene")
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
scRNA$pval.C2_pers =as.numeric(scRNA$pval.C2_pers)
scRNA$qval.C2_pers =as.numeric(scRNA$qval.C2_pers)
res_unique_gene = res %>% tidyr::separate_rows(Gene) %>% unique %>% dplyr::select(id, log2FC.Persister,
                                                                           qval.Persister,Gene)
scRNA_scChIP = left_join(scRNA, res_unique_gene,by =c("Symbol"="Gene"))
top_demethylated <- head(scRNA_scChIP$Symbol[order(scRNA_scChIP$log2FC.Persister,decreasing=F)],n=15)

pdf(file.path(plotdir_sup,"scRNA_scChIP_log2FC_Persister_sc_grouped_vs_C3.pdf"))
sp <- ggscatter(scRNA_scChIP, x = "log2FC.Persister", y = "log2FC.C2_pers",
                add = "reg.line", ylab = "scRNA log2FC", xlab="scChIP K27 log2FC",  # Add regressin line
                add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE,
                label="Symbol",repel=TRUE,
                label.select= top_demethylated) #add size of dot function of initial expression in DMSO
# Add correlation coefficient
print(sp + stat_cor(method = "pearson",label.x=1,label.sep = "\n"))

# Only th
scRNA_scChIP_filt = scRNA_scChIP %>% dplyr::filter(abs(log2FC.Persister) > logFC.th,qval.Persister < qval.th)

sp <- ggscatter(scRNA_scChIP_filt, x = "log2FC.Persister", y = "log2FC.C2_pers",
                add = "reg.line", ylab = "scRNA log2FC", xlab="scChIP K27 log2FC",  # Add regressin line
                add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE,
                label="Symbol",repel=TRUE,
                label.select= top_demethylated) #add size of dot function of initial expression in DMSO
print(sp + stat_cor(method = "pearson",label.x=1,label.sep = "\n"))
                                                   
library(VennDiagram)
persisters = scRNA_scChIP$Symbol[which(scRNA_scChIP$log2FC.C2_pers > log2(3) &
                                           scRNA_scChIP$qval.C2_pers < 0.01)]
K27_dep = scRNA_scChIP$Symbol[which(scRNA_scChIP$log2FC.Persister < - logFC.th &
                                        scRNA_scChIP$qval.Persister < qval.th)]
VennDiagram_2(persisters,
              K27_dep,
              savePNG=F, "", "Venn", c("scRNA over", "K27 under"))
print(intersect(persisters, K27_dep))
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
                                         "_logFC",round(logFC.th,2),".xlsx")), 
        SheetNames = c( paste0("Overexp_in_", gp), paste0("Underexp_in_", gp) ),
        perl = "perl", verbose = FALSE, row.names = FALSE, col.names = TRUE,
        AdjWidth = TRUE, AutoFilter = TRUE, BoldHeaderRow = TRUE, na = "",
        FreezeRow = 1, FreezeCol = 1)
}

run_DA_sc = FALSE
if(run_DA_sc) {
    de_type  = "one_vs_rest"
    scExp_Pers_vs_C3 = differential_analysis_scExp(scExp_Pers_vs_C3,de_type = de_type)
    
    scDiff_grouped = readxl::read_xlsx(file.path(ChromSCape_directory, "Diff_Analysis_Gene_Sets",
                                                 paste0("DA_grouped_Persister_vs_C3.xlsx")),sheet = 3)
    scDiff_grouped$end = as.numeric(scDiff_grouped$end)
    colnames(scDiff_grouped)[1] = "chr"
    colnames(scDiff_grouped) = gsub("Persister","Persister_grouped",colnames(scDiff_grouped))
    scDiff_sc = readxl::read_xlsx(file.path(ChromSCape_directory, "Diff_Analysis_Gene_Sets",
                                            paste0("DA_Pers_vs_C3_logFC_1.xlsx")),sheet = 3)
    
    colnames(scDiff_sc) = gsub("Pers_vs_C3","Persister_notgrouped",colnames(scDiff_sc))
    scDiff_sc$end = as.numeric(scDiff_sc$end)
    
    bulkDiff = readxl::read_xlsx(file.path(maindir,"output","bulk_ChIPseq","MM468","K27_peaks_K27",
                                           "Supervised","Tables","Differential_analysis_Limma.xlsx"), sheet = 3)
    combined_table = left_join(scDiff_sc,scDiff_grouped)
    bulkDiff$end = as.numeric(bulkDiff$end)
    combined_table = left_join(combined_table,bulkDiff)
    colnames(combined_table) = gsub("X5FU2_3_5","bulk_Pers_vs_DMSO",colnames(combined_table))
    combined_table = combined_table %>% dplyr::select(-peak_affectation,-gene_affectation,-K4,-bivalent,-enhancer,-distance,-distance_all,
                                                      -transcript,-transcript_all, -genebody, -genebody_all, -gene)
    combined_table$cdiff.Persister_notgrouped = as.numeric(combined_table$cdiff.Persister_notgrouped)
    combined_table$qval.Persister_notgrouped = as.numeric(combined_table$qval.Persister_notgrouped)
    combined_table$qval.Persister_grouped = as.numeric(combined_table$qval.Persister_grouped)
    WriteXLS(combined_table,file.path(ChromSCape_directory, "Diff_Analysis_Gene_Sets",
                                      paste0("DA_pers_vs_DMSO_sc_bulk.xlsx")))
    
    
    pdf(file.path(plotdir_sup,paste0("Comparison_DA_Pers_vs_DMSO.pdf")),height=5,width=5)
    smoothScatter(y = combined_table$cdiff.Persister_notgrouped,
                        x = combined_table$log2FC.Persister_grouped,colramp=viridis,pch=19,
                  cex=0.25)
    combined_table. = combined_table
    combined_table.$qval.Persister_notgrouped = -log10(combined_table.$qval.Persister_notgrouped)
    combined_table.$qval.Persister_grouped = -log10(combined_table.$qval.Persister_grouped)
    smoothScatter(y = combined_table.$qval.Persister_notgrouped,
                  x = combined_table.$qval.Persister_grouped,colramp=viridis,pch=19,
                  cex=0.25)
    smoothScatter(y = combined_table$log2FC.Persister_grouped,
                  x = combined_table$log2FC.bulk_Pers_vs_DMSO,colramp=viridis,pch=19,
                  cex=0.25)

    combined_table. = combined_table
    combined_table.$qval.Persister_grouped = -log10(combined_table.$qval.Persister_grouped)
    combined_table.$qval.bulk_Pers_vs_DMSO = -log10(combined_table.$qval.bulk_Pers_vs_DMSO)
    smoothScatter(y = combined_table.$qval.Persister_grouped,
                  x = combined_table.$qval.bulk_Pers_vs_DMSO,colramp=viridis,pch=19,
                  cex=0.25)
    
    smoothScatter(y = combined_table.$cdiff.Persister_notgrouped,
                  x = combined_table.$log2FC.bulk_Pers_vs_DMSO,colramp=viridis,pch=19,
                  cex=0.25)
    
    
    v = data.frame(
        "Depleted - sc" = (combined_table$cdiff.Persister_notgrouped < -logFC.th  & combined_table$qval.Persister_notgrouped < qval.th),
        "Depleted - sc (grouped)" = (combined_table$log2FC.Persister_grouped < -logFC.th  & combined_table$qval.Persister_grouped < qval.th),
        "Depleted - bulk" = (combined_table$log2FC.bulk_Pers_vs_DMSO < -1.58  & combined_table$qval.bulk_Pers_vs_DMSO < qval.th)        )

    library(VennDiagram)
    VennDiagram_3(combined_table$Gene[v$Depleted...sc],
                  combined_table$Gene[v$Depleted...sc..grouped.],
                  combined_table$Gene[v$Depleted...bulk], 
                  savePNG=F, "", "Venn", c("Sc", "Sc - Grouped", "Bulk"))
    dev.off()
    
    ### Gene Set Analysis
    scExp_Pers_vs_C3 = gene_set_enrichment_analysis_scExp(scExp_Pers_vs_C3,enrichment_qval = 0.1,
                                                          qval.th = qval.th, cdiff.th = logFC.th)
    
    diff_Pers_vs_C3 = scExp_Pers_vs_C3@metadata$diff$res
    diff_Pers_vs_C3 = diff_Pers_vs_C3[,c(1:4,10:14)]
    colnames(diff_Pers_vs_C3) = gsub("C2","Pers_vs_C3",colnames(diff_Pers_vs_C3))
    
    GSA_Pers_vs_C3 = scExp_Pers_vs_C3@metadata$enr
    
    save(diff_Pers_vs_C3,GSA_Pers_vs_C3,
         file = file.path(ChromSCape_directory, "Diff_Analysis_Gene_Sets",
                          paste0("Pers_vs_C3_",qval.th,"_",logFC.th,"_",de_type,".RData"))) #save the data 
    
    diff_Pers_vs_C3
    diff_Pers_vs_C3$gene = annotK27$gene[match(diff_Pers_vs_C3$ID,gsub(":|-","_",annotK27$ID))]
    
    under_res = diff_Pers_vs_C3 %>% dplyr::filter(cdiff.Pers_vs_C3 < -logFC.th & qval.Pers_vs_C3 < qval.th) %>% 
        dplyr::arrange(qval.Pers_vs_C3) %>% dplyr::select(chr,start,end,,cdiff.Pers_vs_C3,qval.Pers_vs_C3,gene)
    over_res = diff_Pers_vs_C3 %>% dplyr::filter(cdiff.Pers_vs_C3 > logFC.th & qval.Pers_vs_C3 < qval.th) %>% 
        dplyr::arrange(qval.Pers_vs_C3) %>% dplyr::select(chr,start,end,cdiff.Pers_vs_C3,qval.Pers_vs_C3,gene)
    
    WriteXLS(c("over_res","under_res","diff_Pers_vs_C3"),
             ExcelFileName = file.path(
                 ChromSCape_directory,"Diff_Analysis_Gene_Sets",paste0("DA_Pers_vs_C3_logFC_",
                                                                       round(logFC.th,2),".xlsx")),
             SheetNames = c(
                 paste0("Over_",round(logFC.th,2),"_",qval.th,"_n",
                        nrow(over_res)),paste0(
                            "Under_-",round(logFC.th,2),"_",qval.th,
                            "_n",nrow(under_res)),"All"),
             perl = "perl", verbose = FALSE, row.names = FALSE,
             col.names = TRUE, AdjWidth = T, AutoFilter = TRUE,
             BoldHeaderRow = TRUE, na = "", FreezeRow = 1, FreezeCol = 1)
    
    
    # Pathways enriched in depleted loci in Persister
    ChromSCape::table_enriched_genes_scExp(scExp_cf,set = "Underexpressed", cell_cluster = "C5",
                                           enr_class_sel = c("c2_curated","hallmark"))
}

#################################################################################
################ DA & GSA cluster scPersister vs DMSO (C3 + C4) with BULK #######
#################################################################################

scExp_Pers_vs_DMSO = scExp_cf
scExp_Pers_vs_DMSO = scExp_Pers_vs_DMSO[,which(scExp_Pers_vs_DMSO$cell_cluster %in% c("C3","C2","C4"))]
scExp_Pers_vs_DMSO$cell_cluster[which(scExp_Pers_vs_DMSO$cell_cluster %in% c("C3","C4"))] = "C1"

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
resultZ <- read.csv(unzip(file.path(maindir,"input","bulk_ChIPseq","MM468","results_Zerone_merge_hg38_K27_annot.zip"),
                          exdir = tempdir()))
rownames(resultZ) <- paste(resultZ$Chromosome,resultZ$Begin,resultZ$End,sep="_")
resultZ = resultZ[,-c(1,2,3)]
bulk_samp_dmso = c("MM468_2_p15_J113_DMSO_K27", "MM468_2_p16_J113_DMSO_K27",
                   "MM468_3_J77_DMSO_K27", "MM468_5_J67_DMSO_K27", "MM468_dosetest_J3_DMSO_K27")
bulk_samp_persister = c("MM468_2_p15b_J113_5FUR_K27","MM468_2_p15_J113_5FUR_K27",
                        "MM468_3_J77_5FUR_K27", "MM468_5_J67_5FUR_K27")
bulk_RZ = resultZ[,c(bulk_samp_dmso,bulk_samp_persister)]
colnames(bulk_RZ) [which(colnames(bulk_RZ) %in% bulk_samp_dmso)] =
    paste0(colnames(bulk_RZ) [which(colnames(bulk_RZ) %in% bulk_samp_dmso)],"_C1")
colnames(bulk_RZ) [which(colnames(bulk_RZ) %in% bulk_samp_persister)] =
    paste0(colnames(bulk_RZ) [which(colnames(bulk_RZ) %in% bulk_samp_persister)],"_C2")

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
                                   paste0("DA_grouped_Persister_vs_DSMO_with_bulk.xlsx")),
         SheetNames = c(paste0("Over_",round(logFC.th,2),"_",qval.th,"_n",nrow(over_res)),
                        paste0("Under_-",round(logFC.th,2),"_",qval.th,"_n",nrow(under_res)),"All"), perl = "perl", verbose = FALSE, row.names = FALSE, col.names = TRUE, AdjWidth = T, AutoFilter = TRUE, BoldHeaderRow = TRUE, na = "", FreezeRow = 1, FreezeCol = 1)

# Volcano plot
res = readxl::read_xlsx(file.path(ChromSCape_directory, "Diff_Analysis_Gene_Sets",paste0("DA_grouped_Persister_vs_DSMO_with_bulk.xlsx")),sheet = 3)
res = as.data.frame(res)
# colnames(res) = c("seqnames", "start", "end", "id", "cdiff.C2", "pval.C2", "qval.C2","Gene")
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

pdf(file.path(plotdir_sup,"scRNA_scChIP_log2FC_Persister_sc_grouped_vs_DMSO.pdf"))
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
                                         "_logFC",round(logFC.th,2),".xlsx")), 
        SheetNames = c( paste0("Overexp_in_", gp), paste0("Underexp_in_", gp) ),
        perl = "perl", verbose = FALSE, row.names = FALSE, col.names = TRUE,
        AdjWidth = TRUE, AutoFilter = TRUE, BoldHeaderRow = TRUE, na = "",
        FreezeRow = 1, FreezeCol = 1)
}

run_DA_sc = FALSE
if(run_DA_sc) {
    de_type  = "one_vs_rest"
    scExp_Pers_vs_DMSO = differential_analysis_scExp(scExp_Pers_vs_DMSO,de_type = de_type)
    
    ### Gene Set Analysis
    scExp_Pers_vs_DMSO = gene_set_enrichment_analysis_scExp(
        scExp_Pers_vs_DMSO,enrichment_qval = 0.1, qval.th = qval.th, cdiff.th = logFC.th)
    
    diff_Pers_vs_DMSO = scExp_Pers_vs_DMSO@metadata$diff$res
    diff_Pers_vs_DMSO = diff_Pers_vs_DMSO[,c(1:4,10:14)]
    colnames(diff_Pers_vs_DMSO) = gsub("C2","Pers_vs_DMSO",colnames(diff_Pers_vs_DMSO))
    
    GSA_Pers_vs_DMSO = scExp_Pers_vs_DMSO@metadata$enr
    
    save(diff_Pers_vs_DMSO,GSA_Pers_vs_DMSO,
         file = file.path(ChromSCape_directory, "Diff_Analysis_Gene_Sets",
                          paste0("Pers_vs_DMSO_",qval.th,"_",logFC.th,"_",de_type,".RData"))) #save the data 
    
    diff_Pers_vs_DMSO$gene = annotK27$gene[match(diff_Pers_vs_DMSO$ID,gsub(":|-","_",annotK27$ID))]
    
    under_res = diff_Pers_vs_DMSO %>% dplyr::filter(cdiff.Pers_vs_DMSO < -logFC.th & qval.Pers_vs_DMSO < qval.th) %>% 
        dplyr::arrange(qval.Pers_vs_DMSO) %>% dplyr::select(chr,start,end,,cdiff.Pers_vs_DMSO,qval.Pers_vs_DMSO,gene)
    over_res = diff_Pers_vs_DMSO %>% dplyr::filter(cdiff.Pers_vs_DMSO > logFC.th & qval.Pers_vs_DMSO < qval.th) %>% 
        dplyr::arrange(qval.Pers_vs_DMSO) %>% dplyr::select(chr,start,end,cdiff.Pers_vs_DMSO,qval.Pers_vs_DMSO,gene)
    
    WriteXLS(c("over_res","under_res","diff_Pers_vs_DMSO"),
             ExcelFileName = file.path(
                 ChromSCape_directory,"Diff_Analysis_Gene_Sets",paste0("DA_Pers_vs_DMSO_sc_logFC_",
                                                                       round(logFC.th,2),".xlsx")),
             SheetNames = c(
                 paste0("Over_",round(logFC.th,2),"_",qval.th,"_n",
                        nrow(over_res)),paste0(
                            "Under_-",round(logFC.th,2),"_",qval.th,
                            "_n",nrow(under_res)),"All"),
             perl = "perl", verbose = FALSE, row.names = FALSE,
             col.names = TRUE, AdjWidth = T, AutoFilter = TRUE,
             BoldHeaderRow = TRUE, na = "", FreezeRow = 1, FreezeCol = 1)
    
    
    # Pathways enriched in depleted loci in Persister
    ChromSCape::table_enriched_genes_scExp(scExp_Pers_vs_DMSO,set = "Underexpressed", cell_cluster = "C2",
                                           enr_class_sel = c("c2_curated","hallmark"))
    
    # Volcano plot
    res = readxl::read_xlsx(file.path(ChromSCape_directory,"Diff_Analysis_Gene_Sets",paste0("DA_Pers_vs_DMSO_sc_logFC_",
                                                              round(logFC.th,2),".xlsx")),sheet = 3)
    res = as.data.frame(res)
    # colnames(res) = c("seqnames", "start", "end", "id", "cdiff.C2", "pval.C2", "qval.C2","Gene")
    res$qval.Pers_vs_DMSO = as.numeric(res$qval.Pers_vs_DMSO)
    res$pval.Pers_vs_DMSO = as.numeric(res$pval.Pers_vs_DMSO)
    res$cdiff.Pers_vs_DMSO = as.numeric(res$cdiff.Pers_vs_DMSO)
    
    # Anticorrelation with scRNA
    scRNA = readxl::read_xlsx(file.path(maindir, "output","scRNAseq","MM468","Persister",
                                        "Supervised","Tables","Differential_analysis_Limma_logFC_1.58.xlsx"),sheet= 3)
    res_unique_gene = res %>% tidyr::separate_rows(gene) %>% unique %>% dplyr::select(ID, cdiff.Pers_vs_DMSO,
                                                                                      qval.Pers_vs_DMSO,gene)
    
    scRNA_scChIP = left_join(scRNA, res_unique_gene,by =c("Symbol"="gene"))
    top_demethylated <- head(scRNA_scChIP$Symbol[order(scRNA_scChIP$cdiff.Pers_vs_DMSO,decreasing=F)],n=15)
    
    pdf(file.path(plotdir_sup,"scRNA_scChIP_log2FC_Persister_sc_vs_DMSO.pdf"))
    sp <- ggscatter(scRNA_scChIP, x = "cdiff.Pers_vs_DMSO", y = "log2FC.C2_pers",
                    add = "reg.line", ylab = "scRNA log2FC", xlab="scChIP K27 log2FC",  # Add regressin line
                    add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
                    conf.int = TRUE,
                    label="Symbol",repel=TRUE,
                    label.select= top_demethylated) #add size of dot function of initial expression in DMSO
    # Add correlation coefficient
    print(sp + stat_cor(method = "pearson",label.x=1,label.sep = "\n"))
    
    # With only diff chip seq
    scRNA_scChIP_filt = scRNA_scChIP %>% dplyr::filter(abs(cdiff.Pers_vs_DMSO) > logFC.th,
                                                       qval.Pers_vs_DMSO < qval.th)
    sp <- ggscatter(scRNA_scChIP_filt, x = "cdiff.Pers_vs_DMSO", y = "log2FC.C2_pers",
                    add = "reg.line", ylab = "scRNA log2FC", xlab="scChIP K27 log2FC",  # Add regressin line
                    add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
                    conf.int = TRUE,
                    label="Symbol",repel=TRUE,
                    label.select= top_demethylated) #add size of dot function of initial expression in DMSO
    # Add correlation coefficient
    print(sp + stat_cor(method = "pearson",label.x=1,label.sep = "\n"))
    dev.off()
    
    Overexpressed = GSA_Pers_vs_DMSO$Overexpressed[[2]]
    Underexpressed = GSA_Pers_vs_DMSO$Underexpressed[[2]]
    WriteXLS(
        c("Overexpressed", "Underexpressed"),
        ExcelFileName = file.path(ChromSCape_directory, "Diff_Analysis_Gene_Sets",
                                  paste0("Enrichment_test_sc_",gp,"_vs_",refs,
                                         "_logFC",round(logFC.th,2),".xlsx")), 
        SheetNames = c( paste0("Overexp_in_", gp), paste0("Underexp_in_", gp) ),
        perl = "perl", verbose = FALSE, row.names = FALSE, col.names = TRUE,
        AdjWidth = TRUE, AutoFilter = TRUE, BoldHeaderRow = TRUE, na = "",
        FreezeRow = 1, FreezeCol = 1)
}

# 10 k 
run10k = TRUE
if(run10k){
    # 10k mat 
    out <- import_scExp(unzip(file.path(datadir,"Count_Matrices", "MM468_K27_transcripts_10k_1000.zip"),
                              exdir = tempdir()))
    
    # Save raw
    datamatrix_10k = out$datamatrix
    annot_raw_10k = out$annot_raw
    
    mat = NULL
    for(cluster in c("C1","C2")){
        samps = unique(scExp_Pers_vs_DMSO$sample_id[which(scExp_Pers_vs_DMSO$cell_cluster == cluster)])
        tot = table(scExp_Pers_vs_DMSO$sample_id[which(scExp_Pers_vs_DMSO$cell_cluster == cluster)])
        for(sample in names(tot)){
            if(tot[sample]>50){
                cells = scExp_Pers_vs_DMSO$cell_id[which(scExp_Pers_vs_DMSO$cell_cluster == cluster &
                                                             scExp_Pers_vs_DMSO$sample_id == sample)]
                cells = intersect(cells,colnames(datamatrix_10k))
                print(paste0(sample, " - ", length(cells)))
                col = as.matrix(rowSums(datamatrix_10k[,match(cells, colnames(datamatrix_10k))]))
                colnames(col) = paste0(sample,"_",cluster)
                if(is.null(mat)) mat = col else mat = cbind(mat,col)
            } 
        } 
    }
    
    # Read bulk Count matrices
    resultZ <- read.csv(unzip(file.path(maindir,"input","bulk_ChIPseq","MM468","CountTable_bulk_MM468_H3K27me3_TSS.zip"),
                              exdir = tempdir()))
    row = paste(resultZ$Chromosome,resultZ$Begin,resultZ$End,sep="_")
    dups = duplicated(row)
    resultZ = resultZ[!dups,]
    rownames(resultZ) <- row[!dups]
    feature = as.data.frame(resultZ[,c(1,2,3)])
    resultZ = resultZ[,-c(1,2,3)]
    bulk_samp_dmso = c("MM468_2_p15_J113_DMSO_K27", "MM468_2_p16_J113_DMSO_K27",
                       "MM468_3_J77_DMSO_K27", "MM468_5_J67_DMSO_K27")
    bulk_samp_persister = c("MM468_2_p15b_J113_5FUR_K27","MM468_2_p15_J113_5FUR_K27",
                            "MM468_3_J77_5FUR_K27", "MM468_5_J67_5FUR_K27")
    bulk_RZ = resultZ[,c(bulk_samp_dmso,bulk_samp_persister)]
    colnames(bulk_RZ) [which(colnames(bulk_RZ) %in% bulk_samp_dmso)] =
        paste0(colnames(bulk_RZ) [which(colnames(bulk_RZ) %in% bulk_samp_dmso)],"_C1")
    colnames(bulk_RZ) [which(colnames(bulk_RZ) %in% bulk_samp_persister)] =
        paste0(colnames(bulk_RZ) [which(colnames(bulk_RZ) %in% bulk_samp_persister)],"_C2")
    
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
    
    #selection of peaks with at least a log2 RPKM of 1 in one sample
    RZ_sel = mat
    annot = data.frame(sample_id = colnames(mat), cluster = gsub(".*_","",colnames(mat)))
    res <- geco.ChIPseqCompareLimma(mat=RZ_sel,
                                    metadata =annot,
                                    ref=myrefs,
                                    groups=mygps,
                                    featureTab=feature
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
                                       paste0("DA_grouped_Persister_vs_DSMO_with_bulk_10k.xlsx")),
             SheetNames = c(paste0("Over_",round(logFC.th,2),"_",qval.th,"_n",nrow(over_res)),
                            paste0("Under_-",round(logFC.th,2),"_",qval.th,"_n",nrow(under_res)),"All"), perl = "perl", verbose = FALSE, row.names = FALSE, col.names = TRUE, AdjWidth = T, AutoFilter = TRUE, BoldHeaderRow = TRUE, na = "", FreezeRow = 1, FreezeCol = 1)
    
    # Volcano plot
    res = readxl::read_xlsx(file.path(ChromSCape_directory, "Diff_Analysis_Gene_Sets",
                                      paste0("DA_grouped_Persister_vs_DSMO_with_bulk_10k.xlsx")),sheet = 3)
    res = as.data.frame(res)
    # colnames(res) = c("seqnames", "start", "end", "id", "cdiff.C2", "pval.C2", "qval.C2","Gene")
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
    
    pdf(file.path(plotdir_sup,"scRNA_scChIP_log2FC_Persister_sc_grouped_vs_DMSO_10k.pdf"))
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
                                             "_logFC",round(logFC.th,2),"_10k.xlsx")), 
            SheetNames = c( paste0("Overexp_in_", gp), paste0("Underexp_in_", gp) ),
            perl = "perl", verbose = FALSE, row.names = FALSE, col.names = TRUE,
            AdjWidth = TRUE, AutoFilter = TRUE, BoldHeaderRow = TRUE, na = "",
            FreezeRow = 1, FreezeCol = 1)
    }
    
}
#################################################################################
################ DA & GSA cluster scPersister vs DMSO (C3 + C4) (sc - sc)###
#################################################################################

scExp_Pers_vs_DMSO_sc = scExp_cf
scExp_Pers_vs_DMSO_sc = scExp_Pers_vs_DMSO_sc[,which(scExp_Pers_vs_DMSO_sc$cell_cluster %in% c("C3","C2","C4"))]
scExp_Pers_vs_DMSO_sc$cell_cluster[which(scExp_Pers_vs_DMSO_sc$cell_cluster %in% c("C3","C4"))] = "C1"

run_DA_sc = TRUE
if(run_DA_sc) {
    de_type  = "one_vs_rest"
    scExp_Pers_vs_DMSO_sc = differential_analysis_scExp(scExp_Pers_vs_DMSO_sc,de_type = de_type)
    
    ### Gene Set Analysis
    scExp_Pers_vs_DMSO_sc = gene_set_enrichment_analysis_scExp(
        scExp_Pers_vs_DMSO_sc,enrichment_qval = 0.1, qval.th = qval.th, cdiff.th = logFC.th)
    
    diff_Pers_vs_DMSO = scExp_Pers_vs_DMSO_sc@metadata$diff$res
    diff_Pers_vs_DMSO = diff_Pers_vs_DMSO[,c(1:4,10:14)]
    colnames(diff_Pers_vs_DMSO) = gsub("C2","Pers_vs_DMSO",colnames(diff_Pers_vs_DMSO))
    
    GSA_Pers_vs_DMSO = scExp_Pers_vs_DMSO_sc@metadata$enr
    
    save(diff_Pers_vs_DMSO,GSA_Pers_vs_DMSO,
         file = file.path(ChromSCape_directory, "Diff_Analysis_Gene_Sets",
                          paste0("Pers_vs_DMSO_",qval.th,"_",logFC.th,"_",de_type,".RData"))) #save the data 
    
    diff_Pers_vs_DMSO$gene = annotK27$gene[match(diff_Pers_vs_DMSO$ID,gsub(":|-","_",annotK27$ID))]
    
    under_res = diff_Pers_vs_DMSO %>% dplyr::filter(cdiff.Pers_vs_DMSO < -logFC.th & qval.Pers_vs_DMSO < qval.th) %>% 
        dplyr::arrange(qval.Pers_vs_DMSO) %>% dplyr::select(chr,start,end,,cdiff.Pers_vs_DMSO,qval.Pers_vs_DMSO,gene)
    over_res = diff_Pers_vs_DMSO %>% dplyr::filter(cdiff.Pers_vs_DMSO > logFC.th & qval.Pers_vs_DMSO < qval.th) %>% 
        dplyr::arrange(qval.Pers_vs_DMSO) %>% dplyr::select(chr,start,end,cdiff.Pers_vs_DMSO,qval.Pers_vs_DMSO,gene)
    
    WriteXLS(c("over_res","under_res","diff_Pers_vs_DMSO"),
             ExcelFileName = file.path(
                 ChromSCape_directory,"Diff_Analysis_Gene_Sets",paste0("DA_Pers_vs_DMSO_sc_logFC_",
                                                                       round(logFC.th,2),".xlsx")),
             SheetNames = c(
                 paste0("Over_",round(logFC.th,2),"_",qval.th,"_n",
                        nrow(over_res)),paste0(
                            "Under_-",round(logFC.th,2),"_",qval.th,
                            "_n",nrow(under_res)),"All"),
             perl = "perl", verbose = FALSE, row.names = FALSE,
             col.names = TRUE, AdjWidth = T, AutoFilter = TRUE,
             BoldHeaderRow = TRUE, na = "", FreezeRow = 1, FreezeCol = 1)
    
    
    # Pathways enriched in depleted loci in Persister
    ChromSCape::table_enriched_genes_scExp(scExp_Pers_vs_DMSO_sc,set = "Underexpressed", cell_cluster = "C2",
                                           enr_class_sel = c("c2_curated","hallmark"))
    
    # Volcano plot
    res = readxl::read_xlsx(file.path(ChromSCape_directory,"Diff_Analysis_Gene_Sets",paste0("DA_Pers_vs_DMSO_sc_logFC_",
                                                                                            round(logFC.th,2),".xlsx")),sheet = 3)
    res = as.data.frame(res)
    # colnames(res) = c("seqnames", "start", "end", "id", "cdiff.C2", "pval.C2", "qval.C2","Gene")
    res$qval.Pers_vs_DMSO = as.numeric(res$qval.Pers_vs_DMSO)
    res$pval.Pers_vs_DMSO = as.numeric(res$pval.Pers_vs_DMSO)
    res$cdiff.Pers_vs_DMSO = as.numeric(res$cdiff.Pers_vs_DMSO)
    
    # Anticorrelation with scRNA
    scRNA = readxl::read_xlsx(file.path(maindir, "output","scRNAseq","MM468","Persister",
                                        "Supervised","Tables","Differential_analysis_Limma_logFC_1.58.xlsx"),sheet= 3)
    res_unique_gene = res %>% tidyr::separate_rows(gene) %>% unique %>% dplyr::select(ID, cdiff.Pers_vs_DMSO,
                                                                                      qval.Pers_vs_DMSO,gene)
    
    scRNA_scChIP = left_join(scRNA, res_unique_gene,by =c("Symbol"="gene"))
    top_demethylated <- head(scRNA_scChIP$Symbol[order(scRNA_scChIP$cdiff.Pers_vs_DMSO,decreasing=F)],n=15)
    
    pdf(file.path(plotdir_sup,"scRNA_scChIP_log2FC_Persister_sc_vs_DMSO.pdf"))
    sp <- ggscatter(scRNA_scChIP, x = "cdiff.Pers_vs_DMSO", y = "log2FC.C2_pers",
                    add = "reg.line", ylab = "scRNA log2FC", xlab="scChIP K27 log2FC",  # Add regressin line
                    add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
                    conf.int = TRUE,
                    label="Symbol",repel=TRUE,
                    label.select= top_demethylated) #add size of dot function of initial expression in DMSO
    # Add correlation coefficient
    print(sp + stat_cor(method = "pearson",label.x=1,label.sep = "\n"))
    
    # With only diff chip seq
    scRNA_scChIP_filt = scRNA_scChIP %>% dplyr::filter(abs(cdiff.Pers_vs_DMSO) > logFC.th,
                                                       qval.Pers_vs_DMSO < qval.th)
    sp <- ggscatter(scRNA_scChIP_filt, x = "cdiff.Pers_vs_DMSO", y = "log2FC.C2_pers",
                    add = "reg.line", ylab = "scRNA log2FC", xlab="scChIP K27 log2FC",  # Add regressin line
                    add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
                    conf.int = TRUE,
                    label="Symbol",repel=TRUE,
                    label.select= top_demethylated) #add size of dot function of initial expression in DMSO
    # Add correlation coefficient
    print(sp + stat_cor(method = "pearson",label.x=1,label.sep = "\n"))
    dev.off()
    
    Overexpressed = GSA_Pers_vs_DMSO$Overexpressed[[2]]
    Underexpressed = GSA_Pers_vs_DMSO$Underexpressed[[2]]
    WriteXLS(
        c("Overexpressed", "Underexpressed"),
        ExcelFileName = file.path(ChromSCape_directory, "Diff_Analysis_Gene_Sets",
                                  paste0("Enrichment_test_sc_",gp,"_vs_",refs,
                                         "_logFC",round(logFC.th,2),".xlsx")), 
        SheetNames = c( paste0("Overexp_in_", gp), paste0("Underexp_in_", gp) ),
        perl = "perl", verbose = FALSE, row.names = FALSE, col.names = TRUE,
        AdjWidth = TRUE, AutoFilter = TRUE, BoldHeaderRow = TRUE, na = "",
        FreezeRow = 1, FreezeCol = 1)
}

########################## ###################################
# Specific enrichment persister vs dmso
metadata$sample_id
scExp_subset = scExp_cf[,scExp_cf$sample_id %in% c("MM468_5FU6_day33","MM468_DMSO3_day77")]
cells_sampled = subsample_scExp(scExp_subset,400)$cell_id
scExp_subset = scExp_subset[,scExp_subset$cell_id %in% cells_sampled]
scExp_subset$cell_cluster = c(rep("C1",400),rep("C2",400))
scExp_subset = differential_analysis_scExp(scExp_subset,de_type = de_type)
scExp_subset = gene_set_enrichment_analysis_scExp(scExp_subset, enrichment_qval = 0.1, qval.th = qval.th, logFC.th = logFC.th)

ChromSCape::table_enriched_genes_scExp(scExp_subset,set = "Underexpressed", cell_cluster = "C1",
                                       enr_class_sel = c("c2_curated","hallmark"))

diff = scExp_subset@metadata$diff
enr = scExp_subset@metadata$enr
tab = table(scExp_subset$sample_id,scExp_subset$cell_cluster)
save(diff, enr, tab, file = file.path(plotdir_sup, paste0("Persister_vs_Initial_log2FC",logFC.th,".RData"))) #save the data 

## Enrichment
scExp_cf = gene_set_enrichment_analysis_scExp(scExp_cf,enrichment_qval = 0.1, qval.th = qval.th, logFC.th = logFC.th)
save(scExp_cf,file = file.path(ChromSCape_directory, "Diff_Analysis_Gene_Sets",
                               paste0(prefix, "_",nclust,"_",qval.th,"_",logFC.th,"_",de_type,".RData"))) #save the data 


# Pathways enriched in depleted loci in Persister
ChromSCape::table_enriched_genes_scExp(scExp_cf,set = "Underexpressed", cell_cluster = "C5",
                                       enr_class_sel = c("c2_curated","hallmark"))

# Specific enrichment persister vs dmso
metadata$sample_id
scExp_subset = scExp_cf[,scExp_cf$sample_id %in% c("MM468_5FU6_day33","MM468_DMSO3_day77")]
cells_sampled = subsample_scExp(scExp_subset,400)$cell_id
scExp_subset = scExp_subset[,scExp_subset$cell_id %in% cells_sampled]
scExp_subset$cell_cluster = c(rep("C1",400),rep("C2",400))
scExp_subset = differential_analysis_scExp(scExp_subset,de_type = de_type)
scExp_subset = gene_set_enrichment_analysis_scExp(scExp_subset, enrichment_qval = 0.1, qval.th = qval.th, logFC.th = logFC.th)

ChromSCape::table_enriched_genes_scExp(scExp_subset,set = "Underexpressed", cell_cluster = "C1",
                                       enr_class_sel = c("c2_curated","hallmark"))

diff = scExp_subset@metadata$diff
enr = scExp_subset@metadata$enr
tab = table(scExp_subset$sample_id,scExp_subset$cell_cluster)
save(diff, enr, tab, file = file.path(plotdir_sup, paste0("Persister_vs_Initial_log2FC",logFC.th,".RData"))) #save the data 


################################################################################
################ Enrichment in TSS #############################################
################################################################################
# Load bulk + sc grouped peak DA 
res = readxl::read_xlsx(file.path(ChromSCape_directory, "Diff_Analysis_Gene_Sets",
                                  "DA_grouped_Persister_vs_DSMO_with_bulk.xlsx"), sheet = 3)
# Load annotation peak affectation 
annotK27_peak_affectation = read.table( file.path(maindir, "annotation", "annotK27_peak_affectation.tsv"), 
                                        sep = "\t", header = TRUE)
res$peak_affectation = annotK27_peak_affectation$peak_affectation[match(res$id,annotK27_peak_affectation$ID)]
col = data.frame("col"=c("#999999ff","#dededeff","#3b9ab2ff","#78b7c5ff","#ebcc2aff","#e1af00ff"))
sig_diff <- (abs(res$log2FC.Persister)> logFC.th & res$qval.Persister < qval.th)
sig_under <- (res$log2FC.Persister< -logFC.th & res$qval.Persister < qval.th)
sig_over <- (res$log2FC.Persister> logFC.th & res$qval.Persister < qval.th)
l = list("Diff" = sig_diff, "Under" = sig_under, "Over" = sig_over)

for(i in c("Diff","Under","Over")){
    sig = l[[i]]
    png(file.path(plotdir,paste0("enrichment_HQ_",i,".png")),width = 1500,height = 1500,res=300)
    sum_sig = sum(table(res$peak_affectation[sig]))
    sum_all = sum(table(res$peak_affectation))
    cbind(as.data.frame(log2((table(res$peak_affectation[sig])/sum_sig) / (table(res$peak_affectation)/sum_all))),col) %>% ggplot() + geom_bar(aes(x=Var1, y=Freq,fill= Var1) , width = 0.5,stat="identity") + ggtitle(paste0(i," enriched H3K27me3 peaks")) +
        theme_classic() + ggtitle(paste0(i," enriched H3K27me3 peaks")) + scale_fill_manual(values =col$col)  +
        theme(axis.text.x = element_text(angle = 90, vjust =1), legend.position = "None") + ylab("Log2 Enrichment")
    dev.off()
}

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
    mutate(is_significative = ifelse(p <0.01, TRUE,FALSE),
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
    mutate(is_significative = ifelse(p <0.01, TRUE,FALSE),
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