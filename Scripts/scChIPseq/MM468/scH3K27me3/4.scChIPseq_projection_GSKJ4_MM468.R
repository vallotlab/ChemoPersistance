library(here)
source(file.path(here(),"Scripts","functions.R"))
source(file.path(here(),"Scripts","global_var.R"))

# Directories -------------------------------------------------------------
maindir = here()
resdir <- file.path(maindir,"output","scChIPseq","MM468","scH3K27me3"); if(!file.exists(resdir)){dir.create(resdir)}
unsupervised_dir  <- file.path(resdir,"Unsupervised","RData")
input_dir  <- file.path(maindir, "input","scChIPseq","MM468")

# Unsupervised directories
resdir <- file.path(resdir, "Supervised") ; if(!file.exists(resdir)){dir.create(resdir)}
resdir_GSKJ4 = file.path(resdir, "GSKJ4") ; if(!file.exists(resdir_GSKJ4)){dir.create(resdir_GSKJ4)}
dataset_name = "MM468_H3K27me3_peaks"
ref_genome ="hg38"

metadata = data.frame(sample_id = c(
    "MM468_DMSO1_day60", "MM468_DMSO3_day77", "MM468_DMSO5_day131",
    "MM468_5FU1_day33", "MM468_5FU2_day67",
    "MM468_5FU6_day131", "MM468_5FU3_day147", "MM468_5FU2_day171"),
    sample_id_color = c("#dfdfdfff", "#999999ff","#363636",
                        "#118675ff", "#8cc453ff",
                        "#ff5722ff", "#feb40fff", "#fd8508ff"))

##########################
## Projection of GSKJ4 onto the UMAP  
##########################

# Load correlation filtered scExp
scExp_cf = qs::qread(file = file.path(unsupervised_dir, "MM468_H3K27me3_peaks_3000_10000_95_uncorrected_correlation_filtered.qs"))


# Read in GSKJ4 peak matrix
out <- import_scExp(list.files(file.path(input_dir,"Count_Matrices"), pattern = "GSKJ4", full.names = T),
                       remove_pattern = "_H3K27me3_peaks")

# Save raw
datamatrix = out$datamatrix
annot_raw = out$annot_raw

annot_raw$cell_id = gsub("MM468_","MM468_GSKJ4_day91",annot_raw$cell_id)
annot_raw$sample_id = gsub("MM468_","MM468_GSKJ4_day91",annot_raw$sample_id)
colnames(datamatrix) = annot_raw$cell_id

scExp <- create_scExp(datamatrix, annot_raw, remove_zero_cells = T, remove_zero_features = FALSE)

scExp = scExp[match(rownames(scExp_cf),rownames(scExp)),]
scExp = normalize_scExp(scExp, type = "CPM")
scExp = reduce_dims_scExp(scExp, dimension_reductions = c("PCA","UMAP"), n = 50)

# Get Rotations
x = t(normcounts(scExp_cf))
x.means <- Matrix::colMeans(x)
svd.0 <- irlba::irlba(x, center = x.means, nv = 50, work = 3)
pca <- svd.0$u %*% diag(svd.0$d)

# Center with all rows:
combined_mat = cbind(normcounts(scExp_cf), normcounts(scExp))
combined_mat = t(combined_mat)
combined_mat = scale(combined_mat, center = TRUE, scale = FALSE)
norm_mat_GSKJ4 = combined_mat[grep("GSK", rownames(combined_mat)),]

# Center only GSKJ4:
norm_mat_GSKJ4_centered_alone = scale(t(normcounts(scExp)), scale=F)

# Project into scExp_cf PCA
projection_GSKJ4 = norm_mat_GSKJ4 %*% svd.0$v
projection_all = combined_mat %*% svd.0$v
projection_GSKJ4_centered_alone = norm_mat_GSKJ4_centered_alone %*% svd.0$v

pca_with_GSKJ4 = rbind(pca, projection_GSKJ4)

# project test with known label
norm_mat_5FU_day33 = scale(t(normcounts(scExp_cf)),scale =F)
norm_mat_5FU_day33 = norm_mat_5FU_day33[grep("5FU1_day33",rownames(norm_mat_5FU_day33)),]
projection_5FU_day33 = norm_mat_5FU_day33 %*% svd.0$v

# Re-run UMAP to run UMAP.predict
config <- umap::umap.defaults
config$metric <- "cosine"
umap <- umap::umap(pca_with_GSKJ4, config = config, method = "naive")
p = predict(umap, projection_GSKJ4)

png(file.path(resdir_GSKJ4,"UMAP_sample_with_GSKJ4_projection.png"), res=300,height=2000,width=2000)
plot(umap$layout[,1],umap$layout[,2], cex = 1.5, pch =20,
     col = alpha(c(scExp_cf$sample_id_color, rep("purple",nrow(norm_mat_GSKJ4))),0.4),
     xlab = "UMAP1", ylab = "UMAP2")
dev.off()

cor_with_GSKJ4 = coop::pcor(t(pca_with_GSKJ4))
clusters = paste0("C", cutree(hclust(as.dist(1-cor_with_GSKJ4),method =  "ward.D"),4))
names(clusters) = c(rownames(x), rownames(projection_GSKJ4))
scExp_cf$cell_cluster = clusters[match(scExp_cf$cell_id, names(clusters))]
anocol_cluster = annotToCol2(as.data.frame(clusters))
barplot(rep(1,4), col = names(table(anocol_cluster)), names.arg = names(table(anocol_cluster)))
anocol_cluster[which(anocol_cluster[,1] == "#F4B400"),1] = "#3cb371ff" # resistant
anocol_cluster[which(anocol_cluster[,1] == "#4285F4"),1] = "#ffd700ff" # persister
anocol_cluster[which(anocol_cluster[,1] == "#DB4437"),1] = "#ffb5c5ff" # DMSO C2
anocol_cluster[which(anocol_cluster[,1] == "#0F9D58"),1] =  "#668b8bff" # DMSO C4

repartition = table(clusters, gsub("_BC.*","",colnames(cor_with_GSKJ4)))
table(scExp_cf$cell_cluster, scExp_cf$sample_id)[c(2,4),6:8]
repartition = t(cbind(table(scExp_cf$cell_cluster, scExp_cf$sample_id)[c(2,4),6:8], repartition[c(2,4),2]))
rownames(repartition)[4] = "MM468_GSKJ4_day91"
repartition = repartition / rowSums(repartition)

png(file.path(resdir_GSKJ4,"Ratio_GSKJ4_cell_in_each_cluster.png"), res=300,height=2000,width=2000)
(100*repartition) %>% as.data.frame() %>% tibble::rownames_to_column("sample_id") %>% tidyr::gather("cell_cluster","nCell",-"sample_id") %>%
    ggplot(aes(x = sample_id, y = nCell, fill = cell_cluster)) +
    geom_bar(stat="identity", size = 1) +
    scale_fill_manual(values = c("#ffb5c5ff","#668b8bff" )) +
    theme_classic() + theme(text = element_text(size = 15), axis.text.x = element_text(angle = 90)) + xlab("") + ylab("% Cluster")
dev.off()

repartition = table(clusters, gsub("_BC.*","",colnames(cor_with_GSKJ4)))
repartition = t(cbind(table(scExp_cf$cell_cluster, scExp_cf$sample_id)[c(2,4),6:8], repartition[c(2,4),2]))
rownames(repartition)[4] = "MM468_GSKJ4_day91"

png(file.path(resdir_GSKJ4,"Number_GSKJ4_cell_in_each_cluster.png"), res=300,height=2000,width=2000)
(repartition) %>% as.data.frame() %>% tibble::rownames_to_column("sample_id") %>% tidyr::gather("cell_cluster","nCell",-"sample_id") %>%
    ggplot(aes(x = sample_id, y = nCell, fill = cell_cluster)) +
    geom_bar(stat="identity", size = 1) +
    scale_fill_manual(values = c("#ffb5c5ff","#668b8bff" )) +
    theme_classic() + theme(text = element_text(size = 15), axis.text.x = element_text(angle = 90)) + xlab("") + ylab("N cells")
dev.off()
