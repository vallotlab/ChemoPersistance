library(here)
source(file.path(here(),"Scripts","functions.R"))
source(file.path(here(),"Scripts","global_var.R"))

maindir = here()

######VARIABLES
options(width=180)
resdir <- file.path(maindir,"output","scRNAseq","Normal_Mammary_Gland")
PDX <- file.path(maindir,"output","scRNAseq","PDX","Supervised","RData")
MM468 <- file.path(maindir,"output","scRNAseq","MM468","Persister","Supervised","RData")
MM468_PDX <- file.path(maindir,"output","scRNAseq","MM468_PDX")
RDatadir_unsup <- file.path(maindir,"output","scRNAseq","PDX","Unsupervised","RData")
RDataSupdir <-  file.path(resdir,"RData");if(!file.exists(RDataSupdir)){dir.create(RDataSupdir)}


## Import annotation file
load(file.path(MM468_PDX,"Overexpressed_GeneLists_MM468_PDX_1.58.RData"))
load(file.path(MM468_PDX,"common_over_genes_pers_vs_unt_log2FC1.58.RData"))
common_overexpressed_genes

PDX = new.env()
MM468 = new.env()
load(file.path(RDatadir_PDX,"Supervised_res_object_edgeR.Rdata"),PDX)
load(file.path(RDatadir_MM468,"Supervised_res_object_edgeR.Rdata"),MM468)

# files = list.files(file.path(resdir,"Data"))
# lists = NULL
# lists = sapply(files, function(x) 
#   if(is.null(lists)) lists = readSparseCounts(file.path(resdir,"Data",x)) else lists = Matrix::cbind2(lists,readSparseCounts(file.path(resdir,"Data",x)) ))
# 
# tmp = lists[[4]]
# for(i in 5:length(lists)){
#   tmp = Matrix::cbind2(tmp, lists[[i]])
# }
# sapply(lists, ncol)
# RawCounts_NormalMammaryGland = tmp
# 
# rownames(RawCounts_NormalMammaryGland) = gsub('\"',"",rownames(RawCounts_NormalMammaryGland))
# 
# # Select Individual 4 only as strong batch effect is present
# # RawCounts_NormalMammaryGland = RawCounts_NormalMammaryGland[,which(annot$sample_id=="Individual_4")]
# 
# gene_metadata = data.frame("gene_short_name" = rownames(RawCounts_NormalMammaryGland),
#                               "total_counts" = rowSums(RawCounts_NormalMammaryGland)
#                            )
# # annot = data.frame("cell_id" = colnames(RawCounts_NormalMammaryGland),
# #                    "sample_id" = c(rep("Individual_4",4116),rep("Individual_5",6964),rep("Individual_6",6014),rep("Individual_7",7552)),
# #                    "total_counts" = colSums(RawCounts_NormalMammaryGland),
# #                    "VIM" = as.numeric(RawCounts_NormalMammaryGland[which(rownames(RawCounts_NormalMammaryGland) %in% "VIM"),]),
# #                    "KRT17" = as.numeric(RawCounts_NormalMammaryGland[which(rownames(RawCounts_NormalMammaryGland) %in% "KRT17"),]),
# #                    "KRT14" = as.numeric(RawCounts_NormalMammaryGland[which(rownames(RawCounts_NormalMammaryGland) %in% "KRT14"),]),
# #                    "NNMT" = as.numeric(RawCounts_NormalMammaryGland[which(rownames(RawCounts_NormalMammaryGland) %in% "NNMT"),]))
# annot = data.frame("cell_id" = colnames(RawCounts_NormalMammaryGland),
#                    "sample_id" = rep("Individual_4",4116),
#                    "total_counts" = colSums(RawCounts_NormalMammaryGland),
#                    "VIM" = as.numeric(RawCounts_NormalMammaryGland[which(rownames(RawCounts_NormalMammaryGland) %in% "VIM"),]),
#                    "KRT17" = as.numeric(RawCounts_NormalMammaryGland[which(rownames(RawCounts_NormalMammaryGland) %in% "KRT17"),]),
#                    "KRT14" = as.numeric(RawCounts_NormalMammaryGland[which(rownames(RawCounts_NormalMammaryGland) %in% "KRT14"),]),
#                    "NNMT" = as.numeric(RawCounts_NormalMammaryGland[which(rownames(RawCounts_NormalMammaryGland) %in% "NNMT"),]))
# 
# save(RawCounts_NormalMammaryGland, gene_metadata,annot,file = file.path(resdir,"RawCounts_NormalMammaryGland_Ind4.RData"))
load(file.path(resdir,"RData","RData_individual_4","RawCounts_NormalMammaryGland_Ind4.RData"))
cds <- new_cell_data_set(RawCounts_NormalMammaryGland,
                         cell_metadata = annot,
                         gene_metadata  = gene_metadata)
gc()

cds <- preprocess_cds(cds,method='PCA', norm_method='size_only',num_dim=50)
cds <- reduce_dimension(cds, reduction_method = 'UMAP')
cds <-  cluster_cells(cds)
  
plot_cells(cds, color_cells_by = "total_counts")
plot_cells(cds,genes=c("VIM"),group_cells_by="partition",show_trajectory_graph=F)
plot_cells(cds,genes=c("KRT17"),group_cells_by="partition",show_trajectory_graph=F)
plot_cells(cds,genes=c("KRT14"),group_cells_by="partition",show_trajectory_graph=F)
plot_cells(cds,genes=c("NNMT"),group_cells_by="partition",show_trajectory_graph=F)

umap_res <- reducedDims(cds)[[2]]
pca_res <- reducedDims(cds)[[1]]

# save(umap_res,file = file.path(resdir,"umap.RData"))
# save(pca_res,file = file.path(resdir,"pca.RData"))

load(file.path(resdir,"RData","RData_individual_4","umap.RData"))
load(file.path(resdir,"RData","RData_individual_4","pca.RData"))

annot$louvain_partition <- cds@clusters clusters[[1]]$clusters
annot$louvain_partition <- paste0("C",annot$louvain_partition)

NormCounts <- t(t(exprs(cds)) /  pData(cds)[, 'Size_Factor'])
LogCounts <- log(NormCounts+1,2)
rm(NormCounts)
gc()

# save(LogCounts,file=file.path(resdir,"LogCounts.RData"))
load(file=file.path(resdir,"RData","RData_individual_4","LogCounts.RData"))

RawCounts = exprs(cds)

# save(RawCounts,file=file.path(resdir,"RData_individual_4","RawCounts.RData"))
load(file=file.path(resdir,"RData","RData_individual_4","RawCounts.RData"))

for(gene in c("KRT18","KRT14","SLPI","ANKRD30A","CD74","LTF","AGR2",
              "SAA2","KRT23","VIM","ESAM","KRT17","ACTA2","TGFBR2","INHBB",common_overexpressed_genes)){
  annot[,gene] <- LogCounts[gene,annot$cell_id]
}

annot$sample_id = as.factor(annot$sample_id)
annot$cell_type = ""
annot$cell_type[which(annot$louvain_partition == "C2")] = "Basal"
annot$cell_type[which(annot$louvain_partition == "C1")] = "Luminal 1"
annot$cell_type[which(annot$louvain_partition == "C3")] = "Luminal 2"
annot$cell_type[which(annot$louvain_partition == "C4")] = "Other"

# save(annot,gene_metadata,file=file.path(resdir,"gene_cell_annot.RData"))
load(file=file.path(resdir,"RData","RData_individual_4","gene_cell_annot.RData"))

anocol <- geco.annotToCol4(annot[,-1],annotT=annot[,-1],plotLegend=F,
                           categCol=NULL, scale_q = "inferno")

umap_res <- umap_res[annot$cell_id,]
save(anocol,file= file.path(resdir,'anocol.RData'))
load(file= file.path(resdir,"RData_individual_4",'anocol.RData'))

for(j in 1:ncol(anocol))
{
  png(file.path(resdir,"UMAPs",paste0("UMAP_",colnames(anocol)[j],".png")), height=1350,width=1200,res=300) 
  if(class(annot[1,colnames(anocol)[j]])=="numeric"){
    plot((umap_res[annot$cell_id,]), col=alpha(anocol[,j],0.2),pch=20,
         cex=0.4,main=paste0(colnames(anocol)[j],
                             " min=",round(min(annot[,colnames(anocol)[j]]),digits=3),
                             " max=",round(max(annot[,colnames(anocol)[j]]),digits=3)),
         xlab="component 1",ylab="component 2")
    dev.off()
  } else{
    plot((umap_res[annot$cell_id,]), col=alpha(anocol[,j],0.2),pch=20,
         cex=0.4,main=paste0(colnames(anocol)[j]),
         xlab="component 1",ylab="component 2")
    dev.off()
  }
}


