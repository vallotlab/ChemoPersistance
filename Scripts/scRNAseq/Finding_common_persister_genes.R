##### PATHS ##### 
library(here)
source(file.path(here(),"Scripts","global_var.R"))


maindir = file.path(here(),"output","scRNAseq")
resdir = file.path(maindir,"MM468_PDX"); if(!dir.exists(resdir)) dir.create(resdir)
RDataDir_MM468 <- file.path(maindir,"MM468","Persister","Supervised","RData")
RDataDir_PDX <- file.path(maindir,"PDX","Supervised","RData")
MM468_env <- new.env()
PDX_env <- new.env()

## Import annotation file
load(file.path(RDataDir_MM468,"Supervised_res_object_edgeR.Rdata"),MM468_env)
load(file.path(RDataDir_PDX,"Supervised_res_object_edgeR.Rdata"),PDX_env)

log2FC_thresholds <- log2(c(2,3,4))
Signif_threshold <- 0.01

for(log2FC_threshold in log2FC_thresholds){
  MM468_env$common_overexpressed_genes <- MM468_env$my.res[which(
    MM468_env$my.res$log2FC.C2_pers > log2FC_threshold & 
      MM468_env$my.res$qval.C2_pers < Signif_threshold ),]
  PDX_env$common_overexpressed_genes <- PDX_env$my.res[which(
    PDX_env$my.res$log2FC.persister_6_vs_UNT > log2FC_threshold & 
      PDX_env$my.res$qval.persister_6_vs_UNT < Signif_threshold ),]
  
  common_overexpressed_genes <- intersect(MM468_env$common_overexpressed_genes$Symbol,
                                        PDX_env$common_overexpressed_genes$Symbol)
  
  save(common_overexpressed_genes, file=file.path(
    resdir,paste0("common_over_genes_pers_vs_unt_log2FC", round(log2FC_threshold,2), ".RData")))
  
  rownames(PDX_env$common_overexpressed_genes) = PDX_env$common_overexpressed_genes$Symbol
  
  # Log2FC
  ctab = cbind(MM468_env$common_overexpressed_genes[common_overexpressed_genes,c("Symbol","log2FC.C2_pers")],
               PDX_env$common_overexpressed_genes[common_overexpressed_genes,"log2FC.persister_6_vs_UNT"])
  colnames(ctab) = c("Symbol","log2FC.MM468","log2FC.PDX")
  top = order(ctab$log2FC.MM468 + ctab$log2FC.PDX,decreasing = T)[1:5]
  
  png(file.path(resdir, paste0("common_over_genes_pers_vs_unt_log2FC", round(log2FC_threshold,2),".png")),
      width = 2000,height = 2000,res =300)
  sp <- ggscatter(ctab, y = "log2FC.MM468", x = "log2FC.PDX",
                  add = "reg.line",  # Add regressin line
                  add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
                  conf.int = TRUE,
                  label="Symbol",repel=TRUE,
                  label.select= ctab$Symbol[top]) #add size of dot function of initial expression in DMSO
  # Add correlation coefficient
  print(sp + stat_cor(method = "pearson",label.x=1,label.sep = "\n"))
  dev.off()
  
  # Q-value
  ctab = cbind(MM468_env$common_overexpressed_genes[common_overexpressed_genes,c("Symbol","qval.C2_pers")],
               PDX_env$common_overexpressed_genes[common_overexpressed_genes,"qval.persister_6_vs_UNT"])
  colnames(ctab) = c("Symbol","-log10_qval.MM468","-log10_qval.PDX")
  
  ctab$`-log10_qval.MM468` = -log10(ctab$`-log10_qval.MM468`)
  ctab$`-log10_qval.PDX` = -log10(ctab$`-log10_qval.PDX`)
  top = order( ctab$`-log10_qval.MM468` +  ctab$`-log10_qval.PDX`, decreasing = T)[1:5]
  png(file.path(resdir, paste0("common_over_genes_pers_vs_unt_Log10qvalue_log2FC", round(log2FC_threshold,2),".png")),
      width = 2000,height = 2000,res =300)
  sp <- ggscatter(ctab, y = "-log10_qval.MM468", x = "-log10_qval.PDX",
                  add = "reg.line",  # Add regressin line
                  add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
                  conf.int = TRUE,
                  label="Symbol",repel=TRUE,
                  label.select= ctab$Symbol[top]) #add size of dot function of initial expression in DMSO
  # Add correlation coefficient
  print(sp + stat_cor(method = "pearson",label.x=1,label.sep = "\n"))
  dev.off()
  
  print(length(MM468_env$common_overexpressed_genes))
  print(length(PDX_env$common_overexpressed_genes))
  print(length(common_overexpressed_genes))
  
}


## Pathways in common
for(log2FC_threshold in log2FC_thresholds){
  load(file.path(RDataDir_MM468, paste0(
    "Enrichment_test_C2_pers_vs_DMSO_MSigDB_logFC",round(log2FC_threshold,2),".RData")),
    MM468_env)
  load(file.path(RDataDir_PDX, paste0(
    "Enrichment_test_persister_6_vs_UNT_vs_UNT_MSigDB_logFC",round(log2FC_threshold,2),".RData")),
    PDX_env)
  MM468_env$Overexpressed = MM468_env$Overexpressed %>% 
    dplyr::filter(Class %in% c("c2_curated","hallmark","c5_GO"))
  PDX_env$Overexpressed = PDX_env$Overexpressed %>% 
    dplyr::filter(Class %in% c("c2_curated","hallmark","c5_GO"))
  common_overexpressed_pathways <- intersect(MM468_env$Overexpressed$Gene.Set,
                                           PDX_env$Overexpressed$Gene.Set)
  save(common_overexpressed_pathways, file=file.path(
    resdir,paste0("common_over_pathways_pers_vs_unt_log2FC", round(log2FC_threshold,2), ".RData")))
  
  rownames(MM468_env$Overexpressed) = MM468_env$Overexpressed$Gene.Set
  rownames(PDX_env$Overexpressed) = PDX_env$Overexpressed$Gene.Set
  # Q-value
  ctab = cbind(MM468_env$Overexpressed[common_overexpressed_pathways,c("Gene.Set","q-value")],
               PDX_env$Overexpressed[common_overexpressed_pathways,"q-value"])
  colnames(ctab) = c("Gene.Set","-log10_qval.MM468","-log10_qval.PDX")
  
  ctab$`-log10_qval.MM468` = -log10(ctab$`-log10_qval.MM468`)
  ctab$`-log10_qval.PDX` = -log10(ctab$`-log10_qval.PDX`)
  top = order( ctab$`-log10_qval.MM468` +  ctab$`-log10_qval.PDX`, decreasing = T)[1:5]
  png(file.path(resdir, paste0("common_over_genes_pers_vs_unt_qvalue_log2FC", round(log2FC_threshold,2),".png")),
      width = 2000,height = 2000,res =300)
  sp <- ggscatter(ctab, y = "-log10_qval.MM468", x = "-log10_qval.PDX",
                  add = "reg.line",  # Add regressin line
                  add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
                  conf.int = TRUE,
                  label="Gene.Set",repel=TRUE,
                  label.select= ctab$Gene.Set[top]) #add size of dot function of initial expression in DMSO
  # Add correlation coefficient
  print(sp + stat_cor(method = "pearson",label.x=1,label.sep = "\n"))
  dev.off()
  
  print(nrow(MM468_env$Overexpressed))
  print(nrow(PDX_env$Overexpressed))
  print(length(common_overexpressed_pathways))
  
}
