#analysis fold-change MM468
library(here)

maindir= here()
resdir <- file.path(maindir,"output","scRNAseq","MM468")
plotdir <- file.path(resdir,"CompareFC");if(!file.exists(plotdir)){dir.create(plotdir)}

source(file.path(maindir,"Scripts","global_var.R"))

# Persister log2FC
persister = new.env()
load(file.path(resdir, "Persister","Supervised","RData","Supervised_res_object_edgeR.Rdata"), persister)
# UNC log2FC
UNC = new.env()
load(file.path(resdir, "UNC","Supervised","RData","Supervised_res_object_edgeR.Rdata"), UNC)

common = intersect(persister$my.res$Symbol,UNC$my.res$Symbol)
persister$my.res = persister$my.res[match(common,persister$my.res$Symbol),]
UNC$my.res = UNC$my.res[match(common,UNC$my.res$Symbol),]

png(file.path(plotdir,"log2FC_persister_vs_UNC.png"), width = 1200, height = 1200, res = 300)
smoothScatter(UNC$my.res$log2FC.UNC, persister$my.res$log2FC.C2_pers,cex=1, 
              colramp=colorRampPalette(mycolramp),xlab="log2FC UNC vs initial",nbin=300,
              ylab="log2FC 5FU6 persister vs initial",bandwidth=c(0.01,0.01),xlim=c(-5,5),ylim=c(-5,5)) 
dev.off()

png(file.path(plotdir,"log2FC_resistant_vs_UNC.png"), width = 1200, height = 1200, res = 300)
smoothScatter(UNC$my.res$log2FC.UNC, persister$my.res$log2FC.MM468_5FU6_day214,cex=1, 
              colramp=colorRampPalette(mycolramp),xlab="log2FC UNC vs initial",nbin=300,
              ylab="log2FC 5FU6 resistant vs initial",bandwidth=c(0.01,0.01),xlim=c(-5,5),ylim=c(-5,5)) 
dev.off()

