library(here)

maindir= here()
source(file.path(maindir, "Scripts", "global_var.R"))
source(file.path(maindir, "Scripts", "functions.R"))
inputdir = file.path(maindir,"input","bulk_ChIPseq","MM468")
outputdir = file.path(maindir,"output","bulk_ChIPseq","MM468","ChromatinIndexing")
if(!dir.exists(outputdir)) dir.create(outputdir)

# Count matrix on MM468 K27 peaks
mat = read.csv(unzip(file.path(inputdir, "results_ChromatinIndexing_poolK27_MM468_K27_peaks.csv.zip"),
            exdir = tempdir()))

colnames(mat)
rownames(mat) = paste0(mat$Chromosome,":",mat$Begin,"-",mat$End)
mat = mat[,-c(1:3)]

# Select DMSO, UNC & UNC+5FU
mat =  mat[,grep("K27",colnames(mat))]

# Chromatin Indexing stats
ratio = read.csv(file.path(outputdir,"ratio_ip_input_K27.csv"))
rownames(ratio) = ratio$Sample

norm_mat = mat
for(i in 1:ncol(norm_mat)){
  norm_mat[,i] = 10^6 * ratio[i,"value"] * (norm_mat[,i] / sum(norm_mat[,i]))
}

mycolramp <- c("white",viridis(n=4))
png(file.path(outputdir,"log2_Count_UNC_vs_DMSO.png"), width = 1500, height = 1500, res = 300)
smoothScatter(x = log2(norm_mat[,"MM468BC_DMSO_D33_K27"]),y = log2(norm_mat[,"MM468BC_UNC_D33_K27"]),cex=1, 
              colramp=colorRampPalette(mycolramp),xlab="DMSO D33 K27",nbin=1000,
              ylab="UNC D33 K27",bandwidth=c(0.01,0.01), xlim =c(0,13), ylim=c(0,13)) 
abline(a = 0, b=1, lty=2, lw = 2)
dev.off()

png(file.path(outputdir,"log2_Count_UNC_5FU_vs_DMSO.png"), width = 1500, height = 1500, res = 300)
smoothScatter(x = log2(norm_mat[,"MM468BC_DMSO_D33_K27"]),y = log2(norm_mat[,"MM468BC_UNC_5FU_D33_K27"]),cex=1, 
              colramp=colorRampPalette(mycolramp),xlab="DMSO D33 K27",nbin=1000,
              ylab="UNC 5FU D33 K27",bandwidth=c(0.01,0.01), xlim =c(0,13), ylim=c(0,13)) 
abline(a = 0, b=1, lty=2, lw = 2)
dev.off()

png(file.path(outputdir,"log2_Count_GSKJ4_vs_DMSO.png"), width = 1500, height = 1500, res = 300)
smoothScatter(x = log2(norm_mat[,"MM468BC_DMSO_D33_K27"]),y = log2(norm_mat[,"MM468BC_GSKJ4_D33_K27"]),cex=1, 
              colramp=colorRampPalette(mycolramp),xlab="DMSO D33 K27",nbin=1000,
              ylab="UNC GSKJ4 D33 K27",bandwidth=c(0.01,0.01), xlim =c(0,13), ylim=c(0,13)) 
abline(a = 0, b=1, lty=2, lw = 2)
dev.off()

png(file.path(outputdir,"log2_Count_5FU_vs_DMSO.png"), width = 1500, height = 1500, res = 300)
smoothScatter(x = log2(norm_mat[,"MM468BC_DMSO_D33_K27"]),y = log2(norm_mat[,"MM468BC_5FU_D33_K27"]),cex=1, 
              colramp=colorRampPalette(mycolramp),xlab="DMSO D33 K27",nbin=1000,
              ylab="UNC 5FU D33 K27",bandwidth=c(0.01,0.01), xlim =c(0,13), ylim=c(0,13)) 
abline(a = 0, b=1, lty=2, lw = 2)
dev.off()


norm_mat = mat
for(i in 1:ncol(norm_mat)){
  norm_mat[,i] = 10^6 * (norm_mat[,i] / sum(norm_mat[,i]))
}

mycolramp <- c("white",viridis(n=4))
png(file.path(outputdir,"log2_Count_UNC_vs_DMSO_libsize_only.png"), width = 1500, height = 1500, res = 300)
smoothScatter(x = log2(norm_mat[,"MM468BC_DMSO_D33_K27"]),y = log2(norm_mat[,"MM468BC_UNC_D33_K27"]),cex=1, 
              colramp=colorRampPalette(mycolramp),xlab="DMSO D33 K27",nbin=1000,
              ylab="UNC D33 K27",bandwidth=c(0.01,0.01), xlim =c(0,13), ylim=c(0,13)) 
abline(a = 0, b=1, lty=2, lw = 2)
dev.off()

png(file.path(outputdir,"log2_Count_UNC_5FU_vs_DMSO_libsize_only.png"), width = 1500, height = 1500, res = 300)
smoothScatter(x = log2(norm_mat[,"MM468BC_DMSO_D33_K27"]),y = log2(norm_mat[,"MM468BC_UNC_5FU_D33_K27"]),cex=1, 
              colramp=colorRampPalette(mycolramp),xlab="DMSO D33 K27",nbin=1000,
              ylab="UNC 5FU D33 K27",bandwidth=c(0.01,0.01), xlim =c(0,13), ylim=c(0,13)) 
abline(a = 0, b=1, lty=2, lw = 2)
dev.off()

png(file.path(outputdir,"log2_Count_GSKJ4_vs_DMSO_libsize_only.png"), width = 1500, height = 1500, res = 300)
smoothScatter(x = log2(norm_mat[,"MM468BC_DMSO_D33_K27"]),y = log2(norm_mat[,"MM468BC_GSKJ4_D33_K27"]),cex=1, 
              colramp=colorRampPalette(mycolramp),xlab="DMSO D33 K27",nbin=1000,
              ylab="UNC GSKJ4 D33 K27",bandwidth=c(0.01,0.01), xlim =c(0,13), ylim=c(0,13)) 
abline(a = 0, b=1, lty=2, lw = 2)
dev.off()

png(file.path(outputdir,"log2_Count_5FU_vs_DMSO_libsize_only.png"), width = 1500, height = 1500, res = 300)
smoothScatter(x = log2(norm_mat[,"MM468BC_DMSO_D33_K27"]),y = log2(norm_mat[,"MM468BC_5FU_D33_K27"]),cex=1, 
              colramp=colorRampPalette(mycolramp),xlab="DMSO D33 K27",nbin=1000,
              ylab="UNC 5FU D33 K27",bandwidth=c(0.01,0.01), xlim =c(0,13), ylim=c(0,13)) 
abline(a = 0, b=1, lty=2, lw = 2)
dev.off()