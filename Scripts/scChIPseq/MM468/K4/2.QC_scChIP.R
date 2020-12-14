library(here)
source(file.path(here(),"Scripts","functions.R"))
source(file.path(here(),"Scripts","global_var.R"))

maindir = here()
dataset_name = "MM468_H3K4me3_10k_TSS"
ref_genome ="hg38"

datadir = file.path(maindir, "input", "scChIPseq", "MM468")
outdir = file.path(maindir, "output", "scChIPseq", "MM468",dataset_name)
if(!dir.exists(outdir)) dir.create(outdir)
plotdir = file.path(outdir, "Plots"); if(!dir.exists(plotdir)) dir.create(plotdir)
plotdir_qc = file.path(plotdir, "QC"); if(!dir.exists(plotdir_qc)) dir.create(plotdir_qc)

ChromSCape_analyses = file.path(maindir, "output", "scChIPseq","ChromSCape_analyses")
ChromSCape_directory = file.path(ChromSCape_analyses,dataset_name)

# Load data K4 - PEAKS
tmp = tempdir()
unzip(file.path(datadir,"Count_Matrices", "MM468_K4_peaks_1000.zip"),
      exdir = tmp)
list.files(tmp,pattern = ".*.tsv")
# Reading in matrices
out <- import_scExp(c(file.path(tmp,"MM468_DMSO6_D0_K4.tsv"),
                      file.path(tmp,"MM468_5FU6_D60_K4.tsv")))

# Save raw
datamatrix = out$datamatrix
annot_raw = out$annot_raw
load(file.path(ChromSCape_directory,"Filtering_Normalize_Reduce","MM468_H3K4me3_10k_TSS_1600_1_95_uncorrected.RData"))

distribution = list()
for(file in list.files(file.path(datadir,"Raw_Counts","K4"),pattern = "*.count$")) {
    samp = gsub(".count","",file)
    distribution[[samp]] = read.csv(file.path(datadir,"Raw_Counts","K4",file),header=F,sep = " ",col.names = c("barcode","reads_total"))
    distribution[[samp]]$sample_id = rep(samp,nrow(distribution[[samp]]))
}

metadata = data.frame(sample_id = c(
    "MM468_DMSO6_D0_K4", "MM468_5FU6_D60_K4"),
    sample_id_color = c("#dfdfdfff", "#118675ff"))

# Calculate FrIP per cell for MM468
annot_raw$in_peak = Matrix::colSums(datamatrix)

annot_raw$total_reads = 0

for(i in names(distribution)){
    annot_raw$total_reads[which(annot_raw$sample_id == i)] = 
        distribution[[i]]$reads_total[match(
            annot_raw$barcode[which(annot_raw$sample_id == i)],distribution[[i]]$barcode)]
}

annot_raw = annot_raw %>% mutate(FRiP = in_peak / total_reads)
annot_raw$sample_id = factor(annot_raw$sample_id, levels=metadata$sample_id)

png(file.path(plotdir_qc, "FRiP_violin_K4.png"),height=1000,width=1500,res=300)
ggplot(annot_raw) + geom_violin(aes(y=FRiP, x=sample_id,fill=sample_id)) +
    theme_classic() + theme(axis.text.x = element_text(angle=90,vjust=1)) + ylim(c(0,1)) +
    scale_fill_manual(values = metadata$sample_id_color)
dev.off()

png(file.path(plotdir_qc, "FRiP_violin_K4_wolegend.png"),height=1500,width=1500,res=300)
ggplot(annot_raw) + geom_violin(aes(y=FRiP, x=sample_id,fill=sample_id)) +
    theme_classic() + theme(text = element_blank()) + ylim(c(0,1)) +
    scale_fill_manual(values = metadata$sample_id_color)
dev.off()

# UMAP
annot = as.data.frame(colData(scExp))
annot$FrIP = annot_raw$FRiP[match(annot$cell_id, annot_raw$cell_id)]
anocol = geco.unsupervised::geco.annotToCol4(annot[,c("total_counts","FrIP","sample_id")],
                                             annot[,c("total_counts","FrIP","sample_id")],
                                             categCol = NULL,
                                             plotLegend = T,
                                             plotLegendFile = file.path(plotdir_qc,"legend_FRiP.pdf"))

png(file.path(plotdir_qc,"UMAP_FRiP_K27.png"), height=1350,width=1200,res=300)
plot(reducedDim(scExp,"UMAP")[,c(1,2)], col=alpha(anocol[,"FrIP"],0.3),pch=20, cex=0.6,main=paste0("UMAP FRiP MM468"))
dev.off()

