library(here)
source(file.path(here(),"Scripts","functions.R"))
source(file.path(here(),"Scripts","global_var.R"))

maindir = here()

dataset_name = "MM468_H3K27me3_peaks"
ref_genome ="hg38"

datadir = file.path(maindir, "input", "scChIPseq", "MM468")
outdir = file.path(maindir, "output", "scChIPseq", "MM468",dataset_name)
if(!dir.exists(outdir)) dir.create(outdir)
plotdir = file.path(outdir, "Plots"); if(!dir.exists(plotdir)) dir.create(plotdir)
plotdir_qc = file.path(plotdir, "QC"); if(!dir.exists(plotdir_qc)) dir.create(plotdir_qc)

ChromSCape_analyses = file.path(maindir, "output", "scChIPseq","ChromSCape_analyses")
ChromSCape_directory = file.path(ChromSCape_analyses,dataset_name)

# Load data K27
load(file.path(ChromSCape_directory,"scChIP_raw.RData"))
load(file.path(ChromSCape_directory,"correlation_clustering","MM468_H3K27me3_peaks_3000_0_95_uncorrected.RData"))

distribution = list()
for(file in list.files(file.path(datadir,"Raw_Counts","K27"),pattern = "*.count$")) {
    samp = gsub(".count","",file)
    distribution[[samp]] = read.csv(file.path(datadir,"Raw_Counts","K27",file),header=F,sep = " ",col.names = c("barcode","reads_total"))
    distribution[[samp]]$sample_id = rep(samp,nrow(distribution[[samp]]))
}

metadata_K27 = data.frame(sample_id = c(
    "MM468_DMSO1_day60", "MM468_DMSO3_day77", "MM468_DMSO5_day131",
    "MM468_5FU1_day33", "MM468_5FU2_day67",
    "MM468_5FU6_day131", "MM468_5FU3_day147", "MM468_5FU2_day171"),
    sample_id_color = c("#dfdfdfff", "#999999ff","#363636",
                        "#118675ff", "#8cc453ff",
                        "#ff5722ff", "#feb40fff", "#fd8508ff"))
# Calculate FrIP per cell for MM468
annot_raw$in_peak = Matrix::colSums(datamatrix)

annot_raw$total_reads = 0

for(i in names(distribution)){
    annot_raw$total_reads[which(annot_raw$sample_id == i)] = 
        distribution[[i]]$reads_total[match(
            annot_raw$barcode[which(annot_raw$sample_id == i)],distribution[[i]]$barcode)]
}

annot_raw = annot_raw %>% mutate(FRiP = in_peak / total_reads)
annot_raw$sample_id = factor(annot_raw$sample_id, levels=metadata_K27$sample_id)

png(file.path(plotdir_qc, "FRiP_violin_K27.png"),height=1000,width=1500,res=300)
ggplot(annot_raw) + geom_violin(aes(y=FRiP, x=sample_id,fill=sample_id)) +
    theme_classic() + theme(axis.text.x = element_text(angle=90,vjust=1)) + ylim(c(0,1)) +
    scale_fill_manual(values = metadata_K27$sample_id_color)
dev.off()

png(file.path(plotdir_qc, "FRiP_violin_K27_wolegend.png"),height=1500,width=1500,res=300)
ggplot(annot_raw) + geom_violin(aes(y=FRiP, x=sample_id,fill=sample_id)) +
    theme_classic() + theme(text = element_blank()) + ylim(c(0,1)) +
    scale_fill_manual(values = metadata_K27$sample_id_color)
dev.off()

# UMAP
annot = as.data.frame(colData(scExp_cf))
umap_res = as.data.frame(reducedDim(scExp_cf,"UMAP"))
annot$FrIP = annot_raw$FRiP[match(annot$cell_id, annot_raw$cell_id)]
annot = annot[-which(is.infinite(annot$FrIP)),]
anocol = geco.unsupervised::geco.annotToCol4(annot[,c("total_counts","FrIP","sample_id")],
                                             annot[,c("total_counts","FrIP","sample_id")],
                                             plotLegend = T,
                                             plotLegendFile = file.path(plotdir_qc,"legend_FRiP.pdf"),
                                             scale_q = "inferno")

png(file.path(plotdir_qc,"UMAP_FRiP_K27.png"), height=1350,width=1200,res=300)
plot(umap_res[,c(1,2)], col=alpha(anocol[,"FrIP"],0.3),pch=20, cex=0.6,main=paste0("UMAP FRiP MM468"))
dev.off()



