#analysis fold-change MM468
library(here)

maindir= here()
resdir <- file.path(maindir,"output","scRNAseq","MM468")
plotdir <- file.path(resdir,"CompareFC");if(!file.exists(plotdir)){dir.create(plotdir)}

source(file.path(maindir,"Scripts","global_var.R"))
mycolramp <- c("white",viridis(n=4))

# Persister log2FC
persister = new.env()
load(file.path(resdir, "Persister","Supervised","RData","Supervised_res_object_edgeR.Rdata"), persister)
# UNC log2FC
UNC = new.env()
load(file.path(resdir, "UNC","Supervised","RData","Supervised_res_object_edgeR.Rdata"), UNC)

common = intersect(persister$my.res$Symbol,UNC$my.res$Symbol)
persister$my.res = persister$my.res[match(common,persister$my.res$Symbol),]
UNC$my.res = UNC$my.res[match(common,UNC$my.res$Symbol),]

png(file.path(plotdir,"log2FC_persister_vs_UNC.png"), width = 2000, height = 2000, res = 600)
smoothScatter(UNC$my.res$log2FC.UNC, persister$my.res$log2FC.C2_pers,cex=1, 
              colramp=colorRampPalette(mycolramp),xlab="log2FC UNC vs initial",nbin=300,
              ylab="log2FC 5FU6 persister vs initial",bandwidth=c(0.01,0.01),xlim=c(-5,5),ylim=c(-5,5)) 
dev.off()

png(file.path(plotdir,"log2FC_resistant_vs_UNC.png"), width = 2000, height = 2000, res = 600)
smoothScatter(UNC$my.res$log2FC.UNC, persister$my.res$log2FC.MM468_5FU6_day214,cex=1, 
              colramp=colorRampPalette(mycolramp),xlab="log2FC UNC vs initial",nbin=300,
              ylab="log2FC 5FU6 resistant vs initial",bandwidth=c(0.0001,0.0001),xlim=c(-5,5),ylim=c(-5,5)) 
dev.off()

# Load K27 10k TSS
annot10k_K27 <- read.table(
    unz(file.path(
        maindir,"annotation","gencode.v34.transcripts10k_K27.zip"),
        "gencode.v34.transcripts10k_K27.bed"), sep="\t")[,-c(6,9)]
colnames(annot10k_K27) <- c("chr","start","end","transcripts","gene","log2FC_K27","qvalue_K27","K27_status")

annot10k_K27_ovlp = annot10k_K27[which(annot10k_K27$K27_status !="Not Overlapping"),]
annot10k_K27_byGene = annot10k_K27_ovlp %>% group_by(gene) %>%  slice_max(order_by = abs(log2FC_K27), n = 1) # Select only the transcript associated with top differential K27 log2FC 



tab = tab[,c("Symbol","log2FC.C2_pers","qval.C2_pers","log2FC.UNC","qval.UNC")]

tab$Persister ="Non Signif."
tab$Persister[which(tab$log2FC.C2_pers > log2(3) & tab$qval.C2_pers < 0.01)] = "Over"
tab$Persister[which(tab$log2FC.C2_pers < -log2(3) & tab$qval.C2_pers < 0.01)] = "Under"

tab$UNC ="Non Signif."
tab$UNC[which(tab$log2FC.UNC > log2(3) & tab$qval.UNC < 0.01)] = "Over"
tab$UNC[which(tab$log2FC.UNC < -log2(3) & tab$qval.UNC < 0.01)] = "Under"


tab$K27 = "Not differential"
annot10k_K27 = annot10k_K27[which(annot10k_K27$gene %in% tab$Symbol),]
sel = which(tab$Symbol %in% annot10k_K27$gene)
tab$K27[sel] = annot10k_K27$K27_status[match(tab$Symbol[sel],annot10k_K27$gene)]
tab = tab %>% mutate(K27 = ifelse(K27=="Not differential", "Overlapping Not Differential",K27))
tab = tab %>% mutate(K27 = ifelse(K27=="Depleted", "Overlapping Depleted",K27))
tab = tab %>% mutate(K27 = ifelse(K27=="Enriched", "Overlapping Enriched",K27))

WriteXLS::WriteXLS(tab,ExcelFileName = file.path(plotdir,"Table_Persister_UNC_K27.xls"))
