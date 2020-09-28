setwd("~/Desktop/bwg/")
resDir <- "~/Desktop/bwg/"
load("~/Documents/bioinfo/Reference_Transcriptome/hg38/Gencode_hg38_v25.RData")
library(Sushi)
options(stringsAsFactors = F)

###Load data#########
C1 <- read.table("MM468_5_J67_DMSO_K27.bdg",header=F) #MM468_DMSO2
C2 <- read.table("MM468_5_J67_5FUR_K27.bdg",header=F) #MM468_5FUR2
#C3 <- read.table("test2.bdg",header=F)

sample1 <- "DMSO5_d67"
sample2 <- "5FUR5_d67"

col1 <- "#BDBDBD"
col2 <- "#9CCC65"


###Load HUMAN annotation#########
#subGenebed <-  GencodeByGene[GencodeByGene$Gene_biotype %in% c("protein_coding","lincRNA","antisense"),]
subGenebed <-  GencodeByGene[GencodeByGene$Gene_biotype %in% c("protein_coding"),]
genebed <- data.frame(chrom=subGenebed$Chr,start=subGenebed$Start,stop=subGenebed$End,gene=subGenebed$Gene_name,score=".",strand=subGenebed$Strand)
#correct in Gencode v25 and current (v33) TGFB1 coordinnates are not in line with other databases (end of gene)
genebed$start[genebed$gene=="TGFB1"] <- 41330323

###Load MOUSE annotation#########
#mouse_annot <- read.table("/Volumes/Data/OneDrive - HifiBio/Paper_scChIP_transfered_dropbox/Kevin/Annotation/Gencode_mm10_vM8_Genes.bed")
#colnames(mouse_annot) <- c("Chr","Start","End","score1","score2","Strand","Gene_name")
#genebed <- data.frame(chrom=mouse_annot$Chr,start=mouse_annot$Start,stop=mouse_annot$End,gene=mouse_annot$Gene_name,score=".",strand=mouse_annot$Strand)


###Select region#########
region="BMP4"
samples="hg38_MM468_persister"
chrom="chr14"
chromstart=53941000
chromend=53965000
C1_sel <- C1[C1[,1]==chrom,]
C2_sel <- C2[C2[,1]==chrom,]
#C3_sel <- C3[C3[,1]==chrom,]
gl.sub <- genebed[ which(genebed[,"chrom"] == chrom ),]

max_y <- 2

####ChIP-seq plots#######

####Plots with Peaks####

#png(file.path(resDir,paste0(samples,"_",region,".png")),width=1200,height=1400,res=300)
pdf(file.path(resDir,paste0(samples,"_",region,".pdf")),width=6,height=6)

par(mfrow=c(10,1))
#layout(m)
par(cex=0.5)
par(mar = c(1, 6, 3, 1), oma = c(1, 1, 1, 1))

plotBedgraph(C1_sel,chrom, chromstart,chromend,addscale = TRUE,ylab=sample1,color=col1,cex.lab=1,cex.main=2.1,range=c(0,max_y))
plotBedgraph(C2_sel,chrom, chromstart,chromend,addscale = TRUE,ylab=sample2,color=col2,cex.lab=1,cex.main=2.1,range=c(0,max_y))
#plotBedgraph(C3_sel,chrom, chromstart,chromend,addscale=TRUE,ylab="C3",color="#bfbebe",cex.lab=1.3,cex.main=2.1,range=c(0,max_y))


#labelgenome(chrom,chromstart,chromend,n=4,scale="Mb",side=3)
par(mar = c(1, 6, 1, 1),xpd=NA)
plotGenes(gl.sub, chrom,chromstart,chromend,bentline=F,plotgenetype = "arrow",labeltext = T,fontsize=1,labeloffset = 0.6,bheight=0.07, col="grey",)
geneinfo=gl.sub
#plot.new()
par(mar = c(1, 6, 1, 1))
labelgenome(chrom,chromstart,chromend,scale="Kb",cex.axis=1.5,n=3,scalecex = 0.8, chromcex= 0.8)

plot.new()
plot.new()
legend("bottom",legend=c(sample1,sample2),fill=c(col1,col2),text.font=15,cex=1.7,bty="n")

dev.off()


