library(here)
maindir= here()
source(file.path(maindir, "Scripts", "global_var.R"))
source(file.path(maindir, "Scripts", "functions.R"))

inputDir_bulk <- file.path(maindir,"input","bulk_ChIPseq","MM468","BigWigs")
input_ratios <- file.path(maindir,"output","bulk_ChIPseq","MM468","ChromatinIndexing","ratio_ip_input_K27.csv")
resDir <- file.path(maindir,"output","bulk_ChIPseq","Snapshots")
load(file.path(maindir,"annotation","Gencode_hg38_v25.RData"))

print(list.files(inputDir_bulk))

###Load data#########
normalizing_ratios = read.csv(input_ratios)
rownames(normalizing_ratios) = normalizing_ratios$Sample

## 5FU & DMSO6 Single Cell K27
cum_5FU_D33_K27_CI <- rtracklayer::import(file.path(inputDir_bulk, "MM468_5FU_day33_F02_H3K27me3.bw"))
cum_DMSO_D33_K27_CI <- rtracklayer::import(file.path(inputDir_bulk, "MM468_DMSO_day33_E02_H3K27me3.bw"))

## 5FU & DMSO Bulk K27
cum_GSKJ4_D33_K27_CI <- rtracklayer::import(file.path(inputDir_bulk, "MM468_GSKJ4_day33_D02_H3K27me3.bw"))
cum_UNC_D33_K27_CI <- rtracklayer::import(file.path(inputDir_bulk, "MM468_UNC_day33_A02_H3K27me3.bw"))
cum_UNC_5FU_D33_K27_CI <- rtracklayer::import(file.path(inputDir_bulk, "MM468_UNC_5FU_day33_C02_H3K27me3.bw"))

# Normalize by the corresponding ratios IP / Input
cum_5FU_D33_K27_CI$score = cum_5FU_D33_K27_CI$score * normalizing_ratios["5FU_D33_K27","value"]
cum_DMSO_D33_K27_CI$score = cum_DMSO_D33_K27_CI$score * normalizing_ratios["DMSO_D33_K27","value"]
cum_GSKJ4_D33_K27_CI$score = cum_GSKJ4_D33_K27_CI$score * normalizing_ratios["GSKJ4_D33_K27","value"]
cum_UNC_D33_K27_CI$score = cum_UNC_D33_K27_CI$score * normalizing_ratios["UNC_D33_K27","value"]
cum_UNC_5FU_D33_K27_CI$score = cum_UNC_5FU_D33_K27_CI$score * normalizing_ratios["UNC_5FU_D33_K27","value"]
gc()

###Load HUMAN annotation#########
subGenebed <-  GencodeByGene[GencodeByGene$Gene_biotype %in% c("protein_coding"),]
genebed <- data.frame(chrom=subGenebed$Chr,start=subGenebed$Start,stop=subGenebed$End,gene=subGenebed$Gene_name,score=".",strand=subGenebed$Strand)
genebed$chrom = as.character(genebed$chrom)
genebed$gene = as.character(genebed$gene)

regions = list(
    c("COL4A1_2","chr13",110153364,110500000),
    c("COL8A1","chr3",99626127,99800686),
    c("TGFB1","chr19",41330000,41360000),
    c("TGFBR3","chr1",91850000,91950000),
    c("FOXBP1","chr15",59976664,60036032),
    c("KRTKs","chr17",41429737,41604660),
    c("NNMT","chr11",114264830,114343463),
    c("FRMD4A","chr10",13429037,14675173),
    c("KCNK9","chr8",139567430,139736399),
    c("KLKs","chr19",51010000,51030000),
    c("LAMB1","chr7",107939129,108020000),
    c("BMP4","chr14",53942693,53962060),
    c("KRT14","chr17",41580279,41588895),
    c("KRT16","chr17",41607779,41614827),
    c("KRT6A","chr12",52485174,52495397),
    c("KRT6B","chr12",52444651,52454126),
    c("ANXA3","chr4",78549588,78612451),
    c("S100P","chr4",6691839,6699170),
    c("INHBA","chr7",41683114,41705108),
    c("INHBB","chr2",120344142,120353808),
    c("TAGLN","chr11",117197324,117206792),
    c("EPHB6","chr7",142853014,142873093),
    c("CDH11","chr16",65090000,65140000),
    c("MIF","chr22",23892948,23896644),
    c("CST6","chr11",66009991,66015505),
    c("FOXQ1","chr6",1300000,1330000),
    c("BMP6","chr6",7724099,7883728),
    c("NEXN","chr1",77886515,77945893),
    c("KRTDAP","chr19",35485324,35502531),
    c("GDA","chr9",72022248,72344485),
    c("TGFBR2","chr3",30604502,30696141),
    c("TGFB2","chr1",218343284,218446619),
    c("KRT17","chr17",41617440,41626630),
    c("KRT5","chr12",52512575,52522459),
    c("KRT8","chr12",52895187,52951866),
    c("TGFBR1","chr9",99103089,99156191),
    c("FOSL1","chr11",65890136,65902526),
    c("LAMB3","chr1",209612873,209654475),
    c("ABCC4","chr13",95272000,95348900),
    c("ABCC4_large","chr13",95210000,95335000),
    c("ABCA5","chr17",69296235,69346469),
    c("LBH","chr2",30220798,30257330)
)


samples="MM468_ChromatinIndexing_poolK27"
if(!dir.exists(file.path(resDir,samples))) dir.create(file.path(resDir,samples))

for(i in 1:length(regions)){
    #Select region of interest
    region=regions[[i]][1]
    chrom=regions[[i]][2]
    chromstart=as.numeric(regions[[i]][3])
    chromend=as.numeric(regions[[i]][4])
    
    roi = GRanges(chrom,ranges = IRanges(chromstart,chromend), Gene=region)
    
    cat("Doing plot for ",region,"",chrom,":",chromstart,"-",chromend,".\n")
    gl.sub <- genebed[ which(genebed[,"chrom"] == chrom ),]

    ## Chromatin indexing samples
    cum_5FU_D33_K27_CI_tmp = cum_5FU_D33_K27_CI[subjectHits(findOverlaps(roi,cum_5FU_D33_K27_CI)),]
    cum_DMSO_D33_K27_CI_tmp = cum_DMSO_D33_K27_CI[subjectHits(findOverlaps(roi,cum_DMSO_D33_K27_CI)),]
    cum_GSKJ4_D33_K27_CI_tmp = cum_GSKJ4_D33_K27_CI[subjectHits(findOverlaps(roi,cum_GSKJ4_D33_K27_CI)),]
    cum_UNC_D33_K27_CI_tmp = cum_UNC_D33_K27_CI[subjectHits(findOverlaps(roi,cum_UNC_D33_K27_CI)),]
    cum_UNC_5FU_D33_K27_CI_tmp = cum_UNC_5FU_D33_K27_CI[subjectHits(findOverlaps(roi,cum_UNC_5FU_D33_K27_CI)),]
    
    max = round(max(c(cum_5FU_D33_K27_CI_tmp$score,
                      cum_DMSO_D33_K27_CI_tmp$score,
                      cum_GSKJ4_D33_K27_CI_tmp$score,
                      cum_UNC_D33_K27_CI_tmp$score,
                      cum_UNC_5FU_D33_K27_CI_tmp$score
                      )
                    ),1)
    
    pdf(file.path(file.path(resDir,samples),paste0(samples,"_",region,".pdf")),width=8,height=6)
    layout.matrix <- matrix(c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,7), ncol = 1)
    
    layout(mat = layout.matrix,
           heights = c(1), # Heights of the two rows
           widths = c(1)) # Widths of the two columns
    
    par(cex=0.5)
    par(mar = c(0.75, 6, 0, 1), oma = c(1, 1, 1, 1))
    
    #bulk K27
    plotBedgraph(as.data.frame(cum_DMSO_D33_K27_CI_tmp)[,c(1,2,3,6)],
                 chrom, chromstart,chromend,
                 range = c(0,max),
                 addscale = TRUE,
                 ylab="DMSO_K27",
                 color="#afafafff",cex.lab=1,cex.main=2.1)
    
    plotBedgraph(as.data.frame(cum_5FU_D33_K27_CI_tmp)[,c(1,2,3,6)],
                 chrom, chromstart,chromend,
                 range = c(0, max),
                 addscale = TRUE,
                 ylab="5FU_K27",
                 color="#118675ff",cex.lab=1,cex.main=2.1)
    
    plotBedgraph(as.data.frame(cum_UNC_D33_K27_CI_tmp)[,c(1,2,3,6)],
                 chrom, chromstart,chromend,
                 range = c(0,max),
                 addscale = TRUE,
                 ylab="UNC_K27",
                 color="#41c8dcff",cex.lab=1,cex.main=2.1)
    
    plotBedgraph(as.data.frame(cum_UNC_5FU_D33_K27_CI_tmp)[,c(1,2,3,6)],
                 chrom, chromstart,chromend,
                 range = c(0, max),
                 addscale = TRUE,
                 ylab="UNC_5FU_D33_K27",
                 color="#2686f0ff",cex.lab=1,cex.main=2.1)
    
    plotBedgraph(as.data.frame(cum_GSKJ4_D33_K27_CI_tmp)[,c(1,2,3,6)],
                 chrom, chromstart,chromend,
                 range = c(0, max),
                 addscale = TRUE,
                 ylab="GSKJ4_K27",
                 color="#a5151cff",cex.lab=1,cex.main=2.1)
    
    par(mar = c(1, 6, 1, 1),xpd=NA)
    
    plotGenes(gl.sub, chrom,chromstart,chromend,
              bentline=F,plotgenetype = "arrow",
              labeltext = T,labelat = "start",fontsize=1,labeloffset = 0.4,bheight=0.07)
    geneinfo=gl.sub
    
    par(mar = c(1, 6, 1, 1))
    labelgenome(chrom,chromstart,chromend,n=4,scale="Mb",cex.axis=1.5)
    
    dev.off()
}

library(here)
maindir= here()
source(file.path(maindir, "Scripts", "global_var.R"))
source(file.path(maindir, "Scripts", "functions.R"))

inputDir_bulk <- file.path(maindir,"input","bulk_ChIPseq","MM468","BigWigs")
input_ratios <- file.path(maindir,"output","bulk_ChIPseq","MM468","ChromatinIndexing","ratio_ip_input_K27.csv")
resDir <- file.path(maindir,"output","bulk_ChIPseq","Snapshots")
load(file.path(maindir,"annotation","Gencode_hg38_v25.RData"))

print(list.files(inputDir_bulk))

###Load data#########
normalizing_ratios = read.csv(input_ratios)
rownames(normalizing_ratios) = normalizing_ratios$Sample

## 5FU & DMSO6 Single Cell K27
cum_5FU_D33_K27_CI <- rtracklayer::import(file.path(inputDir_bulk, "MM468_5FU_day33_F02_H3K27me3.bw"))
cum_DMSO_D33_K27_CI <- rtracklayer::import(file.path(inputDir_bulk, "MM468_DMSO_day33_E02_H3K27me3.bw"))

## 5FU & DMSO Bulk K27
cum_GSKJ4_D33_K27_CI <- rtracklayer::import(file.path(inputDir_bulk, "MM468_GSKJ4_day33_D02_H3K27me3.bw"))
cum_UNC_D33_K27_CI <- rtracklayer::import(file.path(inputDir_bulk, "MM468_UNC_day33_A02_H3K27me3.bw"))
cum_UNC_5FU_D33_K27_CI <- rtracklayer::import(file.path(inputDir_bulk, "MM468_UNC_5FU_day33_C02_H3K27me3.bw"))

# Normalize by the corresponding ratios IP / Input
cum_5FU_D33_K27_CI$score = cum_5FU_D33_K27_CI$score * normalizing_ratios["5FU_D33_K27","value"]
cum_DMSO_D33_K27_CI$score = cum_DMSO_D33_K27_CI$score * normalizing_ratios["DMSO_D33_K27","value"]
cum_GSKJ4_D33_K27_CI$score = cum_GSKJ4_D33_K27_CI$score * normalizing_ratios["GSKJ4_D33_K27","value"]
cum_UNC_D33_K27_CI$score = cum_UNC_D33_K27_CI$score * normalizing_ratios["UNC_D33_K27","value"]
cum_UNC_5FU_D33_K27_CI$score = cum_UNC_5FU_D33_K27_CI$score * normalizing_ratios["UNC_5FU_D33_K27","value"]
gc()

###Load HUMAN annotation#########
subGenebed <-  GencodeByGene[GencodeByGene$Gene_biotype %in% c("protein_coding"),]
genebed <- data.frame(chrom=subGenebed$Chr,start=subGenebed$Start,stop=subGenebed$End,gene=subGenebed$Gene_name,score=".",strand=subGenebed$Strand)
genebed$chrom = as.character(genebed$chrom)
genebed$gene = as.character(genebed$gene)

regions = list(
    c("COL4A1_2","chr13",110153364,110500000),
    c("COL8A1","chr3",99626127,99800686),
    c("TGFB1","chr19",41330000,41360000),
    c("TGFBR3","chr1",91850000,91950000),
    c("FOXBP1","chr15",59976664,60036032),
    c("KRTKs","chr17",41429737,41604660),
    c("NNMT","chr11",114264830,114343463),
    c("FRMD4A","chr10",13429037,14675173),
    c("KCNK9","chr8",139567430,139736399),
    c("KLKs","chr19",51010000,51030000),
    c("LAMB1","chr7",107939129,108020000),
    c("BMP4","chr14",53942693,53962060),
    c("KRT14","chr17",41580279,41588895),
    c("KRT16","chr17",41607779,41614827),
    c("KRT6A","chr12",52485174,52495397),
    c("KRT6B","chr12",52444651,52454126),
    c("ANXA3","chr4",78549588,78612451),
    c("S100P","chr4",6691839,6699170),
    c("INHBA","chr7",41683114,41705108),
    c("INHBB","chr2",120344142,120353808),
    c("TAGLN","chr11",117197324,117206792),
    c("EPHB6","chr7",142853014,142873093),
    c("CDH11","chr16",65090000,65140000),
    c("MIF","chr22",23892948,23896644),
    c("CST6","chr11",66009991,66015505),
    c("FOXQ1","chr6",1300000,1330000),
    c("BMP6","chr6",7724099,7883728),
    c("NEXN","chr1",77886515,77945893),
    c("KRTDAP","chr19",35485324,35502531),
    c("GDA","chr9",72022248,72344485),
    c("TGFBR2","chr3",30604502,30696141),
    c("TGFB2","chr1",218343284,218446619),
    c("KRT17","chr17",41617440,41626630),
    c("KRT5","chr12",52512575,52522459),
    c("KRT8","chr12",52895187,52951866),
    c("TGFBR1","chr9",99103089,99156191),
    c("FOSL1","chr11",65890136,65902526),
    c("LAMB3","chr1",209612873,209654475),
    c("ABCC4","chr13",95272000,95348900),
    c("ABCC4_large","chr13",95210000,95335000),
    c("ABCA5","chr17",69296235,69346469),
    c("LBH","chr2",30220798,30257330)
)


samples="MM468_ChromatinIndexing_poolK27"
if(!dir.exists(file.path(resDir,samples))) dir.create(file.path(resDir,samples))

for(i in 1:length(regions)){
    #Select region of interest
    region=regions[[i]][1]
    chrom=regions[[i]][2]
    chromstart=as.numeric(regions[[i]][3])
    chromend=as.numeric(regions[[i]][4])
    
    roi = GRanges(chrom,ranges = IRanges(chromstart,chromend), Gene=region)
    
    cat("Doing plot for ",region,"",chrom,":",chromstart,"-",chromend,".\n")
    gl.sub <- genebed[ which(genebed[,"chrom"] == chrom ),]

    ## Chromatin indexing samples
    cum_5FU_D33_K27_CI_tmp = cum_5FU_D33_K27_CI[subjectHits(findOverlaps(roi,cum_5FU_D33_K27_CI)),]
    cum_DMSO_D33_K27_CI_tmp = cum_DMSO_D33_K27_CI[subjectHits(findOverlaps(roi,cum_DMSO_D33_K27_CI)),]
    cum_GSKJ4_D33_K27_CI_tmp = cum_GSKJ4_D33_K27_CI[subjectHits(findOverlaps(roi,cum_GSKJ4_D33_K27_CI)),]
    cum_UNC_D33_K27_CI_tmp = cum_UNC_D33_K27_CI[subjectHits(findOverlaps(roi,cum_UNC_D33_K27_CI)),]
    cum_UNC_5FU_D33_K27_CI_tmp = cum_UNC_5FU_D33_K27_CI[subjectHits(findOverlaps(roi,cum_UNC_5FU_D33_K27_CI)),]
    
    max = round(max(c(cum_5FU_D33_K27_CI_tmp$score,
                      cum_DMSO_D33_K27_CI_tmp$score,
                      cum_GSKJ4_D33_K27_CI_tmp$score,
                      cum_UNC_D33_K27_CI_tmp$score,
                      cum_UNC_5FU_D33_K27_CI_tmp$score
                      )
                    ),1)
    
    pdf(file.path(file.path(resDir,samples),paste0(samples,"_",region,".pdf")),width=8,height=6)
    layout.matrix <- matrix(c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,7), ncol = 1)
    
    layout(mat = layout.matrix,
           heights = c(1), # Heights of the two rows
           widths = c(1)) # Widths of the two columns
    
    par(cex=0.5)
    par(mar = c(0.75, 6, 0, 1), oma = c(1, 1, 1, 1))
    
    #bulk K27
    plotBedgraph(as.data.frame(cum_DMSO_D33_K27_CI_tmp)[,c(1,2,3,6)],
                 chrom, chromstart,chromend,
                 range = c(0,max),
                 addscale = TRUE,
                 ylab="DMSO_K27",
                 color="#afafafff",cex.lab=1,cex.main=2.1)
    
    plotBedgraph(as.data.frame(cum_5FU_D33_K27_CI_tmp)[,c(1,2,3,6)],
                 chrom, chromstart,chromend,
                 range = c(0, max),
                 addscale = TRUE,
                 ylab="5FU_K27",
                 color="#118675ff",cex.lab=1,cex.main=2.1)
    
    plotBedgraph(as.data.frame(cum_UNC_D33_K27_CI_tmp)[,c(1,2,3,6)],
                 chrom, chromstart,chromend,
                 range = c(0,max),
                 addscale = TRUE,
                 ylab="UNC_K27",
                 color="#41c8dcff",cex.lab=1,cex.main=2.1)
    
    plotBedgraph(as.data.frame(cum_UNC_5FU_D33_K27_CI_tmp)[,c(1,2,3,6)],
                 chrom, chromstart,chromend,
                 range = c(0, max),
                 addscale = TRUE,
                 ylab="UNC_5FU_D33_K27",
                 color="#2686f0ff",cex.lab=1,cex.main=2.1)
    
    plotBedgraph(as.data.frame(cum_GSKJ4_D33_K27_CI_tmp)[,c(1,2,3,6)],
                 chrom, chromstart,chromend,
                 range = c(0, max),
                 addscale = TRUE,
                 ylab="GSKJ4_K27",
                 color="#a5151cff",cex.lab=1,cex.main=2.1)
    
    par(mar = c(1, 6, 1, 1),xpd=NA)
    
    plotGenes(gl.sub, chrom,chromstart,chromend,
              bentline=F,plotgenetype = "arrow",
              labeltext = T,labelat = "start",fontsize=1,labeloffset = 0.4,bheight=0.07)
    geneinfo=gl.sub
    
    par(mar = c(1, 6, 1, 1))
    labelgenome(chrom,chromstart,chromend,n=4,scale="Mb",cex.axis=1.5)
    
    dev.off()
}

