library(here)

maindir= here()
source(file.path(maindir, "Scripts", "global_var.R"))
source(file.path(maindir, "Scripts", "functions.R"))

inputDir_sc <- file.path(maindir,"input","scChIPseq","MM468","BigWigs")
inputDir_bulk <- file.path(maindir,"input","bulk_ChIPseq","BigWigs","MM468")
resDir <- file.path(maindir,"output","bulk_ChIPseq","Snapshots")
load(file.path(maindir,"annotation","Gencode_hg38_v25.RData"))

input_ratios <- file.path(maindir,"output","bulk_ChIPseq","MM468","ChromatinIndexing","ratio_ip_input_K4.csv")
normalizing_ratios = read.csv(input_ratios)
rownames(normalizing_ratios) = normalizing_ratios$Sample

print(list.files(inputDir_sc))
print(list.files(inputDir_bulk))

###Load data#########
## 5FU & DMSO6 Single Cell K27
cum_5FU6_D33_K27 <- rtracklayer::import(file.path(inputDir_sc, "MM468_5FU6_day33_K27me3.bw"))
cum_DMSO3_D77_K27 <- rtracklayer::import(file.path(inputDir_sc, "MM468_DMSO3_day77_K27me3.bw"))

## 5FU & DMSO Bulk K27
cum_5FU5_day67_bulk_K27 <- rtracklayer::import(file.path(inputDir_bulk, "MM468_5_J67_5FUR_K27.bw"))
cum_DMSO5_day67_bulk_K27 <- rtracklayer::import(file.path(inputDir_bulk, "MM468_5_J67_DMSO_K27.bw"))

## 5FU & DMSO 6 SingleCell K4  
cum_5FU6_D60_K4 <- rtracklayer::import(file.path(inputDir_sc, "MM468_5FU6_day60_K4me3.bw"))
cum_DMSO6_day0_K4 <- rtracklayer::import(file.path(inputDir_sc, "MM468_DMSO6_day0_K4me3.bw"))

## 5FU & DMSO Bulk K4  
cum_5FU6_day33_bulk_K4 <- rtracklayer::import(file.path(inputDir_bulk, "MM468_5FU_D33_K4_E02.bw"))
cum_DMSO6_day33_bulk_K4 <- rtracklayer::import(file.path(inputDir_bulk, "MM468_DMSO_D33_K4_C02.bw"))

# correct by IP / Input ratio


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
    c("FOXQ1","chr6",1310440,1316758),
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
    c("LAMB3","chr1",209612873,209654475)
)

samples="MM468_5FU6_vs_DMSO6_sc_bulk_K27_K4"
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
    
    ## 5FU & DMSO sc K27
    cum_5FU6_D33_K27_tmp = cum_5FU6_D33_K27[subjectHits(findOverlaps(roi,cum_5FU6_D33_K27)),]
    cum_DMSO3_D77_K27_tmp = cum_DMSO3_D77_K27[subjectHits(findOverlaps(roi,cum_DMSO3_D77_K27)),]
    max1 = round(max(c(cum_5FU6_D33_K27_tmp$score,cum_DMSO3_D77_K27_tmp$score)),3)
    
    ## 5FU & DMSO Bulk K27
    cum_5FU5_day67_bulk_K27_tmp = cum_5FU5_day67_bulk_K27[subjectHits(findOverlaps(roi,cum_5FU5_day67_bulk_K27)),]
    cum_DMSO5_day67_bulk_K27_tmp = cum_DMSO5_day67_bulk_K27[subjectHits(findOverlaps(roi,cum_DMSO5_day67_bulk_K27)),]
    max2 = round(max(c(cum_DMSO5_day67_bulk_K27_tmp$score,cum_5FU5_day67_bulk_K27_tmp$score)),3)
    
    ## 5FU & DMSO 6 SingleCell K4  
    cum_5FU6_D60_K4_tmp = cum_5FU6_D60_K4[subjectHits(findOverlaps(roi,cum_5FU6_D60_K4)),]
    cum_DMSO6_day0_K4_tmp = cum_DMSO6_day0_K4[subjectHits(findOverlaps(roi,cum_DMSO6_day0_K4)),]
    max3 = round(max(c(cum_5FU6_D60_K4_tmp$score,cum_DMSO6_day0_K4_tmp$score)),3)
    
    ## 5FU & DMSO Bulk K4  
    cum_5FU6_day33_bulk_K4_tmp = cum_5FU6_day33_bulk_K4[subjectHits(findOverlaps(roi,cum_5FU6_day33_bulk_K4)),]
    cum_DMSO6_day33_bulk_K4_tmp = cum_DMSO6_day33_bulk_K4[subjectHits(findOverlaps(roi,cum_DMSO6_day33_bulk_K4)),]
    max4 = round(max(c(cum_5FU6_day33_bulk_K4_tmp$score,cum_DMSO6_day33_bulk_K4_tmp$score)),3)
    
    pdf(file.path(file.path(resDir,samples),paste0(samples,"_",region,".pdf")),width=8,height=6)
    layout.matrix <- matrix(c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6,7,7,7,8,8,8,9,9,10), ncol = 1)
    
    layout(mat = layout.matrix,
           heights = c(1), # Heights of the two rows
           widths = c(1)) # Widths of the two columns
    
    par(cex=0.5)
    par(mar = c(0.75, 6, 0, 1), oma = c(1, 1, 1, 1))
    
    #bulk K27
    plotBedgraph(as.data.frame(cum_DMSO5_day67_bulk_K27_tmp)[,c(1,2,3,6)],
                 chrom, chromstart,chromend,
                 range = c(0, max2),
                 addscale = TRUE,
                 ylab="bulk_DMSO_K27",
                 color="gray88",cex.lab=1,cex.main=2.1)
    
    
    plotBedgraph(as.data.frame(cum_5FU5_day67_bulk_K27_tmp)[,c(1,2,3,6)],
                 chrom, chromstart,chromend,
                 range = c(0,max2),
                 addscale = TRUE,
                 ylab="bulk_5FU_K27",
                 color="#009688",cex.lab=1,cex.main=2.1)
    
    
    #sc K27 c("#A6A6A6FF", "#00544CFF")
    plotBedgraph(as.data.frame(cum_DMSO3_D77_K27_tmp)[,c(1,2,3,6)],
                 chrom, chromstart,chromend,
                 range = c(0, max1),
                 addscale = TRUE,
                 ylab="sc_DMSO_K27",
                 color="#A6A6A6FF",cex.lab=1,cex.main=2.1)
    
    plotBedgraph(as.data.frame(cum_5FU6_D33_K27_tmp)[,c(1,2,3,6)],
                 chrom, chromstart,chromend,
                 range = c(0,max1),
                 addscale = TRUE,
                 ylab="sc_5FU_K27",
                 color="#00544CFF",cex.lab=1,cex.main=2.1)
    
    #Bulk K4
    plotBedgraph(as.data.frame(cum_DMSO6_day33_bulk_K4_tmp)[,c(1,2,3,6)],
                 chrom, chromstart,chromend,
                 range = c(0, max4),
                 addscale = TRUE,
                 ylab="bulk_DMSO_K4",
                 color="gray88",cex.lab=1,cex.main=2.1)
    
    plotBedgraph(as.data.frame(cum_5FU6_day33_bulk_K4_tmp)[,c(1,2,3,6)],
                 chrom, chromstart,chromend,
                 range = c(0,max4),
                 addscale = TRUE,
                 ylab="bulk_5FU_K4",
                 color="#009688",cex.lab=1,cex.main=2.1)
    
    #sc K4 c("#A6A6A6FF", "#00544CFF")
    plotBedgraph(as.data.frame(cum_DMSO6_day0_K4_tmp)[,c(1,2,3,6)],
                 chrom, chromstart,chromend,
                 range = c(0, max3),
                 addscale = TRUE,
                 ylab="sc_DMSO_K4",
                 color="#A6A6A6FF",cex.lab=1,cex.main=2.1)
    
    plotBedgraph(as.data.frame(cum_5FU6_D60_K4_tmp)[,c(1,2,3,6)],
                 chrom, chromstart,chromend,
                 range = c(0,max3),
                 addscale = TRUE,
                 ylab="sc_5FU_K4",
                 color="#00544CFF",cex.lab=1,cex.main=2.1)
    
    par(mar = c(1, 6, 1, 1),xpd=NA)
    
    plotGenes(gl.sub, chrom,chromstart,chromend,
              bentline=F,plotgenetype = "arrow",
              labeltext = T,labelat = "start",fontsize=1,labeloffset = 0.4,bheight=0.07)
    geneinfo=gl.sub
    
    par(mar = c(1, 6, 1, 1))
    labelgenome(chrom,chromstart,chromend,n=4,scale="Mb",cex.axis=1.5)
    
    dev.off()
}

