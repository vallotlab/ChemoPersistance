library(here)

maindir= here()
source(file.path(maindir, "Scripts", "global_var.R"))
source(file.path(maindir, "Scripts", "functions.R"))

inputDir <- file.path(maindir,"input","bulk_ChIPseq","MM468", "ChIPreChIP")
outputDir <- file.path(maindir,"output","bulk_ChIPseq","MM468", "ChIPreChIP")
resDir <- file.path(maindir,"output","bulk_ChIPseq","Snapshots")
load(file.path(maindir,"annotation","Gencode_hg38_v25.RData"))

print(list.files(inputDir))

###Load data#########
## K27 -> K4
chiprechip_H3K27me3_H3K4me3_vs_primaryK27 <- rtracklayer::import(file.path(inputDir, "MM468_reChIP_H3K27me3_H3K4me3_vs_primary_H3K27me3.bw"))
chiprechip_H3K27me3_IgG_vs_primaryK27 <- rtracklayer::import(file.path(inputDir, "MM468_reChIP_H3K27me3_IgG_vs_primary_H3K27me3.bw"))

bed_chiprechip_H3K27me3_H3K4me3 <- rtracklayer::import(file.path(outputDir, "summits_K27_K4_peak_significant_0.001_0.15.bed"))

## K4 -> K27
chiprechip_H3K4me3_H3K27me3_vs_primaryK4 <- rtracklayer::import(file.path(inputDir, "MM468_reChIP_H3K4me3_H3K27me3_vs_primary_H3K4me3.bw"))
chiprechip_H3K4me3_IgG_vs_primaryK4 <- rtracklayer::import(file.path(inputDir, "MM468_reChIP_H3K4me3_IgG_vs_primary_H3K4me3.bw"))

bed_chiprechip_H3K4me3_H3K27me3 <- rtracklayer::import(file.path(outputDir, "summits_K4_K27_peak_significant_0.001_4.bed"))
gc()

###Load HUMAN annotation#########
subGenebed <-  GencodeByGene[GencodeByGene$Gene_biotype %in% c("protein_coding"),]
genebed <- data.frame(chrom=subGenebed$Chr,start=subGenebed$Start,stop=subGenebed$End,gene=subGenebed$Gene_name,score=".",strand=subGenebed$Strand)
genebed$chrom = as.character(genebed$chrom)
genebed$gene = as.character(genebed$gene)

K27_K4_persister_genes = readxl::read_xls(file.path(outputDir,"Persister_K27_bivalence.xls"))
K4_K27_persister_genes = readxl::read_xls(file.path(outputDir,"Persister_K4_K27_bivalence.xls"))

K27_K4_persister_genes = K27_K4_persister_genes$Symbol[which(as.logical(K27_K4_persister_genes$bivalent))]
K4_K27_persister_genes = K4_K27_persister_genes$Symbol[which(as.logical(K4_K27_persister_genes$bivalent))]

interesting_genes = union(K4_K27_persister_genes, K27_K4_persister_genes)

regions = list(
    c("CDKN2A","chr9",21982199,22009918),
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
    c("ABCC3","chr17",50620203,50661087),
    c("ABCA5","chr17",69296235,69346469),
    c("LBH","chr2",30220798,30257330),
    c("HLA-A","chr6",29940470,29947884),
    c("HLA-B","chr6",31351866,31359245),
    c("HLA-C","chr6",31266749,31274136),
    c("TAP1_TAP2-PSBM8-9","chr6",32820532,32862767),
    c("NLRC5","chr16",56987485,57085524),
    c("B2M","chr15",44709487,44720159),
    c("TAPBP","chr6",33297695,33316387),
    c("THAP1","chr8",42721854,42896999),
    c("ATF3","chr1",212546419,212646966),
    c("JUN","chr1",58778791,58786047),
    c("JUND","chr19",18277694,18283622),
    c("RUNX1","chr21", 34604667, 36034668)
)


regions = list()
interesting_genes = unique(c(interesting_genes, "KLF4","TP63","FOSL1","BNC1",sapply(regions, function(i) i[1])))
for(gene in interesting_genes){
    regions[[gene]] = c(gene,
                        genebed$chrom[which(genebed$gene == gene)],
                        genebed$start[which(genebed$gene == gene)] - 5000,
                        genebed$stop[which(genebed$gene == gene)] + 5000)
}


samples="ChIPreChIP_compare_K27_K4_K4_K27"
if(!dir.exists(file.path(resDir,samples))) dir.create(file.path(resDir,samples))

for(i in 41:length(regions)){
    #Select region of interest
    region=regions[[i]][1]
    chrom=regions[[i]][2]
    chromstart=as.numeric(regions[[i]][3])
    chromend=as.numeric(regions[[i]][4])
    
    roi = GRanges(chrom,ranges = IRanges(chromstart,chromend), Gene=region)
    
    cat("Doing plot for ",region,"",chrom,":",chromstart,"-",chromend,".\n")
    gl.sub <- genebed[ which(genebed[,"chrom"] == chrom ),]
    
    ## K27 -> K4
    chiprechip_H3K27me3_H3K4me3_vs_primaryK27_tmp = chiprechip_H3K27me3_H3K4me3_vs_primaryK27[subjectHits(findOverlaps(roi,chiprechip_H3K27me3_H3K4me3_vs_primaryK27)),]
    chiprechip_H3K27me3_IgG_vs_primaryK27_tmp = chiprechip_H3K27me3_IgG_vs_primaryK27[subjectHits(findOverlaps(roi,chiprechip_H3K27me3_IgG_vs_primaryK27)),]
    max = round(max(c(chiprechip_H3K27me3_H3K4me3_vs_primaryK27_tmp$score,chiprechip_H3K27me3_IgG_vs_primaryK27_tmp$score)),3)
    min = round(min(c(chiprechip_H3K27me3_H3K4me3_vs_primaryK27_tmp$score,chiprechip_H3K27me3_IgG_vs_primaryK27_tmp$score)),3)
    
    ## K4 -> K27
    chiprechip_H3K4me3_H3K27me3_vs_primaryK4_tmp = chiprechip_H3K4me3_H3K27me3_vs_primaryK4[subjectHits(findOverlaps(roi,chiprechip_H3K4me3_H3K27me3_vs_primaryK4)),]
    chiprechip_H3K4me3_IgG_vs_primaryK4_tmp = chiprechip_H3K4me3_IgG_vs_primaryK4[subjectHits(findOverlaps(roi,chiprechip_H3K4me3_IgG_vs_primaryK4)),]
    max2 = round(max(c(chiprechip_H3K4me3_H3K27me3_vs_primaryK4_tmp$score,chiprechip_H3K4me3_IgG_vs_primaryK4_tmp$score)),3)
    min2 = round(min(c(chiprechip_H3K4me3_H3K27me3_vs_primaryK4_tmp$score,chiprechip_H3K4me3_IgG_vs_primaryK4_tmp$score)),3)
    
    if(is.infinite(max)) max = 1
    if(is.infinite(min)) min = 1
    if(is.infinite(max2)) max2 = 1
    if(is.infinite(min2)) min2 = 1

    if(length(chiprechip_H3K27me3_IgG_vs_primaryK27_tmp) !=0 & length(chiprechip_H3K27me3_H3K4me3_vs_primaryK27_tmp) !=0 |
       length(chiprechip_H3K4me3_H3K27me3_vs_primaryK4_tmp) !=0 & length(chiprechip_H3K4me3_IgG_vs_primaryK4_tmp) !=0 ){
        pdf(file.path(file.path(resDir,samples),paste0(samples,"_",region,".pdf")),width=8,height=6)
        layout.matrix <- matrix(c(1,1,1,2,2,2,3,3,4,4,4,5,5,5,6,6,7,7,8), ncol = 1)
        
        layout(mat = layout.matrix,
               heights = c(1), # Heights of the two rows
               widths = c(1)) # Widths of the two columns
        
        par(cex=0.5)
        par(mar = c(0.75, 6, 0, 1), oma = c(1, 1, 1, 1))
        
        plotBedgraph(as.data.frame(chiprechip_H3K27me3_H3K4me3_vs_primaryK27_tmp)[,c(1,2,3,6)],
                     chrom, chromstart,chromend,
                     range = c(min, max),
                     addscale = TRUE,
                     ylab="H3K27me3_H3K4me3 vs H3K27me3",
                     color="#afafafff",cex.lab=1,cex.main=2.1)
        
        
        plotBedgraph(as.data.frame(chiprechip_H3K27me3_IgG_vs_primaryK27_tmp)[,c(1,2,3,6)],
                     chrom, chromstart,chromend,
                     range = c(min,max),
                     addscale = TRUE,
                     ylab="H3K27me3_IgG vs H3K27me3",
                     color="#afafafff",cex.lab=1,cex.main=2.1)

        plotGenes(as.data.frame(bed_chiprechip_H3K27me3_H3K4me3), chrom,chromstart,chromend,
                  bentline=FALSE,plotgenetype = "box", col="dark red",
                  labeltext = F)
        
        plotBedgraph(as.data.frame(chiprechip_H3K4me3_H3K27me3_vs_primaryK4_tmp)[,c(1,2,3,6)],
                     chrom, chromstart,chromend,
                     range = c(min2, max2),
                     addscale = TRUE,
                     ylab="H3K4me3_H3K27me3 vs H3K4me3",
                     color="#afafafff",cex.lab=1,cex.main=2.1)
        
        
        plotBedgraph(as.data.frame(chiprechip_H3K4me3_IgG_vs_primaryK4_tmp)[,c(1,2,3,6)],
                     chrom, chromstart,chromend,
                     range = c(min2,max2),
                     addscale = TRUE,
                     ylab="H3K4me3_IgG vs H3K4me3",
                     color="#afafafff",cex.lab=1,cex.main=2.1)
        
        plotGenes(as.data.frame(bed_chiprechip_H3K4me3_H3K27me3), chrom,chromstart,chromend,
                  bentline=F,plotgenetype = "box", col="dark red",
                  labeltext = F)
        
        par(mar = c(1, 6, 1, 1),xpd=NA)
        
        
        
        plotGenes(gl.sub, chrom,chromstart,chromend,
                  bentline=F,plotgenetype = "arrow",
                  labeltext = T,labelat = "start",fontsize=1,labeloffset = 0.4,bheight=0.07)
        geneinfo=gl.sub
        
        par(mar = c(1, 6, 1, 1))
        labelgenome(chrom,chromstart,chromend,n=4,scale="Mb",cex.axis=1.5)
        
        dev.off()
    }
    
}

