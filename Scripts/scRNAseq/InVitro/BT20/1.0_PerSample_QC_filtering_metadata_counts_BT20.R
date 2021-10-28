library(here)
source(file.path(here(),"Scripts","functions.R"))
source(file.path(here(),"Scripts","global_var.R"))

maindir = here()

# Directories -------------------------------------------------------------
outdir = file.path(maindir, "output","scRNAseq", "BT20");  if(!dir.exists(outdir)){dir.create(outdir)}
QCdir <- file.path(outdir, "QC"); if(!dir.exists(QCdir)){dir.create(QCdir)}
resdir <- file.path(outdir, "Persister"); if(!dir.exists(resdir)){dir.create(resdir)}
resSUBdir <- file.path(resdir, paste0("Unsupervised")) ; if(!file.exists(resSUBdir)){dir.create(resSUBdir)}
RDatadir <- file.path(resSUBdir,"RData") ; if(!file.exists(RDatadir)){dir.create(RDatadir)}
RDatadirSamples <- file.path(RDatadir,"RData_perSample") ; if(!file.exists(RDatadirSamples)){dir.create(RDatadirSamples)}

# Replace GEO path with the path to the GEO directory downloaded on your computer (scRNA BT20)
GEOdir_scRNA_BT20 = file.path(maindir,"input","scRNAseq","MM468") 

# Load Data ---------------------------------------------------------------
BCmodel <- "BT20"

scRNA = list()
for(file in list.files(GEOdir_scRNA_BT20,pattern = "BT20.*_barcodes.tsv.gz", full.names = TRUE )){
    prefix = gsub("barcodes.tsv.gz", "", file)
    name = gsub("_barcodes.tsv.gz", "", basename(file))
    scRNA[[name]] = DropletUtils::read10xCounts(prefix, type = "prefix")
    gc()
}

# QC, normalization  ------------------------------------------------------
for(i in seq_along(scRNA)){
    umi = scRNA[[i]]
    name = names(scRNA)[i]
    keep_feat <- rowSums(counts(umi)>0)>0 #remove genes that are not expressed in any cells
    umi <- umi[keep_feat,]
    
    mt = which(grepl("^MT-", rowData(umi)[[2]]))
    ercc = which(grepl("ERCC", rowData(umi)[[2]]))
    control_features = rowData(umi)$ID[c(mt,ercc)]
    cellpool = seq_along(rowData(umi)[[2]])[-c(mt,ercc)]
    
    subsets = list("mt" = mt, "ercc" = ercc)
    rb.genes <- rowData(umi)[[1]][grep("^RP[SL]",rowData(umi)[[2]])]
    
    
    umi <- scater::addPerCellQC(umi, subsets = subsets)
    umi <- scater::addPerFeatureQC(umi)
    
    ##Manual quality filters:
    #Distribution of number of total detected features should be normal, get rid of the heavy right tail (<6000)
    umi$filter_by_expr_features <- (umi$detected<8000 & umi$detected>3000)
    #Distribution of number of total number of counts
    umi$filter_by_lib_size <- (umi$sum < 100000)
    #MT contamination
    umi$filter_by_MT_conta <- (umi$subsets_mt_percent<15)
    #Spike-in
    umi$filter_by_spikein <- (umi$subsets_ercc_percent<0.05)
    
    umi$use <- (umi$filter_by_lib_size & umi$filter_by_MT_conta & umi$filter_by_spikein & umi$filter_by_expr_features)
    
    umi.qc <- umi[,umi$use]
    
    #remove MT & ERCC genes from dataset for further analysis
    umi.qc <- umi.qc[-which(rowData(umi.qc)$ID %in% control_features),]
    dim(umi.qc)
    
    #Normalization by library size only
    sizeFactors(umi.qc) <- librarySizeFactors(umi.qc)
    umi.qc <- scater::logNormCounts(umi.qc) #normalized log2 expression values
    umi.qc <- scater::logNormCounts(umi.qc,log = FALSE) #normalized expression values
    
    
    pdf(file.path(QCdir,paste0(name,"_QC_scRNAseq.pdf")))
    scater::plotColData(umi,x="detected",y="subsets_mt_percent")
    scater::plotColData(umi,x="detected",y="subsets_ercc_percent")
    hist(umi$detected,breaks=100,xlim=c(0,10000),main=" Number of detected genes",xlab="NUMBER of GENES")
    abline(v=6000,col="red")
    abline(v=1500,col="red")
    hist(umi$sum,breaks=100, xlim=c(0,300000),main="Library size",xlab="TOTAL COUNTS")
    M <- median(umi$total_counts)
    abline(v=100000,col="red")
    dev.off()
    
    unc <- counts(umi.qc)
    rownames(unc) <- rowData(umi.qc)$Symbol
    
    
    set.seed(100)
    # test <- scran::doubletCells(umi.qc)
    test <- scDblFinder::computeDoubletDensity(umi.qc)
    metadata <- data.frame(barcode= substr(colData(umi.qc)$Barcode,0,16),
                           cell_id= paste0(name, "_", substr(colData(umi.qc)$Barcode,0,16)),
                           sample_id=rep(name,dim(umi.qc)[2]),
                           total_features=colData(umi.qc)$detected,
                           total_counts=colData(umi.qc)$sum,
                           doublet_score=test)
    gc()
    
    colnames(unc) <- metadata$cell_id
    assign(paste0("counts.",name),unc)
    
    assign(paste0("metadata.",name),metadata)
    save(list=c(paste0("counts.",name),paste0("metadata.",name)),
         file=file.path(RDatadirSamples,paste0(name,".RData")))
    rm(list=(paste0("counts.",name)))
}



