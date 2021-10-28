library(here)
source(file.path(here(),"Scripts","functions.R"))
source(file.path(here(),"Scripts","global_var.R"))

# Directories -------------------------------------------------------------
maindir = here()
resdir <- file.path(maindir,"output","scChIPseq","MM468","scH3K27me3"); if(!file.exists(resdir)){dir.create(resdir)}
unsupervised_dir  <- file.path(resdir,"Unsupervised","RData")
input_dir  <- file.path(maindir, "input","scChIPseq","MM468")

# Unsupervised directories
resdir <- file.path(resdir, "Supervised") ; if(!file.exists(resdir)){dir.create(resdir)}
resdir_tables = file.path(resdir,"Tables"); if(!dir.exists(resdir_tables)) dir.create(resdir_tables)
resdir_enrich <- file.path(resdir_tables,"Gene_set_analysis");if(!file.exists(resdir_enrich)){dir.create(resdir_enrich)}

dataset_name = "MM468_H3K27me3_peaks"
ref_genome ="hg38"

metadata = data.frame(sample_id = c(
    "MM468_DMSO1_day60", "MM468_DMSO3_day77", "MM468_DMSO5_day131",
    "MM468_5FU1_day33", "MM468_5FU2_day67",
    "MM468_5FU6_day131", "MM468_5FU3_day147", "MM468_5FU2_day171"),
    sample_id_color = c("#dfdfdfff", "#999999ff","#363636",
                        "#118675ff", "#8cc453ff",
                        "#ff5722ff", "#feb40fff", "#fd8508ff"))

# Run ChromSCape :
set.seed(47)
qval.th = 0.1
logFC.th =1
enrichment_qval = 0.1

# Load correlation filtered scExp
scExp_cf = qs::qread(file = file.path(unsupervised_dir, "MM468_H3K27me3_peaks_3000_10000_95_uncorrected_correlation_filtered.qs"))

# Load gene TSS annotation and create peak to gene association
annotK27 <- rtracklayer::import(file.path(
    maindir,"annotation","MM468_peaks_K27.bed.gz"))
start(annotK27)  = start(annotK27) -1
annotK27$ID = paste0(seqnames(annotK27),"_",start(annotK27),"_",end(annotK27))
annot_10k = read.table(gzfile(file.path(maindir,"annotation",
                                        "gencode.v34.annotation.transcriptTSS_10k.bed.gz")),
                       sep="\t", header = F)
colnames(annot_10k) = c("chr","start","end","transcripts","gene","strand")
annot_10k = as(annot_10k,"GRanges")
annot_10k$ID = paste0(seqnames(annot_10k),"_",start(annot_10k),"_", end(annot_10k))
start(annot_10k) = start(annot_10k) + 4500
end(annot_10k) = start(annot_10k)  + 1000

annotK27$gene = ""
hits <- findOverlaps(annotK27, annot_10k)
agg <- aggregate(annot_10k, hits, gene=paste(gene, collapse = ","))
annotK27$gene[match(subsetByOverlaps(annotK27,annot_10k)$name,annotK27$ID)] = agg$gene

###############################################################################
# Run differential analysis Persister (C3) vs Untreated (C2 & C4) 
# Using single-cell as well as bulk ChIP-seq
# On consensus peak annotation
################################################################################

scExp_Pers_vs_DMSO = scExp_cf
scExp_Pers_vs_DMSO = scExp_Pers_vs_DMSO[,which(scExp_Pers_vs_DMSO$cell_cluster %in% c("C1","C2","C4"))]
scExp_Pers_vs_DMSO$cell_cluster[which(scExp_Pers_vs_DMSO$cell_cluster %in% c("C1"))] = "bis"
scExp_Pers_vs_DMSO$cell_cluster[which(scExp_Pers_vs_DMSO$cell_cluster %in% c("C2","C4"))] = "C1"
scExp_Pers_vs_DMSO$cell_cluster[which(scExp_Pers_vs_DMSO$cell_cluster %in% c("bis"))] = "C2"

# Grouped DA
mat = NULL
for(cluster in c("C1","C2")){
    samps = unique(scExp_Pers_vs_DMSO$sample_id[which(scExp_Pers_vs_DMSO$cell_cluster == cluster)])
    tot = table(scExp_Pers_vs_DMSO$sample_id[which(scExp_Pers_vs_DMSO$cell_cluster == cluster)])
    for(sample in names(tot)){
        if(tot[sample]>50){
            col = as.matrix(rowSums(counts(scExp_Pers_vs_DMSO)[,which(scExp_Pers_vs_DMSO$cell_cluster == cluster &
                                                                          scExp_Pers_vs_DMSO$sample_id == sample)]))
            colnames(col) = paste0(sample,"_",cluster)
            if(is.null(mat)) mat = col else mat = cbind(mat,col)
        } 
    } 
}

# Read bulk Count matrices
resultZ <- read.csv(gzfile(file.path(maindir,"input","bulk_ChIPseq","InVitro","bulk_ChIPseq","Count_Matrices","CountTable_bulk_MM468_H3K27me3_peaks.csv.gz")))
resultZ$ID  <- paste(resultZ$Chromosome,resultZ$Begin,resultZ$End,sep="_")
resultZ = resultZ[,-c(1,2,3)]
bulk_RZ = resultZ
colnames(bulk_RZ) [grep("DMSO",colnames(bulk_RZ))] =
    paste0(colnames(bulk_RZ) [grep("DMSO",colnames(bulk_RZ))],"_C1")
colnames(bulk_RZ) [grep("5FU",colnames(bulk_RZ))] =
    paste0(colnames(bulk_RZ) [grep("5FU",colnames(bulk_RZ))],"_C2")

bulk_RZ = bulk_RZ[match(rownames(mat),bulk_RZ$ID),]
rownames(bulk_RZ) = bulk_RZ$ID
bulk_RZ = bulk_RZ[,-ncol(bulk_RZ)]

mat = cbind(mat, bulk_RZ)
myrefs <- list(
    DMSO = colnames(mat)[grep("C1",colnames(mat))]
)
mygps <- list(
    Persister = colnames(mat)[grep("C2",colnames(mat))]
)

refs <- names(myrefs)
groups <- names(mygps)

#selection of peaks with at least a log2 RPKM of 1 in one sample
RZ_sel = mat
feature = as.data.frame(scExp_Pers_vs_DMSO@rowRanges)
feature = feature[,c(1,2,3)]
annot = data.frame(sample_id = colnames(mat), cluster = gsub(".*_","",colnames(mat)))
res <- geco.ChIPseqCompareLimma(mat=RZ_sel,
                                metadata =annot,
                                ref=myrefs,
                                groups=mygps,
                                featureTab=feature
)

res$Gene = annotK27$gene[match(res$id,gsub(":|-","_",annotK27$ID))]
under_res = res %>% dplyr::filter(log2FC.Persister < -logFC.th & qval.Persister < qval.th) %>% 
    dplyr::arrange(qval.Persister) %>% dplyr::select(id,seqnames,start,end,
                                                     log2FC.Persister,qval.Persister, Gene)
over_res = res %>% dplyr::filter(log2FC.Persister > logFC.th & qval.Persister < qval.th) %>% 
    dplyr::arrange(qval.Persister) %>% dplyr::select(id,seqnames,start,end,
                                                     log2FC.Persister,qval.Persister, Gene)

WriteXLS(c("over_res","under_res","res"), 
         ExcelFileName = file.path(resdir_tables, "DA_grouped_Persister_vs_DSMO_with_bulk_peaks.xlsx"),
         SheetNames = c(paste0("Over_",round(logFC.th,2),"_",qval.th,"_n",nrow(over_res)),
                        paste0("Under_-",round(logFC.th,2),"_",qval.th,"_n",nrow(under_res)),"All"),
         perl = "perl", verbose = FALSE, row.names = FALSE, col.names = TRUE, AdjWidth = T,
         AutoFilter = TRUE, BoldHeaderRow = TRUE, na = "", FreezeRow = 1, FreezeCol = 1)

summaryTab <- geco.summaryCompareedgeR(restab=res,
                                       ref=myrefs,
                                       groups=mygps,
                                       qval.th=qval.th,
                                       fc.th=2,
                                       plotdir=file.path(resdir_tables))

# Anticorrelation with scRNA
scRNA = readxl::read_xlsx(file.path(maindir, "output","scRNAseq","MM468","Persister",
                                    "Supervised","Tables","Differential_analysis_Limma_logFC_1.58.xlsx"),sheet= 3)
res_unique_gene = res %>% tidyr::separate_rows(Gene) %>% unique %>% dplyr::select(id, log2FC.Persister,
                                                                                  qval.Persister,Gene)
scRNA_scChIP = left_join(scRNA, res_unique_gene,by =c("Symbol"="Gene"))
top_demethylated <- head(scRNA_scChIP$Symbol[order(scRNA_scChIP$log2FC.Persister,decreasing=F)],n=15)

pdf(file.path(resdir_tables,"scRNA_scChIP_log2FC_Persister_sc_grouped_vs_DMSO_peaks.pdf"))
sp <- ggscatter(scRNA_scChIP, x = "log2FC.Persister", y = "log2FC.C2_pers",
                add = "reg.line", ylab = "scRNA log2FC", xlab="scChIP K27 log2FC",  # Add regressin line
                add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE,
                label="Symbol",repel=TRUE,
                label.select= top_demethylated) #add size of dot function of initial expression in DMSO
# Add correlation coefficient
print(sp + stat_cor(method = "pearson",label.x=1,label.sep = "\n"))

# With only diff chip seq
scRNA_scChIP_filt = scRNA_scChIP %>% dplyr::filter(abs(log2FC.Persister) > logFC.th,
                                                   qval.Persister < qval.th)
sp <- ggscatter(scRNA_scChIP_filt, x = "log2FC.Persister", y = "log2FC.C2_pers",
                add = "reg.line", ylab = "scRNA log2FC", xlab="scChIP K27 log2FC",  # Add regressin line
                add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE,
                label="Symbol",repel=TRUE,
                label.select= top_demethylated) #add size of dot function of initial expression in DMSO
# Add correlation coefficient
print(sp + stat_cor(method = "pearson",label.x=1,label.sep = "\n"))
dev.off()

# Gene Set Analysis GROUPED
run_GSA = TRUE
if(run_GSA){
    annotbase <- "MSigDB" 
    database <- MSIG.ls ##MSigDB
    Overexpressed  <- Underexpressed <- data.frame()
    data("hg38.GeneTSS")
    GencodeGene = hg38.GeneTSS 
    possibleGenes <- union( unique(unlist(strsplit(res$Gene,split=","))),
                            c(as.character(unique(GencodeGene$gene)),"ELFN2"))
    gp <- groups
    ref <- refs
    print(paste0("Processing ",gp, " vs ", ref, " _ ",annotbase ))
    
    signific <- which(res[,paste("qval",gp,sep=".")] <= qval.th & abs(res[,paste("log2FC",gp,sep=".")]) > logFC.th)
    over <- which(res[,paste("qval",gp,sep=".")] <= qval.th & res[,paste("log2FC",gp,sep=".")] > logFC.th)
    under <- which(res[,paste("qval",gp,sep=".")] <= qval.th & res[,paste("log2FC",gp,sep=".")] < -logFC.th)
    print(paste0("significant = ", length(signific))) ; print(paste0("over = ", length(over))) ; print(paste0("under = ", length(under)))
    
    if(length(over)){
        enrich.test <- geco.enrichmentTest(gene.sets=database,
                                           mylist = unique(unlist(strsplit(res$Gene[over],split=","))),
                                           possibleIds=possibleGenes)
        enrich.test <- data.frame(Gene_set_name=rownames(enrich.test), enrich.test, check.names=FALSE)
        enrich.test <- merge( subset(MSIG.gs, select=-Genes), enrich.test, by.x="Gene.Set", by.y="Gene_set_name", all.y=TRUE, sort=FALSE ) ## Get class of gene set
        enrich.test <- enrich.test[order(enrich.test$`p-value`),]
        Overexpressed  <- enrich.test
    }
    if(length(under)){
        enrich.test <- geco.enrichmentTest(gene.sets=database,
                                           mylist= unique(unlist(strsplit(res$Gene[under],split=","))),
                                           possibleIds=possibleGenes)
        enrich.test <- data.frame(Gene_set_name=rownames(enrich.test), enrich.test, check.names=FALSE)
        enrich.test <- merge( subset(MSIG.gs, select=-Genes), enrich.test, by.x="Gene.Set", by.y="Gene_set_name", all.y=TRUE, sort=FALSE ) ## Get class of gene set
        enrich.test <- enrich.test[order(enrich.test$`p-value`),]
        ind <- which(enrich.test$`q-value`<= 0.1);if(!length(ind)){ind <- 1:20}
        Underexpressed <- enrich.test[ind,]		
    }
    
    WriteXLS(
        c("Overexpressed", "Underexpressed"),
        ExcelFileName = file.path(resdir_enrich, 
                                  paste0("Enrichment_test_",gp,"_vs_",ref,
                                         "_logFC",round(logFC.th,2),"_peaks.xlsx")), 
        SheetNames = c( paste0("Overexp_in_", gp), paste0("Underexp_in_", gp) ),
        perl = "perl", verbose = FALSE, row.names = FALSE, col.names = TRUE,
        AdjWidth = TRUE, AutoFilter = TRUE, BoldHeaderRow = TRUE, na = "",
        FreezeRow = 1, FreezeCol = 1)
}


################################################################################
# Run differential analysis Persister (C3) vs Untreated (C2 & C4) 
# Using single-cell as well as bulk ChIP-seq
# On gene TSS annotation
################################################################################

out <- import_scExp(grep("K4", invert = TRUE, list.files(file.path(input_dir,"Count_Matrices"),
                                                         pattern = "TSS", full.names = T), value = TRUE),
                       remove_pattern = "_H3K27me3_TSS")

# Save raw
datamatrix_10k = out$datamatrix
annot_raw_10k = out$annot_raw
scExp_10k = create_scExp(datamatrix_10k, annot_raw_10k)

# Excluding known CNA loci
exclude_regions = rtracklayer::import(file.path(maindir,"annotation","MM468_identified_CNA.bed.gz"))
scExp_10k = exclude_features_scExp(scExp_10k, exclude_regions)

datamatrix_10k = counts(scExp_10k)
mat = NULL
for(cluster in c("C1","C2")){
    samps = unique(scExp_Pers_vs_DMSO$sample_id[which(scExp_Pers_vs_DMSO$cell_cluster == cluster)])
    tot = table(scExp_Pers_vs_DMSO$sample_id[which(scExp_Pers_vs_DMSO$cell_cluster == cluster)])
    for(sample in names(tot)){
        cells = scExp_Pers_vs_DMSO$cell_id[which(scExp_Pers_vs_DMSO$cell_cluster == cluster &
                                                     scExp_Pers_vs_DMSO$sample_id == sample)]
        cells = intersect(cells,colnames(datamatrix_10k))
        
        if(length(cells)>50){
            col = as.matrix(rowSums(datamatrix_10k[,match(cells, colnames(datamatrix_10k))]))
            print(paste0(sample, " - ", length(cells)))
            colnames(col) = paste0(sample,"_",cluster)
            if(is.null(mat)) mat = col else mat = cbind(mat,col)
        } 
    } 
}

# Read bulk Count matrices
resultZ <- read.csv(gzfile(file.path(maindir,"input","bulk_ChIPseq","InVitro","bulk_ChIPseq","Count_Matrices","CountTable_bulk_MM468_H3K27me3_TSS.csv.gz")))
row = paste0(resultZ$Chromosome, ":", resultZ$Begin, "-",resultZ$End)
dups = duplicated(row)
resultZ = resultZ[!dups,]
rownames(resultZ) <- row[!dups]
feature = as.data.frame(resultZ[,c(1,2,3)])
rownames(feature) = gsub(":|-","_",rownames(feature))
resultZ = resultZ[,-c(1,2,3)]

bulk_RZ = resultZ
rownames(bulk_RZ) = gsub(":|-","_",rownames(bulk_RZ))
colnames(bulk_RZ) [grep("DMSO",colnames(bulk_RZ))] =
    paste0(colnames(bulk_RZ) [grep("DMSO",colnames(bulk_RZ))],"_bulk_C1")
colnames(bulk_RZ) [grep("5FU",colnames(bulk_RZ))] =
    paste0(colnames(bulk_RZ) [grep("5FU",colnames(bulk_RZ))],"_bulk_C2")

bulk_RZ = bulk_RZ[match(rownames(mat),rownames(bulk_RZ)),]
mat = cbind(mat, bulk_RZ)

myrefs <- list(
    DMSO = colnames(mat)[grep("C1",colnames(mat))]
)
mygps <- list(
    Persister = colnames(mat)[grep("C2",colnames(mat))]
)

refs <- names(myrefs)
groups <- names(mygps)

RZ_sel = mat
annot = data.frame(sample_id = colnames(mat), cluster = gsub(".*_","",colnames(mat)))
rownames(annot) = annot$sample_id
res <- geco.ChIPseqCompareLimma(mat = RZ_sel,
                                metadata = annot,
                                ref = myrefs,
                                groups = mygps,
                                featureTab = feature
)

res$Gene = annot_10k$gene[match(res$id,annot_10k$ID)]
under_res = res %>% dplyr::filter(log2FC.Persister < -logFC.th & qval.Persister < qval.th) %>% 
    dplyr::arrange(qval.Persister) %>% dplyr::select(id,Chromosome,Begin,End,
                                                     log2FC.Persister,qval.Persister, Gene)
over_res = res %>% dplyr::filter(log2FC.Persister > logFC.th & qval.Persister < qval.th) %>% 
    dplyr::arrange(qval.Persister) %>% dplyr::select(id,Chromosome,Begin,End,
                                                     log2FC.Persister,qval.Persister, Gene)

WriteXLS(c("over_res","under_res","res"), 
         ExcelFileName = file.path(resdir_tables,
                                   paste0("DA_grouped_Persister_vs_DSMO_with_bulk_TSS.xlsx")),
         SheetNames = c(paste0("Over_",round(logFC.th,2),"_",qval.th,"_n",nrow(over_res)),
                        paste0("Under_-",round(logFC.th,2),"_",qval.th,"_n",nrow(under_res)),"All"), perl = "perl", verbose = FALSE, row.names = FALSE, col.names = TRUE, AdjWidth = T, AutoFilter = TRUE, BoldHeaderRow = TRUE, na = "", FreezeRow = 1, FreezeCol = 1)

# Volcano plot
summaryTab <- geco.summaryCompareedgeR(restab=res,
                                       ref=myrefs,
                                       groups=mygps,
                                       qval.th=qval.th,
                                       fc.th=logFC.th,
                                       plotdir=file.path(resdir_tables))

# Anticorrelation with scRNA
scRNA = readxl::read_xlsx(file.path(maindir, "output","scRNAseq","MM468","Persister",
                                    "Supervised","Tables","Differential_analysis_Limma_logFC_1.58.xlsx"),sheet= 3)
res_unique_gene = res %>% group_by(Gene) %>%
    slice_max(abs(log2FC.Persister)) %>% dplyr::select(id, log2FC.Persister,qval.Persister,Gene)
scRNA_scChIP = left_join(scRNA, res_unique_gene,by =c("Symbol"="Gene"))

top_demethylated <- scRNA_scChIP %>% filter(log2FC.C2_pers > log2(3), scRNA_scChIP$log2FC.Persister < -1) %>% 
    arrange(desc(log2FC.Persister)) %>% dplyr::select(Symbol) %>% unique %>%  head(n=35) 
top_demethylated = top_demethylated$Symbol

options(ggrepel.max.overlaps = 25)
options(max.overlaps = 25)

pdf(file.path(resdir_tables,"scRNA_scChIP_log2FC_Persister_sc_grouped_vs_DMSO_TSS.pdf"))
sp <- ggscatter(scRNA_scChIP, x = "log2FC.Persister", y = "log2FC.C2_pers",
                add = "reg.line", ylab = "scRNA log2FC", xlab="scChIP K27 log2FC",  # Add regressin line
                add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE,
                label="Symbol",repel=TRUE,
                label.select= top_demethylated) #add size of dot function of initial expression in DMSO
# Add correlation coefficient
print(sp + stat_cor(method = "pearson",label.x=1,label.sep = "\n"))

# With only diff chip seq
scRNA_scChIP_filt = scRNA_scChIP %>% dplyr::filter(abs(log2FC.Persister) > logFC.th,
                                                   qval.Persister < qval.th)
sp <- ggscatter(scRNA_scChIP_filt, x = "log2FC.Persister", y = "log2FC.C2_pers",
                add = "reg.line", ylab = "scRNA log2FC", xlab="scChIP K27 log2FC",  # Add regressin line
                add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE,
                label="Symbol",repel=TRUE,
                label.select= top_demethylated, max.overlaps = 25) #add size of dot function of initial expression in DMSO
# Add correlation coefficient
print(sp + stat_cor(method = "pearson",label.x=1,label.sep = "\n"))
dev.off()

run_GSA = TRUE
if(run_GSA){
    annotbase <- "MSigDB" 
    database <- MSIG.ls ##MSigDB
    Overexpressed  <- Underexpressed <- data.frame()
    data("hg38.GeneTSS")
    GencodeGene = hg38.GeneTSS 
    # possibleGenes <- unique(unlist(strsplit(annotK27$gene,split=",")))
    possibleGenes <- unique(union(res$Gene,c(as.character(unique(GencodeGene$gene)),"ELFN2")))
    gp <- groups
    ref <- refs
    print(paste0("Processing ",gp, " vs ", refs, " _ ",annotbase ))
    
    signific <- which(res[,paste("qval",gp,sep=".")] <= qval.th & abs(res[,paste("log2FC",gp,sep=".")]) > logFC.th)
    over <- which(res[,paste("qval",gp,sep=".")] <= qval.th & res[,paste("log2FC",gp,sep=".")] > logFC.th)
    under <- which(res[,paste("qval",gp,sep=".")] <= qval.th & res[,paste("log2FC",gp,sep=".")] < -logFC.th)
    print(paste0("significant = ", length(signific))) ; print(paste0("over = ", length(over))) ; print(paste0("under = ", length(under)))
    
    if(length(over)){
        enrich.test <- geco.enrichmentTest(gene.sets=database,
                                           mylist = unique(unlist(strsplit(res$Gene[over],split=","))),
                                           possibleIds=possibleGenes)
        enrich.test <- data.frame(Gene_set_name=rownames(enrich.test), enrich.test, check.names=FALSE)
        enrich.test <- merge( subset(MSIG.gs, select=-Genes), enrich.test, by.x="Gene.Set", by.y="Gene_set_name", all.y=TRUE, sort=FALSE ) ## Get class of gene set
        enrich.test <- enrich.test[order(enrich.test$`p-value`),]
        Overexpressed  <- enrich.test
    }
    if(length(under)){
        enrich.test <- geco.enrichmentTest(gene.sets=database,
                                           mylist= unique(unlist(strsplit(res$Gene[under],split=","))),
                                           possibleIds=possibleGenes)
        enrich.test <- data.frame(Gene_set_name=rownames(enrich.test), enrich.test, check.names=FALSE)
        enrich.test <- merge( subset(MSIG.gs, select=-Genes), enrich.test, by.x="Gene.Set", by.y="Gene_set_name", all.y=TRUE, sort=FALSE ) ## Get class of gene set
        enrich.test <- enrich.test[order(enrich.test$`p-value`),]
        ind <- which(enrich.test$`q-value`<= 0.1);if(!length(ind)){ind <- 1:20}
        Underexpressed <- enrich.test[ind,]		
    }
    
    WriteXLS(
        c("Overexpressed", "Underexpressed"),
        ExcelFileName = file.path(resdir_tables,
                                  paste0("Enrichment_test_",gp,"_vs_",refs,
                                         "_logFC",round(logFC.th,2),"_TSS.xlsx")), 
        SheetNames = c( paste0("Overexp_in_", gp), paste0("Underexp_in_", gp) ),
        perl = "perl", verbose = FALSE, row.names = FALSE, col.names = TRUE,
        AdjWidth = TRUE, AutoFilter = TRUE, BoldHeaderRow = TRUE, na = "",
        FreezeRow = 1, FreezeCol = 1)
}

################################################################################
################ Enrichment in TSS #############################################
################################################################################

# Load bulk + sc grouped peak DA 
res = readxl::read_xlsx(file.path(resdir_tables,"DA_grouped_Persister_vs_DSMO_with_bulk_peaks.xlsx"), sheet = 3)
# Load annotation peak affectation 
annotK27_peak_affectation = read.table( file.path(maindir, "annotation", "annotK27_peak_affectation.tsv"), 
                                        sep = "\t", header = TRUE)
res$peak_affectation = annotK27_peak_affectation$peak_affectation[match(res$id,annotK27_peak_affectation$ID)]
col = data.frame("col"=c("#999999ff","#dededeff","#3b9ab2ff","#78b7c5ff","#ebcc2aff","#e1af00ff"))
sig_diff <- (abs(res$log2FC.Persister)> logFC.th & res$qval.Persister < qval.th)
sig_under <- (res$log2FC.Persister< -logFC.th & res$qval.Persister < qval.th)
sig_over <- (res$log2FC.Persister> logFC.th & res$qval.Persister < qval.th)
l = list("Diff" = sig_diff, "Under" = sig_under, "Over" = sig_over)

#Enrichment
fisher <- function(a,b,c,d){
    data <- matrix(c(a,b,c,d),ncol=2)
    c(p = fisher.test(data)$p.value)
}

sig = l[["Over"]]
sum_sig = sum(table(res$peak_affectation[sig]))
sum_all = sum(table(res$peak_affectation))
res$peak_affectation = factor(res$peak_affectation, levels =c(
    "intergenic", "intergenic_enhancer", "tss_pc", "genebody_pc", "tss_other", "genebody_other"
))
tab = res[sig,] %>% dplyr::count(peak_affectation,.drop=F) %>%
    left_join(res %>% dplyr::count(peak_affectation,.drop=F), by = "peak_affectation") %>%
    mutate(sum_sig =sum_sig, sum_all =sum_all, Freq = log2( (n.x/sum_sig) / (n.y/sum_all) ),
           color = col$col) %>% rowwise() %>% mutate(p=fisher(n.x,sum_sig-n.x,n.y,sum_all-n.y)) %>% 
    mutate(is_significative = ifelse(p <0.05, TRUE,FALSE),
           Freq = ifelse(is.infinite(Freq),0,Freq),
           fill = ifelse(is_significative,color,"white"))

png(file.path(resdir_tables,paste0("enrichment_HQ_enriched.png")),width = 1000,height = 1000,res=300)
p1 = tab %>% ggplot() + geom_bar(aes(x=peak_affectation, y=Freq),fill= tab$fill, color = tab$color,
                                 stat="identity") + ggtitle(paste0(i," enriched H3K27me3 peaks")) +
    theme_classic() + ggtitle(paste0("Enriched H3K27me3 peaks - ",length(which(sig))," peaks")) + 
    ylab("Log2 Enrichment") + xlab("")  + geom_hline(yintercept = 0, color ="grey")
print(p1 +  theme(axis.text.x = element_text(angle = 90, vjust =1), legend.position = "None"))
dev.off()

sig = l[["Under"]]
sum_sig = sum(table(res$peak_affectation[sig]))
sum_all = sum(table(res$peak_affectation))

tab = res[sig,] %>% dplyr::count(peak_affectation,.drop=F) %>%
    left_join(res %>% dplyr::count(peak_affectation,.drop=F), by = "peak_affectation") %>%
    mutate(sum_sig =sum_sig, sum_all =sum_all, Freq = log2( (n.x/sum_sig) / (n.y/sum_all) ),
           color = col$col) %>% rowwise() %>% mutate(p=fisher(n.x,sum_sig-n.x,n.y,sum_all-n.y)) %>% 
    mutate(is_significative = ifelse(p <0.05, TRUE,FALSE),
           Freq = ifelse(is.infinite(Freq),0,Freq),
           fill = ifelse(is_significative,color,"white"))

png(file.path(resdir_tables,paste0("enrichment_HQ_under.png")),width = 1000,height = 1000,res=300)
p2= tab %>% ggplot() + geom_bar(aes(x=peak_affectation, y=Freq),fill= tab$fill, color = tab$color,
                                stat="identity") + ggtitle(paste0(i," depleted H3K27me3 peaks")) +
    theme_classic() + ggtitle(paste0("Depleted H3K27me3 peaks - ",length(which(sig))," peaks")) +
    ylab("Log2 Enrichment") + xlab("")  + geom_hline(yintercept = 0, color ="grey")
print(p2 + theme(axis.text.x = element_text(angle = 90, vjust =1), legend.position = "None"))
dev.off()

png(file.path(resdir_tables,paste0("enrichment_HQ_over_wolabel.png")),width = 800,height = 1000,res=300)
print(p1 + theme(text=element_blank()))
dev.off()

png(file.path(resdir_tables,paste0("enrichment_HQ_under_wolabel.png")),width = 800,height = 1000,res=300)
print(p2 + theme(text=element_blank()))
dev.off()
