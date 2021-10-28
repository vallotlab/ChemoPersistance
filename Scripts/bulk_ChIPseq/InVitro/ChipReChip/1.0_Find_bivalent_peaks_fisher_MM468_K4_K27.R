library(here)

maindir= here()
resdir <- file.path(maindir,"output","bulk_ChIPseq","MM468","ChIPreChIP")
if(!dir.exists(resdir)) dir.create(resdir)
source(file.path(maindir,"Scripts","global_var.R"))


inputDir =  "/media/pprompsy/Depic_bioinfo_1/InstitutCurie/Documents/Data/results/ChIP_reChIP_D767/MM468/03subsample_bam/"# file.path(maindir, "input","bulk_ChIPseq","MM468","ChIPreChIP")
outputDir = file.path(maindir, "output","bulk_ChIPseq","MM468","ChIPreChIP")

# K27 K4 files
primary_K4 = file.path(inputDir,"MM468BC_primary_H3K4me3_H3K27me3_1.subsample.bam")
chipReChip_K4_K27 = file.path(inputDir,"MM468BC_reChIP_H3K4me3_H3K27me3_1.subsample.bam")
summits_K4_K27_primary = rtracklayer::import(file.path(maindir, "input","bulk_ChIPseq","InVitro","ChIPreChIP", "Peaks","H3K4me3_H3K27me3_primary_summits_merged_filtered_by_TSS.bed"))
summits_K4_K27_primary = unique(summits_K4_K27_primary)

# K27 IgG files
primary_K4_IgG = file.path(inputDir,"MM468BC_primary_H3K4me3_IgG_1.subsample.bam")
chipReChip_K4_IgG= file.path(inputDir,"MM468BC_reChIP_H3K4me3_IgG_1.subsample.bam")
summits_K4_IgG_primary = rtracklayer::import(file.path(maindir, "input","bulk_ChIPseq","InVitro","ChIPreChIP", "Peaks","H3K4me3_IgG_primary_summits_merged_filtered_by_TSS.bed"))
summits_K4_IgG_primary = unique(summits_K4_IgG_primary)

cat("Raw False positive rate is expected to be around ",round(length(summits_K4_IgG_primary)/length(summits_K4_K27_primary),3), ".\n")

# Summits K27 K4
summits_K4_K27_around = summits_K4_K27_primary
summits_K4_K27_close = summits_K4_K27_primary

# Close is +- 10000 around summit
start(summits_K4_K27_close) = start(summits_K4_K27_close) -2500
end(summits_K4_K27_close) = end(summits_K4_K27_close) + 2500
summary(width(summits_K4_K27_close))

# Get count for chip re chip K27 K4
filt <- Rsamtools:: ScanBamParam(which = summits_K4_K27_close)
primary_count_K4_K27_close = Rsamtools::countBam(primary_K4, param = filt)
colnames(primary_count_K4_K27_close) = paste0(colnames(primary_count_K4_K27_close),"_close")

# Get count for primary K27 (K4)
filt <- Rsamtools:: ScanBamParam(which = summits_K4_K27_close)
chiprechip_count_K4_K27_close = Rsamtools::countBam(chipReChip_K4_K27, param = filt)
colnames(chiprechip_count_K4_K27_close) = paste0(colnames(chiprechip_count_K4_K27_close),"_close")


fisher_pvalue_mat = matrix(0,nrow = nrow(chiprechip_count_K4_K27_close))
odd_ratio_mat = matrix(0,nrow = nrow(chiprechip_count_K4_K27_close))
signal_compared_to_close = matrix(0,nrow = nrow(chiprechip_count_K4_K27_close))
rownames(fisher_pvalue_mat) = paste0(chiprechip_count_K4_K27_close$space_close,":",chiprechip_count_K4_K27_close$start_close,"-",chiprechip_count_K4_K27_close$end_close)
rownames(odd_ratio_mat) = paste0(chiprechip_count_K4_K27_close$space_close,":",chiprechip_count_K4_K27_close$start_close,"-",chiprechip_count_K4_K27_close$end_close)

# Normalize by library size
primary_count_K4_K27_close$normalized_records = 10^6 * primary_count_K4_K27_close$records_close / sum(primary_count_K4_K27_close$records_close)
chiprechip_count_K4_K27_close$normalized_records = 10^6 * chiprechip_count_K4_K27_close$records_close / sum(chiprechip_count_K4_K27_close$records_close)

# Calculate the ratio ReChip vs Primary for each peaks
ratio = chiprechip_count_K4_K27_close$normalized_records / primary_count_K4_K27_close$normalized_records
names(ratio) = rownames(fisher_pvalue_mat)

png(file.path(outputDir,"Distribution_peak_ratio_K4_K27.png"),width = 1200, height = 1200, res = 200)
hist(ratio,breaks=150)
dev.off()

# ratio["chr6:1310085-1315544"]
# i = which(names(ratio)=="chr6:1310085-1315544")
# For each peak, retrieve the 10 peaks before & the 10 peaks after
# And test if the given peak ratio is significantly higher than the 20 peaks around
# 
# permutation.test <- function(treatment, outcome, original, n){
#     distribution = sapply(1:n, function(j){
#         samp = sample(treatment, length(treatment), FALSE)
#         print(length(outcome[(samp==1)]))
#         print(length(median(outcome[samp==0])))
#         print("#")
#         outcome[(samp==1)] - median(outcome[samp==0])
#     })
#     result=sum(abs(distribution) >= abs(original))/(n)
#     return(list(result, distribution))
# }

for(i in seq_len(nrow(fisher_pvalue_mat))){
    ratio_i = ratio[i]
    
    index_peaks_around = max(i-30,1):min(i+30,nrow(fisher_pvalue_mat))
    ratios_around = ratio[index_peaks_around]
    
    signal_peaks_around_primary = sum(primary_count_K4_K27_close$records_close[index_peaks_around])
    signal_peaks_around_reChIP = sum(chiprechip_count_K4_K27_close$records_close[index_peaks_around])
    
    # ratios_around = ratios_around[-min(i,16)]
    # ratio_i / median(ratios_around)
    
    # is_peak_of_interest = c(rep(0,min((i-1),30)),1,rep(0,min(nrow(fisher_pvalue_mat) - i,30)))
    # 
    # # Difference in means
    # original <- ratios_around[is_peak_of_interest==1] - median(ratios_around[(is_peak_of_interest==0)])
    # 
    # test1 <- permutation.test(is_peak_of_interest, ratios_around, original, 1000)
    # # hist(test1[[2]], breaks=50, col='grey', main="Permutation Distribution", las=1, xlab='')
    # abline(v=ratio_i, lwd=3, col="red")
    
    contingency_table = data.frame(a= c(chiprechip_count_K4_K27_close$records_close[i],
                                        signal_peaks_around_reChIP),
                                   b = c(primary_count_K4_K27_close$records_close[i],
                                         signal_peaks_around_primary))
    
    fisher_pvalue_mat[i,] = fisher.test(contingency_table,alternative = "greater")$p.value
    
    # fisher_pvalue_mat[i,] = test1[[1]]
    odd_ratio_mat[i,] = ratio_i / median(ratios_around)
    
    if(i%%100==0){print(i)}
}

fisher_pvalue_df = as.data.frame(as.matrix(fisher_pvalue_mat))
length(which(fisher_pvalue_df$V1 < 0.25))
length(which(fisher_pvalue_df$V1 < 0.15))
length(which(fisher_pvalue_df$V1 < 0.1))
length(which(fisher_pvalue_df$V1 < 0.05))
length(which(fisher_pvalue_df$V1 < 0.01))
length(which(fisher_pvalue_df$V1 < 0.005))
length(which(fisher_pvalue_df$V1 < 0.001))

fisher_pvalue_df$odd_ratio = as.numeric(odd_ratio_mat)
summary(fisher_pvalue_df$odd_ratio)

annot_10k = read.table(gzfile(file.path(maindir,"annotation",
                                        "gencode.v34.annotation.transcriptTSS_10k.bed.gz")),
                       sep="\t", header = F)
colnames(annot_10k) = c("chr","start","end","transcripts","gene","strand")
annot_10k = as(annot_10k,"GRanges")

# Intersect genes K4 K27 & persister genes
scRNA = readxl::read_xls(file.path(maindir,"output","scRNAseq","MM468","Persister","ComparisonChIPseq","K27","persister_K27_status.xls"))
depleted_K27 = scRNA$Symbol[which(!scRNA$K27_status %in% c("No K27 Peak","Not differential"))]

# Getting Pvalues for IgG
# IgG
# Summits K27 K4
summits_K4_IgG_close = summits_K4_IgG_primary

# Close is +- 10000 around summit
start(summits_K4_IgG_close) = start(summits_K4_IgG_close) - 2500
end(summits_K4_IgG_close) = end(summits_K4_IgG_close) + 2500

# Get count for chip re chip K27 K4
filt <- Rsamtools:: ScanBamParam(which = summits_K4_IgG_close)
primary_count_K4_IgG_close = Rsamtools::countBam(primary_K4_IgG, param = filt)
colnames(primary_count_K4_IgG_close) = paste0(colnames(primary_count_K4_IgG_close),"_close")

# Get count for primary K27 (K4)
filt <- Rsamtools:: ScanBamParam(which = summits_K4_IgG_close)
chiprechip_count_K4_IgG_close = Rsamtools::countBam(chipReChip_K4_IgG, param = filt)
colnames(chiprechip_count_K4_IgG_close) = paste0(colnames(chiprechip_count_K4_IgG_close),"_close")

fisher_pvalue_IgG_mat = matrix(0,nrow = nrow(chiprechip_count_K4_IgG_close))
odd_ratio_mat = matrix(0,nrow = nrow(chiprechip_count_K4_IgG_close))

rownames(fisher_pvalue_IgG_mat) = paste0(chiprechip_count_K4_IgG_close$space_close,":",chiprechip_count_K4_IgG_close$start_close,"-",chiprechip_count_K4_IgG_close$end_close)
rownames(odd_ratio_mat) = paste0(chiprechip_count_K4_IgG_close$space_close,":",chiprechip_count_K4_IgG_close$start_close,"-",chiprechip_count_K4_IgG_close$end_close)

# Normalize by library size
primary_count_K4_IgG_close$normalized_records = 10^6 * primary_count_K4_IgG_close$records_close / sum(primary_count_K4_IgG_close$records_close)
chiprechip_count_K4_IgG_close$normalized_records = 10^6 * chiprechip_count_K4_IgG_close$records_close / sum(chiprechip_count_K4_IgG_close$records_close)

# Calculate the ratio ReChip vs Primary for each peaks
ratio = chiprechip_count_K4_IgG_close$normalized_records / primary_count_K4_IgG_close$normalized_records
names(ratio) = rownames(fisher_pvalue_IgG_mat)

png(file.path(outputDir,"Distribution_peak_ratio_K4_IgG.png"),width = 1200, height = 1200, res = 200)
hist(ratio,breaks=150)
dev.off()

for(i in  seq_len(nrow(fisher_pvalue_IgG_mat))){
    ratio_i = ratio[i]
    index_peaks_around = max(i-30,1):min(i+30,nrow(fisher_pvalue_IgG_mat))
    ratios_around = ratio[index_peaks_around]
    signal_peaks_around_primary = sum(primary_count_K4_IgG_close$records_close[index_peaks_around])
    signal_peaks_around_reChIP = sum(chiprechip_count_K4_IgG_close$records_close[index_peaks_around])
    
    contingency_table = data.frame(a= c(chiprechip_count_K4_IgG_close$records_close[i],
                                        signal_peaks_around_reChIP),
                                   b = c(primary_count_K4_IgG_close$records_close[i],
                                         signal_peaks_around_primary))
    
    fisher_pvalue_IgG_mat[i,] = fisher.test(contingency_table,alternative = "greater")$p.value
    odd_ratio_mat[i,] = ratio_i / median(ratios_around)
    
    if(i%%100==0){print(i)}
}

fisher_pvalue_IgG_df = as.data.frame(as.matrix(fisher_pvalue_IgG_mat))
fisher_pvalue_IgG_df$odd_ratio = as.numeric(odd_ratio_mat)

length(which(fisher_pvalue_IgG_df$V1 < 0.25))
length(which(fisher_pvalue_IgG_df$V1 < 0.15 ))
length(which(fisher_pvalue_IgG_df$V1 < 0.1))
length(which(fisher_pvalue_IgG_df$V1 < 0.05 ))
length(which(fisher_pvalue_IgG_df$V1 < 0.01))
length(which(fisher_pvalue_df$odd_ratio > 5))

summary(fisher_pvalue_IgG_df$odd_ratio)
summary(fisher_pvalue_df$odd_ratio)

colnames(fisher_pvalue_df)[1] = "fisher.pvalue"
colnames(fisher_pvalue_IgG_df)[1] = "fisher.pvalue"

fisher_pvalue_df$qvalue = p.adjust(fisher_pvalue_df$fisher.pvalue,method = "BH")
fisher_pvalue_IgG_df$qvalue = p.adjust(fisher_pvalue_IgG_df$fisher.pvalue,method = "BH")

tab=data.frame()
for(qval in c(0.25,0.2,0.15,0.1,0.05,0.01,0.005,0.001,0.0001)){
    for(peak_ratio in c(1,2,3,4,5)){
        # select signif peaks for K27 K4
        fisher_pvalue_df$Significant = (fisher_pvalue_df$fisher.pvalue<qval & fisher_pvalue_df$odd_ratio > peak_ratio)
        
        summits_K4_K27_peak_significant = summits_K4_K27_primary[fisher_pvalue_df$Significant]
        summits_K4_K27_peak_significant$gene = ""
        summits_K4_K27_peak_significant$ID = paste0(seqnames(summits_K4_K27_peak_significant),":",
                                                    start(summits_K4_K27_peak_significant),"-",
                                                    end(summits_K4_K27_peak_significant))
        
        hits <- findOverlaps(summits_K4_K27_peak_significant, annot_10k)
        agg <- aggregate(annot_10k, hits, gene=paste(gene, collapse = ","))
        
        summits_K4_K27_peak_significant$gene[match(subsetByOverlaps(summits_K4_K27_peak_significant,annot_10k)$ID,
                                                   summits_K4_K27_peak_significant$ID)] = agg$gene
        
        # select significant peaks for K27 IgG
        fisher_pvalue_IgG_df$Significant = (fisher_pvalue_IgG_df$fisher.pvalue<qval & fisher_pvalue_IgG_df$odd_ratio > peak_ratio)
        
        if(length(which(fisher_pvalue_IgG_df$Significant))>0){
            summits_K4_IgG_peak_significant = summits_K4_IgG_primary[fisher_pvalue_IgG_df$Significant]
            summits_K4_IgG_peak_significant$gene = ""
            summits_K4_IgG_peak_significant$ID = paste0(seqnames(summits_K4_IgG_peak_significant),":",
                                                        start(summits_K4_IgG_peak_significant),"-",
                                                        end(summits_K4_IgG_peak_significant))
            
            hits <- findOverlaps(summits_K4_IgG_peak_significant, annot_10k)
            agg <- aggregate(annot_10k, hits, gene=paste(gene, collapse = ","))
            
            summits_K4_IgG_peak_significant$gene[match(subsetByOverlaps(summits_K4_IgG_peak_significant,annot_10k)$ID,summits_K4_IgG_peak_significant$ID)] = agg$gene
        } else {
            summits_K4_IgG_peak_significant = NULL
        }
        
        bivalent_genes = unique(unlist(strsplit(summits_K4_K27_peak_significant$gene,split=",")))
        cat("qval = ", qval, "; Ratio = ",peak_ratio," - num biv peaks = ", length(which(fisher_pvalue_df$Significant)),
            "; num depleted K27 persister bivalent = ", length(intersect(bivalent_genes, depleted_K27)),
            "; num  persister bivalent = ", length(intersect(bivalent_genes,scRNA$Symbol)),
            "False positive = ",length(summits_K4_IgG_peak_significant),".\n")
        # rtracklayer::export(summits_K4_K27_peak_significant,file.path(outputDir,paste0("summits_K4_K27_peak_significant_",qval,"_",ratio,".bed")))
        tab = rbind(tab, data.frame(
            "p.value" = qval, "ratio" = peak_ratio, "num_bivalent_peaks" = length(which(fisher_pvalue_df$Significant)),
            "num_bivalent_persister" = length(intersect(bivalent_genes,scRNA$Symbol)),
            "num_bivalent_persister_K27" =  length(intersect(bivalent_genes, depleted_K27)),
            "false_positive" = length(summits_K4_IgG_peak_significant)
        ))
    }
       
}


# Selecting qval = 0.001 & ratio = 0.15
qval = 0.001
peak_ratio = 4

fisher_pvalue_df$Significant = (fisher_pvalue_df$fisher.pvalue<qval & fisher_pvalue_df$odd_ratio > peak_ratio)

png(file.path(outputDir,"False_positive_rate_K4_K27.png"), height = 1200, width = 1500, res=300)
ggplot(tab) + geom_point(aes(x = (num_bivalent_peaks), y =false_positive,
                             color = as.character(ratio),
                             size = (p.value) , alpha=0.6)) + scale_size_continuous(
                                 breaks=c(0.001,0.005,0.01,0.05,0.1)
                             ) + theme_classic() + xlab("# Bivalent peaks") + 
    ylab("# False positive")
dev.off()

fisher_pvalue_df$gene = ""

hits <- findOverlaps(summits_K4_K27_close, annot_10k)
agg <- aggregate(annot_10k, hits, gene=paste(gene, collapse = ","))
tmp = subsetByOverlaps(summits_K4_K27_close,annot_10k)
tmp$ID = paste0(seqnames(tmp),":",start(tmp),"-",end(tmp))
fisher_pvalue_df$gene[match(tmp$ID,rownames(fisher_pvalue_df))] = agg$gene

summits_K4_K27_peak_significant = summits_K4_K27_close[fisher_pvalue_df$Significant]
summits_K4_K27_peak_significant$gene = fisher_pvalue_df$gene[fisher_pvalue_df$Significant]
rtracklayer::export(summits_K4_K27_peak_significant,file.path(outputDir,paste0("summits_K4_K27_peak_significant_",qval,"_",peak_ratio,".bed")))

write.csv(fisher_pvalue_df, file.path(outputDir,paste0("fisher_pvalue_K4_K27_peak_",qval,"_",peak_ratio,".csv")))

# summits_K4_K27_peak_significant = rtracklayer::import(file.path(outputDir,paste0("summits_K4_K27_peak_significant_",qval,"_",peak_ratio,".bed")))
# fisher_pvalue_df = read.csv(file.path(outputDir,paste0("fisher_pvalue_K4_K27_peak_",qval,"_",peak_ratio,".csv")))

bivalent_genes = unique(unlist(strsplit(fisher_pvalue_df$gene[which(fisher_pvalue_df$Significant)],split=",")))

cat("qval = ", qval, "; Ratio = ",peak_ratio," - num biv peaks = ", length(which(fisher_pvalue_df$Significant)),
    "; num depleted K27 persister bivalent = ", length(intersect(bivalent_genes, depleted_K27)),
    "; num  persister bivalent = ", length(intersect(bivalent_genes,scRNA$Symbol)),
    " False positive =", length(summits_K4_IgG_peak_significant),".\n")

(intersect(scRNA$Symbol,bivalent_genes))
(intersect(scRNA$Symbol[which(!scRNA$K27_status %in% c("No K27 Peak","Not differential"))], bivalent_genes))

scRNA$bivalent = FALSE
scRNA$bivalent[which(scRNA$Symbol %in% bivalent_genes)] = TRUE
table(scRNA$bivalent)
WriteXLS(scRNA,file.path(outputDir,"Persister_K4_K27_bivalence.xls"))

# Significant for IgG
fisher_pvalue_IgG_df$Significant = (fisher_pvalue_IgG_df$fisher.pvalue<qval & fisher_pvalue_IgG_df$odd_ratio > peak_ratio)
fisher_pvalue_IgG_df$gene = ""
hits <- findOverlaps(summits_K4_IgG_close, annot_10k)
agg <- aggregate(annot_10k, hits, gene=paste(gene, collapse = ","))
tmp = subsetByOverlaps(summits_K4_IgG_close,annot_10k)
tmp$ID = paste0(seqnames(tmp),":",start(tmp),"-",end(tmp))
fisher_pvalue_IgG_df$gene[match(tmp$ID,rownames(fisher_pvalue_IgG_df))] = agg$gene

summits_K4_IgG_peak_significant = summits_K4_IgG_close[fisher_pvalue_IgG_df$Significant]
summits_K4_IgG_peak_significant$gene = fisher_pvalue_IgG_df$gene[fisher_pvalue_IgG_df$Significant]
rtracklayer::export(summits_K4_IgG_peak_significant,file.path(outputDir,paste0("summits_K4_IgG_peak_significant_",qval,"_",peak_ratio,".bed")))
write.csv(fisher_pvalue_IgG_df, file.path(outputDir,paste0("fisher_pvalue_K4_IgG_peak_",qval,"_",peak_ratio,".csv")))

# Gene set enrichment on bivalent genes:
bivalent_genes

database <- MSIG.ls ##MSigDB
data("hg38.GeneTSS")
reflist <- union(unique(hg38.GeneTSS$Gene),bivalent_genes);length(reflist)

enrich.test <- geco.enrichmentTest(gene.sets=database,mylist=bivalent_genes,possibleIds=reflist)
enrich.test <- data.frame(Gene_set_name=rownames(enrich.test), enrich.test, check.names=FALSE)
enrich.test <- merge( subset(MSIG.gs, select=-Genes), enrich.test, by.x="Gene.Set", by.y="Gene_set_name", all.y=TRUE, sort=FALSE ) ## Get class of gene set
enrich.test <- enrich.test[order(enrich.test$`p-value`),]
enrich.test <- enrich.test[order(enrich.test$`p-value`),]
#ind <- which(enrich.test$`q-value`<= 0.1);if(!length(ind)){ind <- 1:20}
#Overexpressed  <- enrich.test[ind,]		}
Bivalent_pathways  <- enrich.test

WriteXLS(
    c("Bivalent_pathways"),
    ExcelFileName = file.path(outputDir,
                              paste0("Enrichment_test_Bivalent_pathways_K4_K27_fisher.xlsx")), 
    SheetNames = "Bivalent_pathways",
    perl = "perl", verbose = FALSE, row.names = FALSE, col.names = TRUE,
    AdjWidth = TRUE, AutoFilter = TRUE, BoldHeaderRow = TRUE, na = "",
    FreezeRow = 1, FreezeCol = 1)
