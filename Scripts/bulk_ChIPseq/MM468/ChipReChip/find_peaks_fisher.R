library(here)

maindir= here()
resdir <- file.path(maindir,"output","bulk_ChIPseq","MM468","ChIPreChIP")
if(!dir.exists(resdir)) dir.create(resdir)
source(file.path(maindir,"Scripts","global_var.R"))


inputDir = file.path(maindir, "input","bulk_ChIPseq","MM468","ChIPreChIP")
outputDir = file.path(maindir, "output","bulk_ChIPseq","MM468","ChIPreChIP")

# K27 K4 files
primary_K27 = file.path(inputDir,"MM468bc_hu_dilution_primaryK27_K4_m10y20_H3K27me3.hg38.sorted.dedup.bam")
chipReChip_K27_K4 = file.path(inputDir,"MM468bc_hu_dilution_m10y20_H3K27me3_H3K4me3.hg38.sorted.dedup.bam")
summits_K27_K4 = rtracklayer::import(file.path(inputDir,"summits_merged_filtered_by_TSS_and_K27.bed"))
summits_K27_K4 = unique(summits_K27_K4)

# K27 IgG files
primary_K27_IgG = file.path(inputDir,"MM468bc_hu_dilution_primaryK27_IgG_m10y20_H3K27me3.hg38.sorted.dedup.bam")
chipReChip_K27_IgG= file.path(inputDir,"MM468bc_hu_dilution_m10y20_H3K27me3_IgG.hg38.sorted.dedup.bam")
summits_K27_IgG = rtracklayer::import(file.path(inputDir,"IgG_summits_merged_filtered_by_TSS_and_K27.bed"))
summits_K27_IgG = unique(summits_K27_IgG)

cat("Raw False positive rate is expected to be around ",round(length(summits_K27_IgG)/length(summits_K27_K4),3), ".\n")

# Summits K27 K4
summits_K27_K4_peak = summits_K27_K4
summits_K27_K4_around = summits_K27_K4
summits_K27_K4_close = summits_K27_K4
# Peak is +- 1000 around summit
start(summits_K27_K4_peak) = start(summits_K27_K4_peak)
end(summits_K27_K4_peak) = end(summits_K27_K4_peak)
# Around is +- 250,000 around summit
start(summits_K27_K4_around) = start(summits_K27_K4_around) - 250000
end(summits_K27_K4_around) = end(summits_K27_K4_around) + 250000
# Close is +- 10000 around summit
start(summits_K27_K4_close) = start(summits_K27_K4_close) - (width(summits_K27_K4_peak) * 5)
end(summits_K27_K4_close) = end(summits_K27_K4_close) + (width(summits_K27_K4_peak) * 5)


# Get count for chip re chip K27 K4
filt <- Rsamtools:: ScanBamParam(which = summits_K27_K4_peak)
primary_count_K27_K4_peak = Rsamtools::countBam(primary_K27, param = filt)
colnames(primary_count_K27_K4_peak) = paste0(colnames(primary_count_K27_K4_peak),"_peak")
filt <- Rsamtools:: ScanBamParam(which = summits_K27_K4_around)
primary_count_K27_K4_around = Rsamtools::countBam(primary_K27, param = filt)
colnames(primary_count_K27_K4_around) = paste0(colnames(primary_count_K27_K4_around),"_around")
filt <- Rsamtools:: ScanBamParam(which = summits_K27_K4_close)
primary_count_K27_K4_close = Rsamtools::countBam(primary_K27, param = filt)
colnames(primary_count_K27_K4_close) = paste0(colnames(primary_count_K27_K4_close),"_close")

primary_count_K27_K4 = cbind(primary_count_K27_K4_peak, primary_count_K27_K4_around, primary_count_K27_K4_close)

# Get count for primary K27 (K4)
filt <- Rsamtools:: ScanBamParam(which = summits_K27_K4_peak)
chiprechip_count_K27_K4_peak = Rsamtools::countBam(chipReChip_K27_K4, param = filt)
colnames(chiprechip_count_K27_K4_peak) = paste0(colnames(chiprechip_count_K27_K4_peak),"_peak")
filt <- Rsamtools:: ScanBamParam(which = summits_K27_K4_around)
chiprechip_count_K27_K4_around = Rsamtools::countBam(chipReChip_K27_K4, param = filt)
colnames(chiprechip_count_K27_K4_around) = paste0(colnames(chiprechip_count_K27_K4_around),"_around")
filt <- Rsamtools:: ScanBamParam(which = summits_K27_K4_close)
chiprechip_count_K27_K4_close = Rsamtools::countBam(chipReChip_K27_K4, param = filt)
colnames(chiprechip_count_K27_K4_close) = paste0(colnames(chiprechip_count_K27_K4_close),"_close")

chiprechip_count_K27_K4 = cbind(chiprechip_count_K27_K4_peak, chiprechip_count_K27_K4_around,chiprechip_count_K27_K4_close)


fisher_pvalue_mat = matrix(0,nrow = nrow(chiprechip_count_K27_K4))
odd_ratio_mat = matrix(0,nrow = nrow(chiprechip_count_K27_K4))
signal_compared_to_close = matrix(0,nrow = nrow(chiprechip_count_K27_K4))
rownames(fisher_pvalue_mat) = paste0(chiprechip_count_K27_K4$space_peak,":",chiprechip_count_K27_K4$start_peak,"-",chiprechip_count_K27_K4$end_peak)
rownames(odd_ratio_mat) = paste0(chiprechip_count_K27_K4$space_peak,":",chiprechip_count_K27_K4$start_peak,"-",chiprechip_count_K27_K4$end_peak)

for(i in seq_len(nrow(fisher_pvalue_mat))){

    contingency_table = data.frame(a= c(chiprechip_count_K27_K4$records_peak[i], chiprechip_count_K27_K4$records_around[i]-chiprechip_count_K27_K4$records_peak[i]),
                                   b = c(primary_count_K27_K4$records_peak[i], primary_count_K27_K4$records_around[i]-primary_count_K27_K4$records_peak[i]))
    
    fisher_pvalue_mat[i,] = fisher.test(contingency_table,alternative = "greater")$p.value
    odd_ratio_mat[i,] = (contingency_table$a[1] / contingency_table$a[2]) / (contingency_table$b[1] / contingency_table$b[2])
    signal_compared_to_close[i,] = chiprechip_count_K27_K4$records_peak[i] / (chiprechip_count_K27_K4$records_peak[i] + chiprechip_count_K27_K4$records_close[i])
}

fisher_pvalue_df = as.data.frame(as.matrix(fisher_pvalue_mat))
fisher_pvalue_df$odd_ratio = as.numeric(odd_ratio_mat)
fisher_pvalue_df$ratio_peak_close = as.numeric(signal_compared_to_close)

annot_10k = read.table(gzfile(file.path(maindir,"annotation",
                                        "gencode.v34.annotation.transcriptTSS_10k.bed.gz")),
                       sep="\t", header = F)
colnames(annot_10k) = c("chr","start","end","transcripts","gene","strand")
annot_10k = as(annot_10k,"GRanges")

# Intersect genes K4 K27 & persister genes
# scRNA = readxl::read_xlsx(file.path(maindir,"output","scRNAseq","MM468","Persister","Supervised","Tables","Differential_analysis_Limma_logFC_1.58.xlsx"))
scRNA = readxl::read_xls(file.path(maindir,"output","scRNAseq","MM468","Persister","ComparisonChIPseq","K27","persister_K27_status.xls"))
depleted_K27 = scRNA$Symbol[which(!scRNA$K27_status %in% c("No K27 Peak","Not differential"))]
# depleted_K27 = scRNA$Symbol[which(scRNA$K27_status %in% c("Depleted FC < -3"))]

# Getting Pvalues for IgG
# IgG
# Summits K27 K4
summits_K27_IgG_peak = summits_K27_IgG
summits_K27_IgG_around = summits_K27_IgG
summits_K27_IgG_close = summits_K27_IgG

# Peak is +- 1000 around summit
start(summits_K27_IgG_peak) = start(summits_K27_IgG_peak)
end(summits_K27_IgG_peak) = end(summits_K27_IgG_peak) 
# Around is +- 250,000 around summit
start(summits_K27_IgG_around) = start(summits_K27_IgG_around) - 250000
end(summits_K27_IgG_around) = end(summits_K27_IgG_around) + 250000
# Close is +- 10000 around summit
start(summits_K27_IgG_close) = start(summits_K27_IgG_close) - (width(summits_K27_IgG_peak) * 5)
end(summits_K27_IgG_close) = end(summits_K27_IgG_close) + (width(summits_K27_IgG_peak) * 5)


# Get count for chip re chip K27 K4
filt <- Rsamtools:: ScanBamParam(which = summits_K27_IgG_peak)
primary_count_K27_IgG_peak = Rsamtools::countBam(primary_K27_IgG, param = filt)
colnames(primary_count_K27_IgG_peak) = paste0(colnames(primary_count_K27_IgG_peak),"_peak")
filt <- Rsamtools:: ScanBamParam(which = summits_K27_IgG_around)
primary_count_K27_IgG_around = Rsamtools::countBam(primary_K27_IgG, param = filt)
colnames(primary_count_K27_IgG_around) = paste0(colnames(primary_count_K27_IgG_around),"_around")
filt <- Rsamtools:: ScanBamParam(which = summits_K27_IgG_close)
primary_count_K27_IgG_close = Rsamtools::countBam(primary_K27_IgG, param = filt)
colnames(primary_count_K27_IgG_close) = paste0(colnames(primary_count_K27_IgG_close),"_close")

primary_count_K27_IgG = cbind(primary_count_K27_IgG_peak,
                              primary_count_K27_IgG_around,
                              primary_count_K27_IgG_close)

# Get count for primary K27 (K4)
filt <- Rsamtools:: ScanBamParam(which = summits_K27_IgG_peak)
chiprechip_count_K27_IgG_peak = Rsamtools::countBam(chipReChip_K27_IgG, param = filt)
colnames(chiprechip_count_K27_IgG_peak) = paste0(colnames(chiprechip_count_K27_IgG_peak),"_peak")
filt <- Rsamtools:: ScanBamParam(which = summits_K27_IgG_around)
chiprechip_count_K27_IgG_around = Rsamtools::countBam(chipReChip_K27_IgG, param = filt)
colnames(chiprechip_count_K27_IgG_around) = paste0(colnames(chiprechip_count_K27_IgG_around),"_around")
filt <- Rsamtools:: ScanBamParam(which = summits_K27_IgG_close)
chiprechip_count_K27_IgG_close = Rsamtools::countBam(chipReChip_K27_IgG, param = filt)
colnames(chiprechip_count_K27_IgG_close) = paste0(colnames(chiprechip_count_K27_IgG_close),"_close")

chiprechip_count_K27_IgG = cbind(chiprechip_count_K27_IgG_peak,
                                 chiprechip_count_K27_IgG_around,
                                 chiprechip_count_K27_IgG_close)

fisher_pvalue_IgG_mat = matrix(0,nrow = nrow(chiprechip_count_K27_IgG))
odd_ratio_mat = matrix(0,nrow = nrow(chiprechip_count_K27_IgG))
signal_compared_to_close = matrix(0,nrow = nrow(chiprechip_count_K27_IgG))
rownames(fisher_pvalue_IgG_mat) = paste0(chiprechip_count_K27_IgG$space_peak,":",chiprechip_count_K27_IgG$start_peak,"-",chiprechip_count_K27_IgG$end_peak)
rownames(odd_ratio_mat) = paste0(chiprechip_count_K27_IgG$space_peak,":",chiprechip_count_K27_IgG$start_peak,"-",chiprechip_count_K27_IgG$end_peak)

for(i in seq_len(nrow(fisher_pvalue_IgG_mat))){
    
    contingency_table = data.frame(a= c(chiprechip_count_K27_IgG$records_peak[i], chiprechip_count_K27_IgG$records_around[i]-chiprechip_count_K27_IgG$records_peak[i]),
                                   b = c(primary_count_K27_IgG$records_peak[i], primary_count_K27_IgG$records_around[i]-primary_count_K27_IgG$records_peak[i]))
    
    fisher_pvalue_IgG_mat[i,] = fisher.test(contingency_table,alternative = "greater")$p.value
    odd_ratio_mat[i,] = (contingency_table$a[1] / contingency_table$a[2]) / (contingency_table$b[1] / contingency_table$b[2])
    signal_compared_to_close[i,] = chiprechip_count_K27_IgG$records_peak[i] / (chiprechip_count_K27_IgG$records_peak[i] + chiprechip_count_K27_IgG$records_close[i])
}

fisher_pvalue_IgG_df = as.data.frame(as.matrix(fisher_pvalue_IgG_mat))
fisher_pvalue_IgG_df$odd_ratio = as.numeric(odd_ratio_mat)
fisher_pvalue_IgG_df$ratio_peak_close = as.numeric(signal_compared_to_close)

colnames(fisher_pvalue_df)[1] = "fisher.pvalue"
colnames(fisher_pvalue_IgG_df)[1] = "fisher.pvalue"

fisher_pvalue_df$qvalue = p.adjust(fisher_pvalue_df$fisher.pvalue,method = "BH")
fisher_pvalue_IgG_df$qvalue = p.adjust(fisher_pvalue_IgG_df$fisher.pvalue,method = "BH")

tab=data.frame()
for(pval in c(0.1,0.05,0.01,0.005,0.001)){
    for(ratio in c(0.1,0.15,0.2,0.25)){
        # select signif peaks for K27 K4
        fisher_pvalue_df$Significant = (fisher_pvalue_df$fisher.pvalue<qval & fisher_pvalue_df$ratio_peak_close > ratio)
        
        summits_K27_K4_peak_significant = summits_K27_K4_peak[fisher_pvalue_df$Significant]
        summits_K27_K4_peak_significant$gene = ""
        summits_K27_K4_peak_significant$ID = paste0(seqnames(summits_K27_K4_peak_significant),":",
                                                    start(summits_K27_K4_peak_significant),"-",
                                                    end(summits_K27_K4_peak_significant))
        
        hits <- findOverlaps(summits_K27_K4_peak_significant, annot_10k)
        agg <- aggregate(annot_10k, hits, gene=paste(gene, collapse = ","))
        
        summits_K27_K4_peak_significant$gene[match(subsetByOverlaps(summits_K27_K4_peak_significant,annot_10k)$ID,
                                                   summits_K27_K4_peak_significant$ID)] = agg$gene
        
        # select significant peaks for K27 IgG
        fisher_pvalue_IgG_df$Significant = (fisher_pvalue_IgG_df$fisher.pvalue<pval & fisher_pvalue_IgG_df$ratio_peak_close > ratio)
        
        if(length(which(fisher_pvalue_IgG_df$Significant))>0){
        summits_K27_IgG_peak_significant = summits_K27_IgG_peak[fisher_pvalue_IgG_df$Significant]
        summits_K27_IgG_peak_significant$gene = ""
        summits_K27_IgG_peak_significant$ID = paste0(seqnames(summits_K27_IgG_peak_significant),":",
                                                    start(summits_K27_IgG_peak_significant),"-",
                                                    end(summits_K27_IgG_peak_significant))
        
        hits <- findOverlaps(summits_K27_IgG_peak_significant, annot_10k)
        agg <- aggregate(annot_10k, hits, gene=paste(gene, collapse = ","))
        
        summits_K27_IgG_peak_significant$gene[match(subsetByOverlaps(summits_K27_IgG_peak_significant,annot_10k)$ID,summits_K27_IgG_peak_significant$ID)] = agg$gene
        } else {
            summits_K27_IgG_peak_significant = NULL
        }
        
        bivalent_genes = unique(unlist(strsplit(summits_K27_K4_peak_significant$gene,split=",")))
        cat("qval = ", pval, "; Ratio = ",ratio," - num biv peaks = ", length(which(fisher_pvalue_df$Significant)),
              "; num depleted K27 persister bivalent = ", length(intersect(bivalent_genes, depleted_K27)),
              "; num  persister bivalent = ", length(intersect(bivalent_genes,scRNA$Symbol)),
              "False positive = ",length(summits_K27_IgG_peak_significant),".\n")
        # rtracklayer::export(summits_K27_K4_peak_significant,file.path(outputDir,paste0("summits_K27_K4_peak_significant_",qval,"_",ratio,".bed")))
        tab = rbind(tab, data.frame(
            "p.value" = pval, "ratio" = ratio, "num_bivalent_peaks" = length(which(fisher_pvalue_df$Significant)),
            "num_bivalent_persister" = length(intersect(bivalent_genes,scRNA$Symbol)),
            "num_bivalent_persister_K27" =  length(intersect(bivalent_genes, depleted_K27)),
            "false_positive" = length(summits_K27_IgG_peak_significant)
        ))
    }
}


# Selecting qval = 0.001 & ratio = 0.15
qval = 0.001
ratio = 0.15

fisher_pvalue_df$Significant = (fisher_pvalue_df$qvalue<qval & fisher_pvalue_df$ratio_peak_close > ratio)

png(file.path(outputDir,"False_positive_rate.png"), height = 1200, width = 1500, res=300)
ggplot(tab) + geom_point(aes(x = (num_bivalent_peaks), y =false_positive,
                             color = as.character(ratio),
                             size = (p.value) , alpha=0.6)) + scale_size_continuous(
                                 breaks=c(0.001,0.005,0.01,0.05,0.1)
                             ) + theme_classic() + xlab("# Bivalent peaks") + 
    ylab("# False positive")
dev.off()


fisher_pvalue_df$gene = ""

hits <- findOverlaps(summits_K27_K4_peak, annot_10k)
agg <- aggregate(annot_10k, hits, gene=paste(gene, collapse = ","))
tmp = subsetByOverlaps(summits_K27_K4_peak,annot_10k)
tmp$ID = paste0(seqnames(tmp),":",start(tmp),"-",end(tmp))
fisher_pvalue_df$gene[match(tmp$ID,rownames(fisher_pvalue_df))] = agg$gene

summits_K27_K4_peak_significant = summits_K27_K4_peak[fisher_pvalue_df$Significant]
summits_K27_K4_peak_significant$gene = fisher_pvalue_df$gene[fisher_pvalue_df$Significant]
rtracklayer::export(summits_K27_K4_peak_significant,file.path(outputDir,paste0("summits_K27_K4_peak_significant_",qval,"_",ratio,".bed")))

write.csv(fisher_pvalue_df, file.path(outputDir,paste0("fisher_pvalue_K27_K4_peak_",qval,"_",ratio,".csv")))

bivalent_genes = unique(unlist(strsplit(summits_K27_K4_peak_significant$gene,split=",")))

cat("qval = ", qval, "; Ratio = ",ratio," - num biv peaks = ", length(which(fisher_pvalue_df$Significant)),
    "; num depleted K27 persister bivalent = ", length(intersect(bivalent_genes, depleted_K27)),
    "; num  persister bivalent = ", length(intersect(bivalent_genes,scRNA$Symbol)),
    " False positive =", length(summits_K27_IgG_peak_significant),".\n")

(intersect(scRNA$Symbol,bivalent_genes))
(intersect(scRNA$Symbol[which(!scRNA$K27_status %in% c("No K27 Peak","Not differential"))], bivalent_genes))
scRNA$bivalent = FALSE
scRNA$bivalent[which(scRNA$Symbol %in% bivalent_genes)] = TRUE
table(scRNA$bivalent)
WriteXLS(scRNA,file.path(outputDir,"Persister_K27_bivalence.xls"))

# Significant for IgG
fisher_pvalue_IgG_df$Significant = (fisher_pvalue_IgG_df$qvalue<qval & fisher_pvalue_IgG_df$ratio_peak_close > ratio)
fisher_pvalue_IgG_df$gene = ""
hits <- findOverlaps(summits_K27_IgG_peak, annot_10k)
agg <- aggregate(annot_10k, hits, gene=paste(gene, collapse = ","))
tmp = subsetByOverlaps(summits_K27_IgG_peak,annot_10k)
tmp$ID = paste0(seqnames(tmp),":",start(tmp),"-",end(tmp))
fisher_pvalue_IgG_df$gene[match(tmp$ID,rownames(fisher_pvalue_IgG_df))] = agg$gene

summits_K27_IgG_peak_significant = summits_K27_IgG_peak[fisher_pvalue_IgG_df$Significant]
summits_K27_IgG_peak_significant$gene = fisher_pvalue_IgG_df$gene[fisher_pvalue_IgG_df$Significant]
rtracklayer::export(summits_K27_IgG_peak_significant,file.path(outputDir,paste0("summits_K27_K4_peak_significant_",qval,"_",ratio,".bed")))

write.csv(fisher_pvalue_IgG_df, file.path(outputDir,paste0("fisher_pvalue_K27_IgG_peak_",qval,"_",ratio,".csv")))

# Gene set enrichment on bivalent genes:
bivalent_genes

database <- MSIG.ls ##MSigDB
reflist <- union(unique(hg38.GeneTSS$gene),bivalent_genes);length(reflist)

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
                              paste0("Enrichment_test_Bivalent_pathways.xlsx")), 
    SheetNames = "Bivalent_pathways",
    perl = "perl", verbose = FALSE, row.names = FALSE, col.names = TRUE,
    AdjWidth = TRUE, AutoFilter = TRUE, BoldHeaderRow = TRUE, na = "",
    FreezeRow = 1, FreezeCol = 1)
