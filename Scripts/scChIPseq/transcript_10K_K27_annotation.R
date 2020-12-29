# Create 10k transcript TSS K27 annotation 
library(here)
maindir = here()
source(file.path(maindir,"Scripts","global_var.R"))
annotDir = file.path(maindir, "annotation")

# Load Transcript 10k annotation (52138 10k regions around transcript TSS)
annot_transcript_10k = read.table(gzfile(file.path(annotDir, "gencode.v34.annotation.transcriptTSS_10k.bed.gz")), sep ="\t")
colnames(annot_transcript_10k) <- c("chr","start","end","transcripts","gene","strand")
annot_transcript_10k = as(annot_transcript_10k,"GRanges")

# Take gtf gencode v34 as "TSS" annotation : take only 1bp of TSS
annot_transcript_TSS = read.table(gzfile(file.path(maindir,"annotation","gencode.v34.annotation_raw.bed.gz")), sep="\t")
colnames(annot_transcript_TSS) = c("chr","start","end","transcripts","gene","strand")
annot_transcript_TSS = as(annot_transcript_TSS,"GRanges")

# Take TSS +- 500bp
annot_transcript_TSS. = annot_transcript_TSS
end(annot_transcript_TSS.) = ifelse(as.character(annot_transcript_TSS@strand)=="+", start(annot_transcript_TSS) + 500, end(annot_transcript_TSS) + 500)
start(annot_transcript_TSS.) = ifelse(as.character(annot_transcript_TSS@strand)=="+", start(annot_transcript_TSS) - 500, end(annot_transcript_TSS) + 500)
annot_transcript_TSS = annot_transcript_TSS.

# Load K27 consensus peaks (9568 variable peak sizes around transcript TSS)
annot_K27 = rtracklayer::import(file.path(annotDir,"MM468_peaks_K27.bed.gz"))
start(annot_K27) = start(annot_K27) - 1

# Annotate transcripts with info: is this TSS overlapping a K27 peak ?
hits = findOverlaps(annot_transcript_TSS, annot_K27,)
annot_transcript_TSS$K27_peak = FALSE
annot_transcript_TSS$K27_peak[queryHits(hits)] = TRUE
table(annot_transcript_TSS$K27_peak)

# Go back to 10k region around TSS:
annot_transcript_10k_split_by_transcript = as.data.frame(annot_transcript_10k)
annot_transcript_10k_split_by_transcript$row.id = seq_len(nrow(annot_transcript_10k_split_by_transcript))
annot_transcript_10k_split_by_transcript = annot_transcript_10k_split_by_transcript %>% tidyr::separate_rows(transcripts, sep=",")

transcript_to_remove = setdiff(annot_transcript_TSS$transcripts, annot_transcript_10k_split_by_transcript$transcripts)
annot_transcript_TSS = annot_transcript_TSS[-which(annot_transcript_TSS$transcripts %in% transcript_to_remove),]

# Check if the two tables have same transcripts 
all.equal(annot_transcript_TSS$transcripts, annot_transcript_10k_split_by_transcript$transcripts)

# Check that transcript are unique
length(unique(annot_transcript_10k_split_by_transcript$transcripts)) == nrow(annot_transcript_10k_split_by_transcript) # == length(annot_transcript_TSS)

annot_transcript_10k_split_by_transcript$K27_peak = annot_transcript_TSS$K27_peak

annot_transcript_10k_split_by_transcript

# If one of the transcripts mapping back to the 10k approximate transcript is overlapping, put TRUE in K27_peak
agg <- aggregate(annot_transcript_10k_split_by_transcript[,c("row.id","K27_peak")],
                 list(hits = annot_transcript_10k_split_by_transcript$row.id), any)
annot_transcript_10k$K27_peak = agg$K27_peak

annot_transcript_10k$transcripts[grep("CNKSR2",annot_transcript_10k$transcripts)]
annot_transcript_10k$K27_peak[grep("CNKSR2",annot_transcript_10k$transcripts)]

table(annot_transcript_10k$K27_peak)

# Load bulk K27 10k-based differential analysis
DiffAnalysis_K27_10k = read.csv(
    file.path(maindir,"output","bulk_ChIPseq","MM468","K27_transcripts_10k",
              "Supervised","Tables", "Supervised_analysis_Limma.csv"), sep=";")[,c(1,2,3,6,7,8,10,11)]
DiffAnalysis_K27_10k$log2FC.X5FU2_3_5 = as.numeric(gsub(",",".",DiffAnalysis_K27_10k$log2FC.X5FU2_3_5))
DiffAnalysis_K27_10k$pval.X5FU2_3_5 = as.numeric(gsub(",",".",DiffAnalysis_K27_10k$pval.X5FU2_3_5))
DiffAnalysis_K27_10k$qval.X5FU2_3_5 = as.numeric(gsub(",",".",DiffAnalysis_K27_10k$qval.X5FU2_3_5))

annot_transcript_10k_df = as.data.frame(annot_transcript_10k)
colnames(annot_transcript_10k_df)[c(1,7)] = c("chr","Gene")
annot_transcript_10k_df = dplyr::left_join(annot_transcript_10k_df, DiffAnalysis_K27_10k)

# Add log2FC & qvalue information to Transcripts TSS associated with a significantly differential K27 peak
annot_transcript_10k$log2FC_K27 = annot_transcript_10k_df$log2FC.X5FU2_3_5
annot_transcript_10k$qval_K27 = annot_transcript_10k_df$qval.X5FU2_3_5

# Filter out too lowly covered 10k loci in Untreated:
# Load bulk matrices K27 in 10k annot
resultZ <- read.csv(unzip(file.path(maindir,"input","bulk_ChIPseq","MM468","results_MM468_all_transcripts10k_v2.zip"),
                          exdir = tempdir()))

checkRange = as(setNames(resultZ[,1:3],c("chr",'start',"end")),"GRanges")
all.equal(ranges(checkRange),ranges(annot_transcript_10k))

#Preparing raw count and normalized count matrices
RZ <- (resultZ[,-c(1:3,grep(pattern = "input",colnames(resultZ)))]) #RZ matrix is the raw count matrix
untreated_samples = c("MM468_2_p15_J113_DMSO_K27", "MM468_2_p16_J113_DMSO_K27", "MM468_3_J77_DMSO_K27", "MM468_5_J67_DMSO_K27")
RZ = RZ[,untreated_samples]

#create RPKM matrix by normalizing per library size and peak length
RZ_RPKM <- (10^9*t(t((RZ+1)/10001)/colSums(RZ)))

annot_transcript_10k$RPM = log2(rowMeans(RZ_RPKM))
    
annot10k_K27_ovlp = as.data.frame(annot_transcript_10k)[which(annot_transcript_10k$K27_peak),]
annot10k_K27_byGene = annot10k_K27_ovlp %>% group_by(gene)  %>%
    slice_max(order_by = abs(log2FC_K27), n = 1) # Select only the transcript associated with top differential K27 log2FC 

hist(log2(rowMeans(RZ_RPKM)), breaks=150)
threshold = log2(3)
abline(v = threshold)
annot_transcript_10k. = annot_transcript_10k
too_low = which(log2(rowMeans(RZ_RPKM)) < threshold)
annot_transcript_10k.$K27_peak[too_low] = FALSE

annot10k_K27_ovlp. = as.data.frame(annot_transcript_10k.)[which(annot_transcript_10k.$K27_peak),]
annot10k_K27_byGene. = annot10k_K27_ovlp. %>% group_by(gene)  %>%
    slice_max(order_by = abs(log2FC_K27), n = 1) # Select only the transcript associated with top differential K27 log2FC 

annot_transcript_10k$K27_peak[too_low] = FALSE

# Categorize each 10k TSS into:
# - Not overlapping K27 peak
# - Overlapping a non significantly differential K27 peak
# - Overlapping a significantly enriched K27 peak
# - Overlapping a significantly depleted K27 peak
annot_transcript_10k$K27_status = ifelse(!annot_transcript_10k$K27_peak, "Not Overlapping",
                                         ifelse(is.nan(annot_transcript_10k$log2FC_K27), "Not differential",
                                         ifelse(annot_transcript_10k$log2FC_K27 > log2(3) & annot_transcript_10k$qval_K27 < 0.1, "Enriched",
                                                ifelse(annot_transcript_10k$log2FC_K27 < -log2(3) & annot_transcript_10k$qval_K27 < 0.1, "Depleted",
                                                       "Not differential"))))
table(annot_transcript_10k$K27_status)

depleted10k = annot_transcript_10k[which(annot_transcript_10k$K27_status == "Depleted"),]

annot_transcript_10k[grep("COL12A1",annot_transcript_10k$gene),]
depleted10k[grep("COL12A1",depleted10k$transcripts),]

write.table(as.data.frame(annot_transcript_10k)[,-c(4,5)], sep="\t", file = file.path(annotDir,"gencode.v34.transcripts10k_K27.tsv"))
write.table(as.data.frame(annot_transcript_10k)[,-c(4,5)], sep="\t",
            file = file.path(annotDir,"gencode.v34.transcripts10k_K27.bed"),
            col.names = F,row.names = F, quote = F)

# Annotate K27 peak annotation with correct gene information
annot10k_K27 = read.table(file.path(annotDir,"gencode.v34.transcripts10k_K27.tsv"), sep = "\t")
# Keep only transcipts overlapping K27 peaks
annot10k_K27_ovlp = annot10k_K27[which(annot10k_K27$K27_peak),]
annot10k_K27_byGene = annot10k_K27_ovlp %>% group_by(gene) %>%  slice_max(order_by = abs(log2FC_K27), n = 1) # Select only the transcript associated with top differential K27 log2FC 
annot10k_K27_byGene = as(annot10k_K27_byGene,"GRanges")

annot_K27_gene = annot_K27
hits = findOverlaps(annot_K27_gene, annot10k_K27_byGene)
agg <- aggregate(annot10k_K27_byGene,hits, gene = paste(gene, collapse=","))

annot_K27_gene$gene = NA
annot_K27_gene$gene[unique(queryHits(hits))] = agg$gene

annot_K27_gene[grep("COL12A1",annot_K27_gene$gene)]
write.table(as.data.frame(annot_K27_gene)[,-c(4,5,6)], sep="\t", file = file.path(annotDir,"MM468_peaks_K27_gene.bed"),
            col.names = F, row.names = F, quote = F)
