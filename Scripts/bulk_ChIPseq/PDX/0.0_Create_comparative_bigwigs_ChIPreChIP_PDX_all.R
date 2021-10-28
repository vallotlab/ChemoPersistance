library(here)

maindir= here()
source(file.path(maindir, "Scripts", "global_var.R"))
source(file.path(maindir, "Scripts", "functions.R"))

inputDir_coverage <- file.path(maindir,"input","bulk_ChIPseq","PDX","ChIPreChIP","BigWigs_Compare","Raw")
outputDir_compare <- file.path(maindir,"input","bulk_ChIPseq","PDX","ChIPreChIP","BigWigs_Compare")
resDir <- file.path(maindir,"output","bulk_ChIPseq","Snapshots")

K4_H3K27me3 = list.files(inputDir_coverage, pattern = "K27")
K4_IgG = list.files(inputDir_coverage, pattern = "IgG")

models = gsub("_.*","",reChIP_K4_IgG)

names(K4_H3K27me3) = models
names(K4_IgG) = models

library(GenomicRanges)
library(ChromSCape)
data("hg38.chromosomes")
canonical_chr <- hg38.chromosomes
canonical_chr$start = 1
canonical_chr <- as(canonical_chr, "GRanges")
bins <- unlist(GenomicRanges::tileGenome(setNames(GenomicRanges::width(canonical_chr), 
                                                  GenomicRanges::seqnames(canonical_chr)), tilewidth = 30))
library(BiocParallel)
BPPARAM = bpparam()
BPPARAM$workers = 2
BiocParallel::register(BPPARAM = BPPARAM)
gc()
for(i in models){
    K4_K27_gr = rtracklayer::import(file.path(inputDir_coverage, K4_H3K27me3[[i]]))
    K4_IgG_gr = rtracklayer::import(file.path(inputDir_coverage, K4_IgG[[i]]))
    
    print(i)
    # K4 K27
    bins. = bins
    hits <- GenomicRanges::findOverlaps(bins., K4_K27_gr)
    bins.$score= 0
    bins.$score[queryHits(hits)] = K4_K27_gr$score[subjectHits(hits)]
    gc()
    print("Running scores...")
    score = ChromSCape:::smoothBin(bins.$score, nb_bins = 10)
    bins.$score = score
    bins. = bins.[-which(bins.$score==0)]
    rtracklayer::export.bw(bins., file.path(outputDir_compare,paste0(i,"_reChIP_H3K4me3_H3K27me3_vs_primary_H3K4me3.bw")))
    print("Done K27 !")
    
    gc()
    # K4 IgG
    bins. = bins
    gc()
    hits <- GenomicRanges::findOverlaps(bins., K4_IgG_gr)
    bins.$score= 0
    bins.$score[queryHits(hits)] = K4_IgG_gr$score[subjectHits(hits)]
    gc()
    BiocParallel::register(BPPARAM = BPPARAM)
    print("Running scores...")
    score = ChromSCape:::smoothBin(bins.$score, nb_bins = 10)
    bins.$score = score
    bins. = bins.[-which(bins.$score==0)]
    rtracklayer::export.bw(bins., file.path(outputDir_compare,paste0(i,"_reChIP_H3K4me3_IgG_vs_primary_H3K4me3.bw")))
    gc()
    print("Done IgG !")
}
