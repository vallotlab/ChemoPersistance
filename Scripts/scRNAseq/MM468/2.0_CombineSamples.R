library(ccRemover)
options(stringsAsFactors = F)
expName <- "MM468"

# Directories -------------------------------------------------------------

resdir <- "~/Desktop/scRNAseq_data_local/Results_MM468/"
resSUBdir <- file.path(resdir, paste0("Unsupervised_persister")) ; if(!file.exists(resSUBdir)){dir.create(resSUBdir)}
RDatadir <- file.path(resSUBdir,"RData") ; if(!file.exists(RDatadir)){dir.create(RDatadir)}
RDatadirSamples <- file.path(RDatadir,"RData_perSample") ; if(!file.exists(RDatadirSamples)){dir.create(RDatadirSamples)}

setwd(RDatadirSamples)
for(f in list.files(RDatadirSamples,full.names=TRUE)) load(f)

colnames(metadata.MM468_DMSO5_day67)[1] <- "BC_10x"
metadata.MM468_DMSO5_day67$cons_BC_lenti <- NA
metadata.MM468_DMSO5_day67$Numi_per_base <- NA
metadata.MM468_DMSO5_day67$Prop_match_cons <- NA
metadata.MM468_DMSO5_day67$side <- NA
metadata.MM468_DMSO5_day67$kept <- NA
metadata.MM468_DMSO5_day67$in_lib <- NA
metadata.MM468_DMSO5_day67$seq_rep1 <- NA
metadata.MM468_DMSO5_day67$seq_rep2 <- NA


# Combining samples -------------------------------------------------------



all_annots <-grep("metadata",names(.GlobalEnv),value=TRUE)
annot <- do.call("rbind",mget(all_annots))
rownames(annot) <- annot$cell_id

all_counts <-grep("counts.",names(.GlobalEnv),value=TRUE)

NAMES <- vector()
for (i in 1: length(all_counts)){
  NAMES <- c(NAMES, rownames(get(all_counts[i])))
}
Names <- unique(NAMES)


Signal <- matrix(0,nrow = length(Names),ncol=dim(annot)[1])
rownames(Signal) <- Names
colnames(Signal) <- annot$cell_id


SAMPLES <- names(table(annot$sample_id))

for( i in 1:length(SAMPLES)){
  counts.sample <- get(paste0("counts.",SAMPLES[i]))
  Signal[rownames(counts.sample),as.vector(annot$cell_id[annot$sample_id %in% SAMPLES[i]])] <- as.matrix(counts.sample)
}

#if gene names start with hg19
if(substr(Names[1],0,4)=="hg19") {
  Names_short <- sub("hg19_","",Names)
  gene_metadata <- data.frame(Symbol=Names_short, gene_short_name=Names_short, cell_cycle=NA)
} else {
  gene_metadata <- data.frame(Symbol=Names, gene_short_name=Names, cell_cycle=NA)
  
}

#gene_metadata <- data.frame(Symbol=Names, gene_short_name=Names, cell_cycle=NA)
gene_metadata$cell_cycle <- (gene_metadata$Symbol %in% human_cell_cycle_genes$HGNC.symbol)
rownames(gene_metadata) <- Names

save(annot,Signal,gene_metadata, file=file.path(RDatadir,paste0(expName,".RData")))
