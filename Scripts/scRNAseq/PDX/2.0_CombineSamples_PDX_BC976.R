rm(list=ls())
library(here)
source(file.path(here(),"Scripts","functions.R"))
source(file.path(here(),"Scripts","global_var.R"))
library(ccRemover)

expName <- "PDX_BC408"

maindir = here()

# Directories -------------------------------------------------------------
outdir = file.path(maindir, "output","scRNAseq", "PDX_BC408");  if(!dir.exists(outdir)){dir.create(outdir)}
resSUBdir <- file.path(outdir, paste0("Unsupervised")) ; if(!file.exists(resSUBdir)){dir.create(resSUBdir)}
RDatadir <- file.path(resSUBdir,"RData") ; if(!file.exists(RDatadir)){dir.create(RDatadir)}
RDatadirSamples <- file.path(RDatadir,"RData_perSample") ; if(!file.exists(RDatadirSamples)){dir.create(RDatadirSamples)}

for(f in list.files(RDatadirSamples,full.names=TRUE)) load(f)

# Combining samples -------------------------------------------------------
all_annots <- grep("metadata",names(.GlobalEnv),value=TRUE)
annot <- do.call("rbind",mget(all_annots))
rownames(annot) <- annot$cell_id

all_counts <-grep("counts.",names(.GlobalEnv),value=TRUE)

NAMES <- vector()
for (i in 1: length(all_counts)){
  NAMES <- c(NAMES, rownames(get(all_counts[i])))
}
#if gene names start with hg19
Names_short <- sub("hg19_","",NAMES)
Names <- unique(Names_short)

SAMPLES <- names(table(annot$sample_id))

gene_metadata <- data.frame(Symbol=Names, gene_short_name=Names, cell_cycle=NA)

count_list = list()
for(i in 1:length(SAMPLES)){
  counts.sample <- get(paste0("counts.",SAMPLES[i]))
  rownames(counts.sample) = sub("hg19_","",rownames(counts.sample))
  missing = setdiff(Names, rownames(counts.sample))
  mat_missing = as(matrix(0, nrow=length(missing), ncol = ncol(counts.sample),
                          dimnames = list(missing,colnames(counts.sample))),"dgCMatrix")
  counts.sample. = rbind(counts.sample, mat_missing)
  count_list[[SAMPLES[i]]] <- counts.sample.[match(Names,rownames(counts.sample.)),]
  rm(counts.sample., counts.sample)
  gc()
}

Signal <- do.call("cbind",count_list)
gc()
if(length(which(rowSums(Signal) == 0))>0) Signal = Signal[-which(rowSums(Signal)==0),]
if(length(which(colSums(Signal) == 0))>0) Signal = Signal[,-which(colSums(Signal)==0)]
gc()

#gene_metadata <- data.frame(Symbol=Names, gene_short_name=Names, cell_cycle=NA)
gene_metadata$cell_cycle <- (gene_metadata$Symbol %in% human_cell_cycle_genes$HGNC.symbol)
rownames(gene_metadata) <- gene_metadata$Symbol
gene_metadata = gene_metadata[match(rownames(Signal), rownames(gene_metadata)),]
annot = annot[match(colnames(Signal), rownames(annot)),]
save(annot, Signal, gene_metadata, file=file.path(RDatadir,paste0(expName,".RData")))


