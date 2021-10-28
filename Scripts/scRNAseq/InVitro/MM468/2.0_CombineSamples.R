rm(list = ls()); gc()

library(here)
source(file.path(here(),"Scripts","functions.R"))
source(file.path(here(),"Scripts","global_var.R"))
library(ccRemover)

expName <- "MM468_UNC_5FU"

maindir = here()

# Directories -------------------------------------------------------------
outdir = file.path(maindir, "output","scRNAseq", "MM468");  if(!dir.exists(outdir)){dir.create(outdir)}
resdir <- file.path(outdir, "UNC_5FU"); if(!dir.exists(resdir)){dir.create(resdir)}
resSUBdir <- file.path(resdir, paste0("Unsupervised")) ; if(!file.exists(resSUBdir)){dir.create(resSUBdir)}
RDatadir <- file.path(resSUBdir,"RData") ; if(!file.exists(RDatadir)){dir.create(RDatadir)}
RDatadirSamples <- file.path(maindir, "output","scRNAseq", "MM468","Persister","Unsupervised","RData", "RData_perSample")  ; if(!file.exists(RDatadirSamples)){dir.create(RDatadirSamples)}

for(f in list.files(RDatadirSamples,full.names=TRUE, include.dirs = F)) load(f)

# Combining samples -------------------------------------------------------
all_annots <- grep("metadata",names(.GlobalEnv),value=TRUE)
annot <- do.call("rbind", mget(all_annots))
rownames(annot) <- annot$cell_id

all_counts <-grep("counts.",names(.GlobalEnv),value=TRUE)

NAMES <- vector()
for (i in 1: length(all_counts)){
  NAMES <- c(NAMES, rownames(get(all_counts[i])))
}
Names <- unique(NAMES)

gc()
SAMPLES <- names(table(annot$sample_id))

gene_metadata <- data.frame(Symbol=Names, gene_short_name=Names, cell_cycle=NA)

count_list = list()
for(i in 1:length(SAMPLES)){
  counts.sample <- get(paste0("counts.",SAMPLES[i]))
  missing = setdiff(Names, rownames(counts.sample))
  mat_missing = as(matrix(0, nrow=length(missing), ncol = ncol(counts.sample),
                          dimnames = list(missing,colnames(counts.sample))),"dgCMatrix")
  counts.sample. = rbind(counts.sample, mat_missing)
  count_list[[SAMPLES[i]]] <- counts.sample.[match(Names,rownames(counts.sample.)),]
  rm(counts.sample., counts.sample)
  gc()
}

Signal <- do.call("cbind",count_list)
if(length(which(rowSums(Signal) == 0))>0) Signal = Signal[-which(rowSums(Signal)==0),]
if(length(which(colSums(Signal) == 0))>0) Signal = Signal[,-which(colSums(Signal)==0)]
gc()

gene_metadata$cell_cycle <- (gene_metadata$Symbol %in% human_cell_cycle_genes$HGNC.symbol)
rownames(gene_metadata) <- Names

gene_metadata = gene_metadata[match(rownames(Signal), rownames(gene_metadata)),]
annot = annot[match(colnames(Signal), rownames(annot)),]
dim(Signal)
dim(annot)
dim(gene_metadata)

samples = unique(annot$sample_id)
samples_test_TwoOfEach = c("MM468_chemonaive", "MM468_chemonaive_2", "MM468_5FU7_day21", "MM468_5FU1_day33", "MM468_5FU7_UNC_day21", "MM468_5FU1_UNC_day33")
samples_test_All = c("MM468_chemonaive", "MM468_chemonaive_2", "MM468_5FU7_day21", "MM468_5FU1_day33", "MM468_5FU3_day50","MM468_5FU2_day67", "MM468_5FU7_UNC_day21", "MM468_5FU1_UNC_day33")

Signal_safe = Signal
annot_safe = annot

# Two of Each
cells = annot_safe$cell_id[which(annot_safe$sample_id %in% samples_test_TwoOfEach)]
annot = annot_safe[match(cells, annot_safe$cell_id),]
Signal = Signal_safe[,match(cells, colnames(Signal_safe))]
dim(Signal)
dim(annot)
table(annot$sample_id)
save(annot,Signal,gene_metadata, file=file.path(RDatadir,paste0(expName,"_TwoOfEach.RData")))

# All
cells = annot_safe$cell_id[which(annot_safe$sample_id %in% samples_test_All)]
annot = annot_safe[match(cells, annot_safe$cell_id),]
Signal = Signal_safe[,match(cells, colnames(Signal_safe))]
dim(Signal)
dim(annot)
table(annot$sample_id)
save(annot,Signal,gene_metadata, file=file.path(RDatadir,paste0(expName,"_All.RData")))

