#Enrichment scores (p and qvalue) need to be kept  all gene lists for log2FC 1.58 for persister_6 and C2
# regenerate Overexpressed analysis
library(here)
source(file.path(here(),"Scripts","global_var.R"))
maindir = here()
log2FC_threshold = log2(3)
outdir <- file.path(maindir,"output","scRNAseq","MM468_PDX")
if(!dir.exists(outdir)) dir.create(outdir)

mycolramp <- c("white",viridis(n=4))

Overexpressed_PDX_genesets  <- readxl::read_xlsx(file.path(maindir,"output","scRNAseq","PDX","Supervised",
                                              "Tables","Gene_set_analysis","Enrichment_test_persister_vs_UNT_vs_UNT_MSigDB_logFC1.58.xlsx"),
                                    sheet = 1)

colnames(Overexpressed_PDX_genesets) <- paste0(colnames(Overexpressed_PDX_genesets),"_PDX")
colnames(Overexpressed_PDX_genesets)[1] <- "Gene.Set"

##########################
#MM468
##########################
Overexpressed_MM468_genesets  <- readxl::read_xlsx(file.path(maindir,"output","scRNAseq","MM468","Persister","Supervised",
                                                           "Tables","Gene_set_analysis","Enrichment_test_C2_pers_vs_DMSO_MSigDB_logFC1.58.xlsx"),
                                                 sheet = 1)
colnames(Overexpressed_MM468_genesets) <- paste0(colnames(Overexpressed_MM468_genesets),"_MM468")
colnames(Overexpressed_MM468_genesets)[1] <- "Gene.Set"

#######################
##Comparison Gene Lists
#######################

Both_genelist <- inner_join(Overexpressed_MM468_genesets,Overexpressed_PDX_genesets,by="Gene.Set")
Both_genelist$log10_qvalue_MM468 <- -log10(Both_genelist$`q-value_MM468`)
Both_genelist$log10_qvalue_PDX <- -log10(Both_genelist$`q-value_PDX`)

correlation <- cor.test(Both_genelist$log10_qvalue_MM468,Both_genelist$log10_qvalue_PDX) 
lmper <- lm(Both_genelist$log10_qvalue_PDX~Both_genelist$log10_qvalue_MM468) 

png("-log10 q value PDX & MM468 Persister gene lists.1.png",width=1500,height=1500,res=300)
smoothScatter(Both_genelist$log10_qvalue_MM468,Both_genelist$log10_qvalue_PDX,nbin=600,
              colramp=colorRampPalette(mycolramp), xlab="-log10 qvalue MM468",ylab="-log10 qvalue PDX",
              main=paste0("r=",round(correlation$estimate,2),
                          " pvalue<2.2e-16"),cex=0.8)
    
abline(lmper)
dev.off()

sel <- Both_genelist$log10_qvalue_MM468>7 & Both_genelist$log10_qvalue_PDX>7
text(Both_genelist$log10_qvalue_MM468[sel],Both_genelist$log10_qvalue_PDX[sel],
     labels=Both_genelist$Gene.Set[sel],cex=0.5)

WriteXLS(Both_genelist,ExcelFileName = file.path(outdir,"Persister_PDX_MM468_GeneLists-2.xls"))


#######################
##Comparison Gene FC
#######################

Both_gene <- inner_join(my.res.MM468,my.res.PDX)
correlation <- cor.test(Both_gene$log2FC.C2_pers,Both_gene$log2FC.persister_6_vs_UNT)
png("log2FC PDX & MM468 Persister.png",width=1500,height=1500,res=300)
smoothScatter(Both_gene$log2FC.C2_pers,Both_gene$log2FC.persister_6_vs_UNT,nbin=300,
              colramp=colorRampPalette(mycolramp), xlab="log2FC MM468",ylab="log2FC PDX",
              main=paste0("r=",round(correlation$estimate,2),
                          " pvalue<2.2e-16"),cex=0.8)
dev.off()
WriteXLS(Both_gene, ExcelFileName = "Persister_PDX_MM468_log2FC.xls")


# All 

# PDX
# All 
load(file.path(maindir, "output","scRNAseq", "PDX", "Supervised",
               "RData", "Supervised_res_object_edgeR.Rdata"))

rownames(my.res) <- my.res$Symbol
annotbase <- "MSigDB"
database <- MSIG.ls ##MSigDB

refs <- "UNT"
groups <- "persister_vs_UNT"

All  <- data.frame()
reflist <- unique(my.res$Symbol);length(reflist)

i=1
gp <- groups[i]
if (length(refs)>1) {ref <- refs[i]} else {ref <- refs[1]}
print(paste0("Processing ",gp, " vs ", ref, " _ ",annotbase," for logFC = ",round(log2FC_threshold,2)))

over <- which(my.res[,paste("qval",gp,sep=".")] <= 0.01 & my.res[,paste("log2FC",gp,sep=".")] > log2FC_threshold)

if(length(over)){
    enrich.test <- geco.enrichmentTest(gene.sets=database,mylist=my.res$Symbol[over],possibleIds=reflist)
    enrich.test <- data.frame(Gene_set_name=rownames(enrich.test), enrich.test, check.names=FALSE)
    enrich.test <- merge( subset(MSIG.gs, select=-Genes), enrich.test, by.x="Gene.Set", by.y="Gene_set_name", all.y=TRUE, sort=FALSE ) ## Get class of gene set
    
    enrich.test <- enrich.test[order(enrich.test$`p-value`),]
    #ind <- which(enrich.test$`q-value`<= 0.1);if(!length(ind)){ind <- 1:20}
    #Overexpressed  <- enrich.test[ind,]
    Overexpressed  <- enrich.test
    colnames(  )[c(5,6)] <- c("pvalue","qvalue")
}

Overexpressed_PDX_genesets <- Overexpressed
colnames(Overexpressed_PDX_genesets) <- paste0(colnames(Overexpressed_PDX_genesets),"_PDX")
colnames(Overexpressed_PDX_genesets)[1] <- "Gene.Set"

# MM468 

load(file.path(maindir, "output","scRNAseq", "MM468", "Persister", "Supervised",
               "RData", "Supervised_res_object_edgeR.Rdata"))

rownames(my.res) <- my.res$Gene
annotbase <- "MSigDB" 
database <- MSIG.ls ##MSigDB

refs <- "DMSO"
groups <- "C2_pers"

All  <- data.frame()
reflist <- unique(my.res$Symbol);length(reflist)

i=1
gp <- groups[i]
if (length(refs)>1) {ref <- refs[i]} else {ref <- refs[1]}
print(paste0("Processing ",gp, " vs ", ref, " _ ",annotbase," for logFC = ",round(log2FC_threshold,2)))

over <- which(my.res[,paste("qval",gp,sep=".")] <= 0.01 & my.res[,paste("log2FC",gp,sep=".")] > log2FC_threshold)

if(length(over)){
    enrich.test <- geco.enrichmentTest(gene.sets=database,mylist=my.res$Symbol[over],possibleIds=reflist)
    enrich.test <- data.frame(Gene_set_name=rownames(enrich.test), enrich.test, check.names=FALSE)
    enrich.test <- merge( subset(MSIG.gs, select=-Genes), enrich.test, by.x="Gene.Set", by.y="Gene_set_name", all.y=TRUE, sort=FALSE ) ## Get class of gene set
    
    enrich.test <- enrich.test[order(enrich.test$`p-value`),]
    #ind <- which(enrich.test$`q-value`<= 0.1);if(!length(ind)){ind <- 1:20}
    #Overexpressed  <- enrich.test[ind,]
    Overexpressed  <- enrich.test
    colnames( Overexpressed )[c(5,6)] <- c("pvalue","qvalue")
}

Overexpressed_MM468_genesets <- Overexpressed
colnames(Overexpressed_MM468_genesets) <- paste0(colnames(Overexpressed_MM468_genesets),"_MM468")
colnames(Overexpressed_MM468_genesets)[1] <- "Gene.Set"

All_pathways <- inner_join(Overexpressed_MM468_genesets,Overexpressed_PDX_genesets, by = "Gene.Set")
All_pathways$log10_qvalue_MM468 <- -log10(All_pathways$qvalue_MM468)
All_pathways$log10_qvalue_PDX <- -log10(All_pathways$qvalue_PDX)

correlation <- cor.test(All_pathways$log10_qvalue_MM468,All_pathways$log10_qvalue_PDX) 
lmper <- lm(All_pathways$log10_qvalue_PDX~All_pathways$log10_qvalue_MM468) 
print(correlation$p.value)

png(file.path(outdir,"-log10 q value PDX & MM468 Persister gene lists.1.png"),width=1500,height=1500,res=300)
smoothScatter(All_pathways$log10_qvalue_MM468,All_pathways$log10_qvalue_PDX,nbin=600,
              colramp=colorRampPalette(mycolramp), xlab="-log10 qvalue MM468",ylab="-log10 qvalue PDX",
              main=paste0("r=",round(correlation$estimate,2),
                          " pvalue<2.2e-16"),cex=0.8)

abline(lmper)
dev.off()
