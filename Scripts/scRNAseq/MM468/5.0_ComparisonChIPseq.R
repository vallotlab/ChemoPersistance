QCdir <- "~/Desktop/scRNAseq_data_local/QC/"
resdir <- "~/Desktop/scRNAseq_data_local/Results_MM468/"
resSUBdir <- file.path(resdir, paste0("ComparisonChIPseq")) ;if(!file.exists(resSUBdir)){dir.create(resSUBdir)}


#Comparing scRNAseq datasets and ChIPseq datasets

load("~/Desktop/scRNAseq_data_local/Results_MM468/Supervised/RData/Supervised_res_object_edgeR.Rdata")
res.scRNA <- my.res

load("~/Google Drive/MM468_bulk/ChIP/results/Supervised/RData/Supervised_res_Limma_object.Rdata")
res.ChIP <- res
res.ChIPgene <- res.ChIP[res.ChIP$peak_affectation %in% c("tss_pc"),]
sig.ChIP <- res.ChIPgene[abs(res.ChIPgene$log2FC.X5FU2_3_5)>1.5 & res.ChIPgene$qval.X5FU2_3_5<0.1,]

trial <- function(x) {
  indice <- grep(pattern = x, res.ChIPgene$gene_affectation)
  res.ChIPgene$ID[indice]
}

test <- sapply(res.scRNA$Symbol,FUN=trial)

res.scRNA$ChIP <- 0
for(i in 1:length(res.scRNA$Symbol)){
  if(length(test[[i]])!=0) res.scRNA$ChIP[i] <- c(test[[i]])
}

#Within expressed genes in the model system, n=14924 genes, 2083 have a H3K27me3 peak over their tss or gene_body
res.scRNA$peak_affectation <- res.ChIP$peak_affectation[match(res.scRNA$ChIP,res.ChIP$ID)]
res.scRNA$log2FC_ChIP <- res.ChIP$log2FC.X5FU2_3_5[match(res.scRNA$ChIP,res.ChIP$ID)]
res.scRNA$qvalue_ChIP <- res.ChIP$qval.X5FU2_3_5[match(res.scRNA$ChIP,res.ChIP$ID)]
res.scRNA$bivalent <- res.ChIP$bivalent[match(res.scRNA$ChIP,res.ChIP$ID)]

#sig.scRNA <- res.scRNA[abs(res.scRNA$log2FC.persisterall)>1 & res.scRNA$qval.persisterall<0.01,]

sig.peaks <- res.ChIPgene$ID[abs(res.ChIPgene$log2FC.X5FU2_3_5)>1 & res.ChIPgene$qval.X5FU2_3_5<0.1]

#Study the correlation between expression and changes in H3K27me3 enrichment on peaks that are diff enriched
subset_int <- res.scRNA[res.scRNA$ChIP %in% sig.peaks,]

top_demethylated <- head(subset_int$Symbol[order(subset_int$log2FC.persisterall,decreasing=T)],n=15)

pdf(file.path(resSUBdir,"Log2FCscRNAseq_for_signigicantPeaks.pdf"),height=5,width=5)
sp <- ggscatter(subset_int, y = "log2FC.persisterall", x = "log2FC_ChIP",
                add = "reg.line",  # Add regressin line
                add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE,
                label="Symbol",repel=TRUE,
                label.select= top_demethylated) #add size of dot function of initial expression in DMSO
# Add correlation coefficient
sp + stat_cor(method = "pearson",label.x=1,label.sep = "\n")

sp <- ggscatter(subset_int, y = "log2FC.persisterall", x = "log2FC_ChIP",
                add = "reg.line",  # Add regressin line
                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE,
                label="Symbol",repel=TRUE,
                label.select= top_demethylated) #add size of dot function of initial expression in DMSO
# Add correlation coefficient
sp + stat_cor(method = "pearson",label.x=1,label.sep = "\n")
sp <- ggscatter(subset_int, y = "log2FC.persisterall", x = "log2FC_ChIP",
                add = "reg.line",  # Add regressin line
                add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE) #add size of dot function of initial expression in DMSO
# Add correlation coefficient
sp + stat_cor(method = "pearson",label.x=1,label.sep = "\n")
dev.off()

# sp + geom_density_2d()
# 
# sp + stat_density_2d(aes(fill = ..level..), geom = "polygon")
# 



#Evaluate how much of diff expression is linked to H3K27me3

library(dplyr)
res.scRNA$decile <- ntile(res.scRNA$log2FC.persisterall, 10)  
figure <- res.scRNA %>% group_by(decile) %>% summarise(ChIP=length(which(log2FC_ChIP<(-2) & qvalue_ChIP<0.1)),ChIP_res=length(which(log2FC_ChIP<(-2) & qvalue_ChIP<0.1)), expression=max(log2FC.persisterall))


pdf(file.path(resSUBdir,"ExpressionDecileforSignificantDepletedPeaks.pdf"),height=5,width=5)

ggplot(data=figure, aes(x=decile , y=ChIP, fill=decile)) +
  geom_rect(data=figure,aes(xmin=decile-0.5,xmax=decile+0.5,ymin=-Inf,ymax=Inf),fill=alpha(colorRampPalette(c("royalblue","white","indianred1"))(10),0.8))+ 
  geom_point(data=figure, aes(x=decile, y=ChIP)) + theme_classic() + theme(legend.position = "none") + xlab("log2FC persister vs DMSO decile") + ylab("number of significantly depleted peaks")  + geom_text(aes(x=decile,y=40,label=round(expression,2)))

dev.off()


