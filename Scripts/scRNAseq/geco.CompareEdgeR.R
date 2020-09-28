geco.CompareedgeRRLT <- function (dataMat = NULL, annot = NULL, ref = NULL, groups = NULL, 
          featureTab = NULL, norm_method = "RLE") 
{
  res <- featureTab
  for (k in 1:length(groups)) {
    print(paste("Comparing", names(ref)[min(c(k, length(ref)))], 
                "versus", names(groups)[k]))
    if (length(ref) == 1) {
      refsamp <- ref[[1]]
    }
    else {
      refsamp <- ref[[k]]
    }
    gpsamp <- groups[[k]]
    annot. <- annot[c(refsamp, gpsamp), 1:2]
    annot.$Condition <- c(rep("ref", length(refsamp)), rep("gpsamp", 
                                                           length(gpsamp)))
    mat. <- as.data.frame(dataMat)[, c(refsamp,gpsamp)]
    edgeRgroup <- c(rep(1, length(refsamp)), rep(2, length(gpsamp)))
    y <- DGEList(counts = mat., group = edgeRgroup)
    y <- calcNormFactors(y, method = norm_method)
  
      
    NormCounts <- 100000*t(t(mat.) /  (y$samples$lib.size*y$samples$norm.factors))
    LogdataMat <- log(NormCounts+1,2)
    rm(NormCounts)
    
    design <- model.matrix(~edgeRgroup)
    y <- estimateDisp(y, design)
    fit <- glmFit(y, design)
    comp <- glmLRT(fit)
    pvals <- comp$table
    colnames(pvals) <- c("log2FC.gpsamp", "logCPM.gpsamp", 
                         "LR.gpsamp", "pval.gpsamp")
    pvals$qval.gpsamp <- p.adjust(pvals$pval.gpsamp, method = "BH")
    
    
    
    Fraction.exp.gpsamp <- apply(dataMat,1,function(x) length(which(x>0 & names(x) %in% as.character(gpsamp) )))
    Fraction.exp.ref <- apply(dataMat,1,function(x) length(which(x>0 & names(x) %in% as.character(refsamp) )))
    
    #Problem to calculate logCPM in expressing cells, two way selection on dataMat and NormCounts
    # NormCounts <- edgeR::cpm(y, log=TRUE)
    # logCPM.exp.gpsamp <-  apply(dataMat,1,function(row,df2)
    #                          {
    #                           apply(df2,1,function(row1,row2) {mean(row1[names(which(row2>0 & names(row2) %in% as.character(gpsamp) ))]) },row)
    #                         },NormCounts)
    # 
    #sinon use LogCounts 
    logCount.exp.gpsamp <- apply(LogdataMat,1,function(x) mean(x[which(x>0 & names(x) %in% as.character(gpsamp) )]))
    logCount.exp.ref <- apply(LogdataMat,1,function(x) mean(x[which(x>0 & names(x) %in% as.character(refsamp) )]))
    
    
    
    res <- data.frame(res, pvals[, c("log2FC.gpsamp", "logCPM.gpsamp", 
                                     "pval.gpsamp", "qval.gpsamp")],Fraction.exp.gpsamp,Fraction.exp.ref,logCount.exp.gpsamp,logCount.exp.ref)
    colnames(res) <- sub("ref", names(ref)[min(c(k, length(ref)))], 
                         sub("gpsamp", names(groups)[k], colnames(res)))
  }
  res
}


