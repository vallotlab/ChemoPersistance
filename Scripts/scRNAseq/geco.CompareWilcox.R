geco.CompareWilcox <- function (dataMat = NULL, annot = NULL, ref = NULL, groups = NULL, 
          featureTab = NULL, logvalues = FALSE) 
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
    mat. <- dataMat[, c(as.character(refsamp), as.character(gpsamp))]
    # function(x) wilcox.test(as.numeric(x[as.character(refsamp)]), 
    #                         as.numeric(x[as.character(gpsamp)]))
    testWilc <- apply(dataMat, 1, function(x) wilcox.test(as.numeric(x[as.character(refsamp)]), 
                                                          as.numeric(x[as.character(gpsamp)])))
    pval.gpsamp <- unlist(lapply(testWilc, function(x) x$p.value))
    qval.gpsamp <- p.adjust(pval.gpsamp, method = "BH")
    Count.gpsamp <- apply(dataMat, 1, function(x) mean(x[as.character(gpsamp)]))
    
    
  
   
    
    Numgp <- length(gpsamp)
    Numref <- length(refsamp)
    
    Fraction.exp.gpsamp <- apply(dataMat,1,function(x) length(which(x>0 & names(x) %in% as.character(gpsamp) ))/Numgp)
    Fraction.exp.refsamp <- apply(dataMat,1,function(x) length(which(x>0 & names(x) %in% as.character(refsamp) ))/Numref)
    Exp.level.gpsamp <- apply(dataMat, 1, function(x) mean(x[which(x>0 & names(x) %in% as.character(gpsamp) )]))
    Exp.level.refsamp <- apply(dataMat, 1, function(x) mean(x[which(x>0 & names(x) %in% as.character(refsamp) )]))
    
    #Compare level of expression in expressing cells only
    
    testWilc_exp <- apply(dataMat[which(Fraction.exp.gpsamp>0 & Fraction.exp.refsamp >0),], 1, function(x) wilcox.test(as.numeric(x[which(x>0 & names(x) %in% as.character(gpsamp) )]), 
                                                          as.numeric(x[which(x>0 & names(x) %in% as.character(refsamp) )])))
    pval.exp.gpsamp <- unlist(lapply(testWilc_exp, function(x) x$p.value))
    qval.exp.gpsamp <- p.adjust(pval.exp.gpsamp, method = "BH")
    
    if (logvalues) {
      log2FC.gpsamp <- apply(dataMat, 1, function(x) (mean(x[as.character(gpsamp)]) - 
                                                        mean(x[as.character(refsamp)])))
    }
    else {
      log2FC.gpsamp <- apply(dataMat, 1, function(x) log(mean(x[as.character(gpsamp)])/mean(x[as.character(refsamp)]), 
                                                         2))
    }
    res <- data.frame(res, Count.gpsamp, log2FC.gpsamp, pval.gpsamp, 
                      qval.gpsamp, Fraction.exp.gpsamp,Fraction.exp.refsamp,Exp.level.gpsamp,Exp.level.refsamp)
    
    res$pval.exp.gpsamp <- NA
      res$qval.exp.gpsamp <- NA
      res$pval.exp.gpsamp[which(Fraction.exp.gpsamp>0 & Fraction.exp.refsamp >0)]  <-  pval.exp.gpsamp
      res$qval.exp.gpsamp[which(Fraction.exp.gpsamp>0 & Fraction.exp.refsamp >0)]  <-  qval.exp.gpsamp
      
    colnames(res) <- sub("ref", names(ref)[min(c(k, length(ref)))], 
                         sub("gpsamp", names(groups)[k], colnames(res)))
  }
  res
}
