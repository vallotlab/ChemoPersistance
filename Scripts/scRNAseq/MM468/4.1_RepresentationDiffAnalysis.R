#analysis fold-change MM468


load("Supervised_res_Wilcox_object.Rdata")


#Compare diff log2FC across conditions


my.res$log2FC_expressing_5FU6 <- my.res$Exp.level.MM468_5FU6_day33-my.res$Exp.level.DMSOsamp
my.res$log2FC_expressing_UNC <- my.res$Exp.level.MM468_UNC_day33-my.res$Exp.level.DMSOsamp
my.res$log2FC_number_exp_UNC <- log(my.res$Fraction.exp.MM468_UNC_day33/my.res$Fraction.exp.DMSO,2)
my.res$log2FC_number_exp_5FU6 <- log(my.res$Fraction.exp.MM468_5FU6_day33/my.res$Fraction.exp.DMSO,2)

mycolramp <- c("white",viridis(n=4))

test <- cor.test(my.res$log2FC.MM468_UNC_day33, my.res$log2FC.MM468_5FU6_day33)
png()
smoothScatter(my.res$log2FC.MM468_UNC_day33, my.res$log2FC.MM468_5FU6_day33,cex=1, 
              nbin=300,colramp=colorRampPalette(mycolramp),xlim=c(-5,5),ylim=c(-5,5),xlab="log2FC UNC vs initial",
              ylab="log2FC 5FU6 persister vs initial",bandwidth=c(0.01,0.01)) 

smoothScatter(my.res$log2FC.MM468_UNC_day33, my.res$log2FC.MM468_5FU6_day214,cex=1, 
              colramp=colorRampPalette(mycolramp),xlab="log2FC UNC vs initial",nbin=300,
              ylab="log2FC 5FU6 resistant vs initial",bandwidth=c(0.01,0.01),xlim=c(-5,5),ylim=c(-5,5)) 

smoothScatter(my.res$log2FC.MM468_UNC_day33, my.res$log2FC.MM468_5FU3_day202,cex=1, 
              colramp=colorRampPalette(mycolramp),xlab="log2FC UNC vs initial",nbin=300,
              ylab="log2FC 5FU3 resistant vs initial",bandwidth=c(0.01,0.01),xlim=c(-5,5),ylim=c(-5,5)) 

smoothScatter(my.res$log2FC.MM468_UNC_day33, my.res$log2FC.MM468_5FU5_day171,cex=1, 
              colramp=colorRampPalette(mycolramp),xlab="log2FC UNC vs initial",nbin=300,
              ylab="log2FC 5FU5 resistant vs initial",bandwidth=c(0.01,0.01),xlim=c(-5,5),ylim=c(-5,5)) 

sel <- my.res$Symbol %in% c("KRT14","KRT17","KRT16","TGFB1","TGFBR1")
text(my.res$log2FC.MM468_UNC_day33[sel], my.res$log2FC.MM468_5FU3_day202[sel],labels=my.res$Symbol[sel])


abline(v=1.5)
abline(h=1.5)



smoothScatter(my.res$log2FC_number_exp_UNC, my.res$log2FC_number_exp_5FU6,cex=1, 
              colramp=colorRampPalette(mycolramp),xlab="log2FC UNC vs initial",nbin=300,
              ylab="log2FC 5FU6 resistant vs initial",bandwidth=c(0.01,0.01),xlim=c(-5,5),ylim=c(-5,5)) 

sel <- my.res$Symbol %in% c("KRT14","KRT17","KRT16","TGFB1","TGFBR1")
sel <- my.res$log2FC.MM468_UNC_day33>3 & my.res$log2FC.MM468_5FU6_day33>3
text(my.res$log2FC_number_exp_UNC[sel], my.res$log2FC_number_exp_5FU6[sel],labels=my.res$Symbol[sel])


smoothScatter(my.res$log2FC.MM468_UNC_day33, my.res$log2FC.MM468_5FU6_day214,cex=1, 
              nbin=300,colramp=colorRampPalette(mycolramp),xlim=c(-4,4),ylim=c(-4,4),xlab="log2FC UNC vs initial",
              ylab="log2FC 5FU6 persister vs initial")
abline(v=1)
abline(h=1)

smoothScatter(my.res$log2FC.MM468_UNC_day33, my.res$log2FC.MM468_DMSO3_day50,cex=1, 
              nbin=300,colramp=colorRampPalette(mycolramp),xlim=c(-4,4),ylim=c(-4,4),xlab="log2FC UNC vs initial",
              ylab="log2FC 5FU6 persister vs initial")
abline(v=1)
abline(h=1)


ggplot(my.res, aes(log2FC.MM468_UNC_day33, log2FC.MM468_5FU6_day214)) + 
  geom_point(shape=16, size=0.25, show.legend = FALSE) +
  stat_bkde2d(aes(fill=..level..), geom="polygon") +
  scale_fill_viridis() + theme_bw()

ggplot(my.res, aes(log2FC.MM468_UNC_day33, log2FC.MM468_5FU6_day214)) + 
geom_point(shape=16, size=0.25, show.legend = FALSE, color="red") +
  scale_fill_viridis()

ggplot(my.res, aes(log2FC.MM468_UNC_day33, log2FC.MM468_5FU6_day33)) + 
  geom_point(shape=16, size=0.25) +
  stat_bkde2d(bandwidth=c(0.15, 0.15), geom="polygon",aes(fill=..level..)) +
  
  scale_fill_viridis() +
  theme_bw() +
  theme(panel.grid=element_blank())


my.res$dens_col <- densCols(my.res$log2FC.MM468_UNC_day33, my.res$log2FC.MM468_5FU6_day33, nbin=80,colramp=viridis_pal())
ggplot() + geom_point(data=my.res, 
           aes(x=log2FC.MM468_UNC_day33, y=log2FC.MM468_5FU6_day33, color=dens_col),
           size=0.6, alpha=1/4) + scale_color_identity()

mycolramp <- c("white",viridis(n=5))
smoothScatter(my.res$log2FC.MM468_UNC_day33, my.res$log2FC.MM468_5FU6_day33,cex=1, 
              nbin=300,colramp=colorRampPalette(mycolramp),xlim=c(-4,4),ylim=c(-4,4))

smoothScatter(my.res$Fraction.exp.MM468_UNC_day33, my.res$Exp.level.MM468_UNC_day33,cex=1, 
              nbin=300,colramp=colorRampPalette(mycolramp))

smoothScatter(my.res$Fraction.exp.MM468_DMSO3_day50, my.res$Exp.level.MM468_DMSO3_day50,cex=1, 
              nbin=300,colramp=colorRampPalette(mycolramp),ylim=c(0,10))


smoothScatter(my.res$Fraction.exp.MM468_5FU6_day33, my.res$Exp.level.MM468_5FU6_day33,cex=1, 
              nbin=300,colramp=colorRampPalette(mycolramp), ylim=c(0,10))             

