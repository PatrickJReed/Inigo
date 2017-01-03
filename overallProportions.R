## All samples
p <- apply(X = tpmProxC, MARGIN = 1, FUN = propExp)
m <- apply(X = tpmProxC, MARGIN = 1, FUN = meanNoZero)
tmp <- data.frame(detection = p, mean = m)

#Each subgroup HC
for (y in levels(as.factor(metaProxC$Mouse_condition))){
  for (f in levels(as.factor(metaProxC$FOS))){
  for (x in levels(as.factor(metaProxC$Brain_Region))){
    samples <- metaProxC[metaProxC$Brain_Region == x &  metaProxC$FOS == f & metaProxC$Mouse_condition == y & metaProxC$alignable >  500000 & metaProxC$Smartseq2_RT_enzyme_used == "ProtoscriptII" ,"Sample_ID"]#|
    dat <- tpmProxC[, samples]
    nm1 <- paste(x,y,f,"p", sep = ".")
    tmp[,nm1] <- apply(X = dat, MARGIN = 1, FUN = propExp)
    nm2 <- paste(x,y,f,"m", sep = ".")
    tmp[,nm2] <- apply(X = dat, MARGIN = 1, FUN = meanNoZero)
    nm3 <- paste(x,y,f,"100", sep = ".")
    tmp[,nm3] <- tmp[,nm1] > 0.9
  }
}
}
tmp$color.group <- as.factor(apply(tmp[,c(grep("100",colnames(tmp)))],MARGIN = 1,FUN = sum))


ggplot(tmp, aes(mean,detection,colour = color.group))+
  geom_point(alpha =0.1)+
  theme_bw()+
  labs(title = "Gene Detection")+
  theme(text = element_text(size=18),panel.border = element_rect(colour=c("black"),size=2),
        axis.ticks = element_line(size=1.5),legend.position = "none")


samples <- metaProxC[  metaProxC$Smartseq2_RT_enzyme_used == "ProtoscriptII" ,"Sample_ID"]#|

met <- metaProxC[samples,]
tmp2 <- melt(t(met[,c("totalreads", "alignable")]))
tmp2$filtered <- tmp2$value < 500000

ggplot(tmp2[tmp2$X1 == "alignable",],aes(reorder(X2,value),value,color = filtered))+
  geom_bar(stat = 'identity')+
  ylab("Aligned Reads")+
  labs(title = "Sample Filtering")

hist(metaProxC$alignable / metaProxC$totalreads, xlab=c("Percent Reads Aligned"), main = "Single-nuclei Alignment")
abline(v =  mean(na.exclude(metaProxC$alignable / metaProxC$totalreads)),col=c("Red"))

ggplot(metaProxC, aes(Brain_Region, alignable/totalreads, fill = FOS))+
  geom_boxplot()+
  theme_bw()+
  labs(title = "Alignment Based on Condition")+
  theme(text = element_text(size=18),panel.border = element_rect(colour=c("black"),size=2),
        axis.ticks = element_line(size=1.5),legend.position = "none")

############################
##
genes <- rownames(dg[dg[,"f"] < 0.05,] )
tmp2 <- tmp[genes,]
min <- 7
dg.genes2 <- rownames(tmp2[tmp2$DG.EE.F.m > 10 &  tmp2$DG.EE.N.m < min & tmp2$CA1.EE.F.m <= min & tmp2$CA3_other_negs.EE.F.m <= min & 
       tmp2$HDG.EE.F.m <= min ,])
dg.only <- RES[[9]][genes2,]
dg.only <- dg.only[order(dg.only$PValue),]
