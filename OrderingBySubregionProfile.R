######
samples <- metaProxC[metaProxC$Mouse_condition == "HC" & metaProxC$Brain_Region == "CA3_other_negs" & metaProxC$alignable >  500000 & metaProxC$Smartseq2_RT_enzyme_used == "ProtoscriptII" ,"Sample_ID"]#
#metaProxC$CTIP2 == "N" & metaProxC$PROX1 == "N" & metaProxC$FOS == "N" & metaProxC$Mouse_condition == "HC" & metaProxC$alignable >  500000 & metaProxC$Smartseq2_RT_enzyme_used == "ProtoscriptII"  ,"Sample_ID"]
dat <- tpmProxC[, samples]
met <- metaProxC[samples,]
colnames(dat) <- paste(met$Brain_Region, met$FOS,c(1:ncol(dat)),sep = ".")
genes <- c("Pcp4","Reln","Bcl2","Rfx3","Man1a2","Sema3e")
dat <- dat[genes,]
for (g in genes){
  dat[g, ] <- scale(as.numeric(dat[g,]))[,1]
}
tmp <- melt(t(dat))
### Get "phenotype"
samples <- samples
a.1 <- rank(tpmProxC["Pcp4",samples],ties.method = "average" )
a.2 <- rank(tpmProxC["Reln",samples],ties.method = "average" )
a.3 <- rank(tpmProxC["Sst",samples],ties.method = "average" )
a.4 <- rank(tpmProxC["Bcl2",samples],ties.method = "average" )
a.5 <- rank(tpmProxC["Rfx3",samples],ties.method = "average" )
a.6 <- rank(tpmProxC["Man1a2",samples],ties.method = "average" )
a.7 <- rank(tpmProxC["Sema3e",samples],ties.method = "average" )

a <- rbind(a.1,a.2,a.3,a.4,a.5,a.6,a.7)
b <- vector()
for (i in 1:ncol(a)){
  b <- c(b, which(a[,i] == max(a[,i]))[1])
}

### Order samples and genes
tmp$geneorder <- rep(c(1:length(genes)),each = ncol(dat))
h <- cutree(hclust((dist(t(dat)))),k = 10)
tmp$sampleorder <- rank(as.numeric(b),ties.method = "first")

#plot
ggplot(tmp, aes(reorder(X1, sampleorder), reorder(X2,geneorder), fill = value)) + 
  geom_tile()+
  scale_fill_gradient(high = "red",low= "blue")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

