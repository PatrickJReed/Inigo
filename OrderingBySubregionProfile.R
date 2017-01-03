######
samples <- metaProxC[metaProxC$Mouse_condition == "HC" & metaProxC$Brain_Region == "CA3_other_negs" & metaProxC$alignable >  500000 & metaProxC$Smartseq2_RT_enzyme_used == "ProtoscriptII" ,"Sample_ID"]#
#metaProxC$CTIP2 == "N" & metaProxC$PROX1 == "N" & metaProxC$FOS == "N" & metaProxC$Mouse_condition == "HC" & metaProxC$alignable >  500000 & metaProxC$Smartseq2_RT_enzyme_used == "ProtoscriptII"  ,"Sample_ID"]
dat <- tpmProxC[, samples]
met <- metaProxC[samples,]
colnames(dat) <- paste(met$Brain_Region, met$FOS,c(1:ncol(dat)),sep = ".")
genes <- c("Pcp4","Reln","Sst","Bcl2","Rfx3","Man1a2","Sema3e","Dusp6","Thsd4","Mpped1","Map3k15","Cdh24")
dat <- dat[genes,]
for (g in genes){
  dat[g, ] <- scale(as.numeric(dat[g,]))[,1]
}
tmp <- melt(t(dat))
### Get "phenotype"
a <- vector()
for (g in genes){
  a <- rbind(a,as.numeric(rank(tpmProxC[g,samples],ties.method = "average" )))
}

b <- vector()
for (i in 1:ncol(a)){
  b <- c(b, which(a[,i] == max(a[,i]))[1])
}

### Order samples and genes
tmp$geneorder <-rep(1:length(genes),each = ncol(dat))
h <- cutree(hclust((dist(t(dat)))),k = 10)
tmp$sampleorder <- rank(as.numeric(b),ties.method = "first")

#plot

ggplot(tmp, aes(reorder(X1, sampleorder), (X2), fill = value)) + 
  geom_tile()+
  scale_fill_gradient(high = "red",low= "blue")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

