overlap <- read.table(as.matrix("~/Documents/SalkProjects/ME/ShortLongSingature/raw/celltype_overlap_melt.txt"),header = TRUE)

ggplot(overlap[overlap$X1 != "IN" & overlap$X1 != "CA23" & overlap$X2 != "IN" & overlap$X2 != "CA23",], aes(reorder(X1,X1order), reorder(X2,x2order), fill = value))+
  geom_tile()+
  scale_fill_gradient(high = "red", low = "blue")+
  xlab("Cell Type with Higher Expression")+
  ylab("Cell Type with Lower Expression")+
  labs(title=c("Differential Expression Between Cell Types"))


samples <- metaProxC[metaProxC$Mouse_condition == "HC" & metaProxC$FOS == "N" & metaProxC$alignable >  500000 & metaProxC$Smartseq2_RT_enzyme_used == "ProtoscriptII" ,"Sample_ID"]#
#metaProxC$CTIP2 == "N" & metaProxC$PROX1 == "N" & metaProxC$FOS == "N" & metaProxC$Mouse_condition == "HC" & metaProxC$alignable >  500000 & metaProxC$Smartseq2_RT_enzyme_used == "ProtoscriptII"  ,"Sample_ID"]
tmp <- dat <- tpmProxC[, samples]
met <- metaProxC[samples,]

colnames(dat) <- met$Brain_Region

cor <- as.matrix(cor(na.exclude(dat)))
heatmap(cor,scale = "col")
