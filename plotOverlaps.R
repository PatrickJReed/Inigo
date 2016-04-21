#Celltype overlap
overlap <- read.table(as.matrix("~/Documents/SalkProjects/ME/ShortLongSingature/raw/celltype_overlap_melt.txt"),header = TRUE)

ggplot(overlap[overlap$X1 != "IN" & overlap$X1 != "CA23" & overlap$X2 != "IN" & overlap$X2 != "CA23",], aes(reorder(X1,X1order), reorder(X2,x2order), fill = value))+
  geom_tile()+
  scale_fill_gradient(high = "red", low = "blue")+
  xlab("Cell Type with Higher Expression")+
  ylab("Cell Type with Lower Expression")+
  labs(title=c("Differential Expression Between Cell Types"))



#Activity overlap
overlap <- read.table(as.matrix("~/Documents/SalkProjects/ME/ShortLongSingature/raw/celltypeactivity_overlap_melt.txt"),header = TRUE)

ggplot(overlap, aes(reorder(X2,X2order), reorder(X1,X1order), fill = value))+
  geom_tile()+
  scale_fill_gradient(high = "red", low = "blue")+
  xlab("FOS")+
  ylab("Cell Type")+
  labs(title=c("Differential Expression Between FOS States"))






samples <- metaProxC[metaProxC$Mouse_condition == 'EE' &   metaProxC$FOS == "F" &metaProxC$alignable >  500000 & metaProxC$Smartseq2_RT_enzyme_used == "ProtoscriptII" ,"Sample_ID"]#|
# metaProxC$subgroup == "CA3" & metaProxC$Mouse_condition == "HC" & metaProxC$alignable >  500000 & metaProxC$Smartseq2_RT_enzyme_used == "ProtoscriptII","Sample_ID"]
dat <- tpmProxC[, samples]
met <- metaProxC[match(samples,metaProxC$Sample_ID),]
colnames(dat) <- met$Brain_Region
COR <- cor(na.exclude(dat))
h <- heatmap(COR)
a <- COR[,73:93]
b <- rowMeans(a)
a[abs(a) > 0.5]

nm <- c("X151221_34","X151221_39")
metaProxC[nm,]
tpmProxC["Dcx",nm]




