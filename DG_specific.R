####
# Distance between FOS+ and FOS-
####
#DG
samples <- rownames(metaProxC[metaProxC$Subgroup2 == "DG"&
                                metaProxC$Mouse_condition == "EE" & metaProxC$FOS != "L" &metaProxC$outliers == "in",])
dat <- na.exclude(countProxC[, samples])
dat <- dat[rowSums(dat) > 0,]
met <- metaProxC[match(samples,metaProxC$Sample_ID),]
res.DG <- exact(dat, group = (met$FOS == "F"), Pair = c(TRUE,FALSE))
#CA1
samples <- rownames(metaProxC[metaProxC$Subgroup2 == "CA1"&
                                metaProxC$Mouse_condition == "EE" & metaProxC$FOS != "L" &metaProxC$outliers == "in",])
dat <- na.exclude(countProxC[, samples])
dat <- dat[rowSums(dat) > 0,]
met <- metaProxC[match(samples,metaProxC$Sample_ID),]
res.CA1 <- exact(dat, group = (met$FOS == "F"), Pair = c(TRUE,FALSE))
#CA1b
samples <- rownames(metaProxC[metaProxC$Subgroup2 == "CA1b"&
                                metaProxC$Mouse_condition == "EE" & metaProxC$FOS != "L" &metaProxC$outliers == "in",])
dat <- na.exclude(countProxC[, samples])
dat <- dat[rowSums(dat) > 0,]
met <- metaProxC[match(samples,metaProxC$Sample_ID),]
res.CA1b <- exact(dat, group = (met$FOS == "F"), Pair = c(TRUE,FALSE))
#CA3
samples <- rownames(metaProxC[metaProxC$Subgroup2 == "CA3"&
                                metaProxC$Mouse_condition == "EE" & metaProxC$FOS != "L" &metaProxC$outliers == "in",])
dat <- na.exclude(countProxC[, samples])
dat <- dat[rowSums(dat) > 0,]
met <- metaProxC[match(samples,metaProxC$Sample_ID),]
res.CA3 <- exact(dat, group = (met$FOS == "F"), Pair = c(TRUE,FALSE))
#IN
samples <- rownames(metaProxC[metaProxC$Subgroup2 == "IN"&
                                metaProxC$Mouse_condition == "EE" & metaProxC$FOS != "L" &metaProxC$outliers == "in",])
dat <- na.exclude(countProxC[, samples])
dat <- dat[rowSums(dat) > 0,]
met <- metaProxC[match(samples,metaProxC$Sample_ID),]
res.IN <- exact(dat, group = (met$FOS == "F"), Pair = c(TRUE,FALSE))
#VIP
samples <- rownames(metaProxC[metaProxC$Subgroup2 == "VIP"&
                                metaProxC$Mouse_condition == "EE" & metaProxC$FOS != "L" &metaProxC$outliers == "in",])
dat <- na.exclude(countProxC[, samples])
dat <- dat[rowSums(dat) > 0,]
met <- metaProxC[match(samples,metaProxC$Sample_ID),]
res.VIP <- exact(dat, group = (met$FOS == "F"), Pair = c(TRUE,FALSE))


##########
#
df <- data.frame(count = c(sum(res.DG$logFC > 0 & res.DG$f < 0.05),
                           sum(res.CA1$logFC > 0 & res.CA1$f < 0.05),
                           sum(res.CA1b$logFC > 0 & res.CA1b$f < 0.05),
                           sum(res.CA3$logFC > 0 & res.CA3$f < 0.05),
                           sum(res.IN$logFC > 0 & res.IN$f < 0.05),
                           sum(res.VIP$logFC > 0 & res.VIP$f < 0.05),
                           
                           sum(res.DG$logFC < 0 & res.DG$f < 0.05),
                           sum(res.CA1$logFC < 0 & res.CA1$f < 0.05),
                           sum(res.CA1b$logFC < 0 & res.CA1b$f < 0.05),
                           sum(res.CA3$logFC < 0 & res.CA3$f < 0.05),
                           sum(res.IN$logFC < 0 & res.IN$f < 0.05),
                           sum(res.VIP$logFC < 0 & res.VIP$f < 0.05)),
              celltype = c("DG","CA1","CA1b","CA3","IN","VIP"),
              FOS = c(rep(c("N","F"),each = 6))
)


df$FOS <- factor(df$FOS, levels = c("N","F"))
ggplot(df[df$celltype != "CA1b" & df$celltype!= "CA3",], aes(celltype,count, fill = FOS))+
  geom_bar(stat = 'identity', position = "dodge")+
  scale_fill_manual(values = c("blue","red"))+
  ylab("DifExp Genes, count")+
  labs(title = "Differential Expression")+
  theme_bw()+
  theme(text=element_text(size=20))+
  theme(panel.border = element_rect(colour=c("black"),size=2),
        axis.ticks = element_line(size=1.5))


###########
#
dg.genes <- rownames(res.DG[res.DG$logFC < 0 & res.DG$f < 0.05 ,])
ca1.genes <- rownames(res.CA1[res.CA1$logFC < 0 & res.CA1$f < 0.05,])
in.genes <- rownames(res.IN[res.IN$logFC < 0 & res.IN$f < 0.05,])
vip.genes <- rownames(res.VIP[res.VIP$logFC < 0 & res.VIP$f < 0.05,])


a <- table(c(rep(dg.genes,5),
       ca1.genes,
       in.genes,
       vip.genes))

genes <- names(a[a==5])
b <- res.DG[genes,]
b <- b[order(b$f),]
######
# Difexp versus 
######
samples <- rownames(metaProxC[metaProxC$Mouse_condition == "EE" &  metaProxC$Subgroup2 != "DG" & metaProxC$FOS == "N" & metaProxC$Subgroup2 != "HDG" & metaProxC$outliers == "in"|
                                metaProxC$Mouse_condition == "EE" &  metaProxC$Subgroup2 == "DG" & metaProxC$FOS == "F" & metaProxC$outliers == "in",])
dat <- na.exclude(countProxC[, samples])
dat <- dat[rowSums(dat) > 0,]
met <- metaProxC[match(samples,metaProxC$Sample_ID),]
res.DGonly <- exact(dat, group = (met$FOS == "F"), Pair = c(TRUE,FALSE))

res2 <- res.DGonly[genes,]
dg.iegspecific <- rownames(na.exclude(res2[res2$logFC < 0 & res2$f < 0.05,]))
dg.quiet <- rownames(na.exclude(res2[res2$f > 0.05,]))

###### 
# Identify genes that are already on in other cell types
###### 
samples <- rownames(metaProxC[metaProxC$Mouse_condition == "EE" &  metaProxC$Subgroup2 != "DG" & metaProxC$FOS == "N" & metaProxC$Subgroup2 != "HDG" & metaProxC$outliers == "in"|
 metaProxC$Mouse_condition == "EE" &  metaProxC$Subgroup2 == "DG" & metaProxC$FOS == "F" & metaProxC$outliers == "in",])
dat <- na.exclude(tpmProxC[, samples])
dat <- dat[rowSums(dat) > 0,]
met <- metaProxC[match(samples,metaProxC$Sample_ID),]


prop <- vector()
for (i in unique(met$Subgroup2)){
  samples2 <- rownames(met[met$Subgroup2 == i,])
  tmp2 <- apply(X = dat[,samples2], MARGIN = 1, FUN = propExp)
  prop <- cbind(prop, as.numeric(tmp2))
}
colnames(prop) <- unique(met$Subgroup2)
rownames(prop) <- rownames(dat)

#
prop2 <- prop[genes,]
prop3 <- melt(t(prop2[,-c(3)]))
prop3$DG <- rep(prop2[,3],each = 5)

ggplot(prop3[prop3$X1 !="CA3",], aes(DG, value, colour = X1))+
  geom_point()

##
dat2 <- as.matrix(dat[genes,])
colnames(dat2) <- paste(met$Subgroup2, met$FOS)
heatmap(dat2)
