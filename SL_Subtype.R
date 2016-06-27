library(cluster)
####################
# Identifying Subgroups by Multivariate Clustering
####################
#Step1 : All FOS- Homecage
####################
# Saved metaProxC with the Subtypes labels here as backup
#write.table(x = c("metaProxC"),file = "~/Documents/SalkProjects/ME/ShortLongSingature/SLSig_resTables/metaProxC_062716.txt",quote = FALSE)
####

a <- t.hc
for(g in genes){
  a[,g] <- NA
  a[,g] <- as.numeric(tpmProxC[g,rownames(t.hc)])
}
rownames(a) <- paste(a$Brain_Region, c(1:nrow(a)),sep = "")
plot(hclust(dist(a[-c(25),c(1:3,19:20,39:73)])))
rownames(a) <- rownames(t.hc)
k <- cutree(hclust(daisy(a[-c(25),c(1:3,19:20,39:73)])),k = 5)
metaProxC$Subgroup <- as.character(metaProxC$Brain_Region)
metaProxC[names(k),"Subgroup"] <- as.numeric(k)

table(metaProxC[names(k),"Brain_Region"], k)

metaProxC[metaProxC$Subgroup == 1,"Subgroup"] <- "CA1"
metaProxC[metaProxC$Subgroup == 2,"Subgroup"] <- "Negs"
metaProxC[metaProxC$Subgroup == 3,"Subgroup"] <- "Negs.2"
metaProxC[metaProxC$Subgroup == 4,"Subgroup"] <- "Vip"
metaProxC[metaProxC$Subgroup == 5,"Subgroup"] <- "DG"

####################
#Step2 : CA1
####################
rownames(a) <- paste(a$Brain_Region, c(1:nrow(a)),sep = "")
a$Subgroup <- metaProxC[rownames(t.hc),"Subgroup"]
plot(hclust(daisy(a[a$Subgroup == "CA1",c(1:3,19:20,39:73)])))
k <- cutree(hclust(daisy(a[a$Subgroup == "CA1",c(1:3,19:20,39:73)])),k = 3)

metaProxC[names(k[k==1]),"Subgroup"] <- "CA1.1"
metaProxC[names(k[k==2]),"Subgroup"] <- "CA1.2"

####################
#Step2 : Negs Major Group
####################
#rownames(a) <- paste(a$Brain_Region, c(1:nrow(a)),sep = "")
a$Subgroup <- metaProxC[rownames(t.hc),"Subgroup"]
plot(hclust(daisy(a[a$Subgroup == "Negs",c(1:3,19:20,39:73)])))
rownames(a) <- rownames(t.hc)
k <- cutree(hclust(daisy(a[a$Subgroup == "Negs",c(1:3,19:20,39:73)])),k = 2)

metaProxC[names(k[k==1]),"Subgroup"] <- "Negs.1.1"
metaProxC[names(k[k==2]),"Subgroup"] <- "Negs.1.2"


#########################
## Negs
#########################
samples <- metaProxC[ metaProxC$Mouse_condition == "HC" & metaProxC$Brain_Region == "CA3_other_negs" & metaProxC$FOS == "N" & metaProxC$Context1 == "none" & metaProxC$outliers == "in" ,"Sample_ID"]
#metaProxC$Brain_Region == "CA3_other_negs" & metaProxC$Mouse_condition == "HC" & metaProxC$alignable >  100000 & metaProxC$Smartseq2_RT_enzyme_used == "ProtoscriptII","Sample_ID"]
dat <- na.exclude(tpmProxC[, samples])
met <- metaProxC[match(samples,metaProxC$Sample_ID),]

colnames(dat) <- paste(met$Brain_Region,met$FOS,c(1:ncol(dat)),sep = ".")
plot(hclust(dist(t(dat[genes,]))))
k <- cutree(hclust(dist(t(dat[genes,]))),2)#cutree(hclust(daisy([-c(25),c(1:3,19:20,39:73)])),k = 5)

table(metaProxC[names(k),"Brain_Region"], k)

metaProxC[names(k[k==1]),"Subgroup2"] <- "Neg.1"
metaProxC[names(k[k==2]),"Subgroup2"] <- "Neg.2"
