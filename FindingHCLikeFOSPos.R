#### Induction versus partition
HighExp <- function(g){
  return(as.numeric(dat[g,]) > as.numeric(Mean[g]))
}
LowExp <- function(g){
  return(as.numeric(dat[g,]) < as.numeric(Mean[g]))
}
# 1) Find genes higher in FOS+
a <- RES[["DG"]]
genes <- rownames(a[a$logFC > 0 & a$f < 0.001,])
# 2) Find homecage neurons with high expression of a majority of these genes
samples <- metaProxC[metaProxC$Mouse_condition == "HC" &  metaProxC$Subgroup2 == "DG" & metaProxC$FOS == "N" &  metaProxC$EE_ArcGroup != "Unk"  & metaProxC$Context1 == "none" & metaProxC$outliers == "in","Sample_ID"]
dat <- na.exclude(tpmProxC[, samples])
#Calcuate mean expression for each gene
Mean <- apply(dat[genes, ] ,1, meanNoZero)
#Find cells with high exp in homecage
samples <- metaProxC[ metaProxC$Mouse_condition == "HC" &  metaProxC$Subgroup2 == "DG" & metaProxC$FOS == "N" & metaProxC$Subgroup2!= "HDG" & metaProxC$Subgroup2!= "CA3" & metaProxC$Subgroup2!= "IN" &  metaProxC$EE_ArcGroup != "Unk" & metaProxC$Subgroup2 != "CA2" & metaProxC$Context1 == "none" & metaProxC$outliers == "in" ,"Sample_ID"]#
dat <- na.exclude(tpmProxC[, samples])
Raw <- data.frame(row.names = genes, do.call("rbind",(lapply(genes,HighExp))))
colnames(Raw) <- colnames(dat)
#####
#low exp
genes2 <- rownames(a[a$logFC < 0 & a$f < 0.001,])
#Calcuate mean expression for each gene
Mean <- apply(dat[genes2, ] ,1, meanNoZero)
#Find cells with low exp in homecage
Raw3 <- data.frame(row.names = genes2, do.call("rbind",(lapply(genes2,LowExp))))
colnames(Raw3) <- colnames(dat)
Raw <- rbind(Raw, Raw3)

TopSamples <- colnames(dat)[(as.numeric(colSums(Raw))/length(c(genes,genes2)))> 0.6]
# 3) Find genes that are on in all of these cells
Raw2 <-  Raw[,TopSamples]
TopGenes <- rownames(Raw2[rowSums(Raw2) > 8,])
#TopGenes.FOSP <- TopGenes
# 4) 



samples <- metaProxC[metaProxC$Mouse_condition == "HC" &  metaProxC$Subgroup2 == "DG" & metaProxC$FOS == "N" &  metaProxC$EE_ArcGroup != "Unk"  & metaProxC$Context1 == "none" & metaProxC$outliers == "in" |
                       metaProxC$Mouse_condition == "EE" &  metaProxC$Subgroup2 == "DG" & metaProxC$FOS != "L" &  metaProxC$EE_ArcGroup != "Unk"  & metaProxC$Context1 == "none" & metaProxC$outliers == "in" ,"Sample_ID"]
