##########
## Error Rate
##########
library(pcaMethods)
library(edgeR)
library(scatterplot3d)
library(fdrtool)
library(ggplot2)
library(randomForest)
### other functions
exact <- function(dat, variable, Pair){
  ########################
  ### Pre-process Data
  ########################
  cds <- DGEList(dat,group=variable)
  # Filter out really lowly detected reads (*from EdgeR tutorial)
  cds <- cds[rowSums(1e+06 * cds$counts/expandAsMatrix(cds$samples$lib.size, dim(cds)) > 1) >= 3, ]
  cds <- calcNormFactors( cds )
  ########################
  ### Estimate parameters from data
  ########################
  cds <- estimateCommonDisp( cds )
  cds <- estimateTagwiseDisp( cds )
  cds <- estimateTrendedDisp(cds)
  ## Run the test
  de.cmn <- exactTest( cds , pair = Pair)
  ########################
  ### Format Results
  ########################
  res <- as.data.frame(de.cmn$table)
  res <- res[order(res$PValue),]
  #f <- fdrtool(x=res$PValue,statistic="pvalue",plot=FALSE)
  res$f <- p.adjust(res$PValue, method = "fdr")
  return(res)
}
##############
#Get Samples
##############
samples <- metaProxC[ metaProxC$Mouse_condition == "HC" & metaProxC$Subgroup != "Unk" & metaProxC$FOS == "N" & metaProxC$Context1 == "none" & metaProxC$outliers == "in" ,"Sample_ID"]#
dat <- na.exclude(countProxC[, samples])
met <- metaProxC[match(samples,metaProxC$Sample_ID),]
##############
#Get CA1,CA3,VIP,DG,and IN genes
##############
res.ca1 <- exact(dat, met$Subgroup2 == "CA1", c(FALSE,TRUE))
genes.ca1 <- rownames(res.ca1[res.ca1$f < 0.05,])[1:50]
res.ca3 <- exact(dat, met$Subgroup2 == "CA3", c(FALSE,TRUE))
genes.ca3 <- rownames(res.ca3[res.ca3$f < 0.05,])[1:50]
res.vip <- exact(dat, c(met$Subgroup == "Vip" | met$Subgroup == "VIP.2"), c(FALSE,TRUE))
genes.vip <- rownames(res.vip[res.vip$f < 0.05,])[1:50]
res.dg <- exact(dat, c(met$Subgroup2 == "DG"  | met$Subgroup2 == "DG.2") , c(FALSE,TRUE))
genes.dg <- rownames(res.dg[res.dg$f < 0.05,])[1:50]
res.in <- exact(dat, c(met$Subgroup2 == "IN.1" | met$Subgroup2 == "IN_Erbb4" | met$Subgroup2 == "IN.Ptger1"), c(FALSE,TRUE))
genes.in <- rownames(res.in[res.in$f < 0.05,])[1:50]

#############
# Random Forest
#############
samples <- metaProxC[ metaProxC$Mouse_condition == "HC" & metaProxC$Subgroup != "Unk" & metaProxC$FOS == "N" & metaProxC$Context1 == "none" & metaProxC$outliers == "in" ,"Sample_ID"]#
dat <- na.exclude(countProxC[, samples])
met <- metaProxC[match(samples,metaProxC$Sample_ID),]
#
genes <- unique(c(genes.ca1, genes.vip,genes.ca3, genes.dg, genes.in))
dat2 <- dat[rownames(dat), ]
p <- apply(dat2, 1, rawExp, 2)
dat2 <- t(dat2[p > 6,])
dat2 <- data.frame(dat2)
dat2$Subgroup <- (met$Subgroup2)
#dat2[met$Subgroup2 == "CA3","Subgroup"] <- "CA3" 
dat2$Subgroup <- as.factor(dat2$Subgroup)
#dat2$treatment <- met$Treatment
rf.model <- randomForest(Subgroup ~ . , dat2)
gini <- rf.model$importance
gini <- gini[order(gini,decreasing=TRUE),]

