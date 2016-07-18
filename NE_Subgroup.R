################
## ID EE subtypes using genes from HC subtypes
################
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
Test <- function(term = "CA1",Dat, Met){
  if (term == "Neg"){
    group <- Met$Subgroup == "CA3" | Met$Subgroup == "IN"
  }else{
    group <- Met$Subgroup == term
  }
  if(term == "Sub"){
    group <- Met$Subgroup2 == "IN" & Met$Brain_Region == "CA1"
  }
  Pair <- levels(as.factor(as.character(group)))
  res <- exact(Dat, group, Pair)
  a <- apply(tpmProxC[,colnames(dat)[group == TRUE]],1,rawExp,1)
  b <- apply(tpmProxC[,colnames(dat)[group == FALSE]],1,rawExp,1)
  res$a <- a[rownames(res)]
  res$b <- b[rownames(res)]
  genes <- rownames(res[res$logFC > 0 & res$f < 0.01 & res$a > (sum(group) * (1/2)) &  res$b < (sum(group == FALSE)/2),])
  
  return(genes)
  
}
SubsampleEdgeR <- function(samples.2, samples.1, condition = "DG" ){
  genes <- vector()
  for (i in 1:20){
    samples <- c(samples.2, sample(x = samples.1,size = length(samples.2),replace = FALSE))
    dat <- na.exclude(countProxC[, samples])
    dat <- dat[rowSums(dat) > 0,]
    met <- metaProxC[match(samples,metaProxC$Sample_ID),]
    ###################
    #Assign groups
    ###################
    group <- met$Subgroup == condition
    Pair <- levels(as.factor(as.character(group)))
    ###################
    # Test genes
    ###################
    #!!!Run exact test
    res <- exact(dat, group, Pair)
    genes <- c(genes, rownames(res[res$logFC > 0 & res$f < 0.05,]))
  }
  
  genes2 <- as.data.frame(table(genes))
  genes2 <- genes2[order(genes2$Freq,decreasing = TRUE),]
  return(genes2)
}

#################
## Dataset for major groups
#################
samples <- rownames(metaProxC[  metaProxC$Context1 == "none" & metaProxC$outliers == "in" & metaProxC$Subgroup != "Unk" & metaProxC$FOS == "N" & metaProxC$Mouse_condition == "HC" ,])
dat <- na.exclude(countProxC[, samples])
dat <- dat[rowSums(dat) > 0,]
met <- metaProxC[match(samples,metaProxC$Sample_ID),]
### Get genes for major groups
ca1.genes <- Test(term = "CA1",Dat = dat, Met = met)
dg.genes <- Test(term = "DG",Dat = dat, Met = met)
vip.genes <- Test(term = "VIP",Dat = dat, Met = met)
neg.genes <- Test(term = "Neg",Dat = dat, Met = met)
ca3.genes <- Test(term = "CA3",Dat = dat, Met = met)
in.genes <- Test(term = "IN",Dat = dat, Met = met)

marker.genes <- c("Calb2","Cacng5","Fibcd1","Iyd","Coch","Wfs1","Dcn")
#
all.genes <- unique(c(ca1.genes,dg.genes,vip.genes,ca3.genes,in.genes,neg.genes))
#################
## Predict based on these genes
#################
#################
## Note: The MAJOR groups of CA1, DG, VIP, and Negs, need to be assessed separately from Negs alone
#################
# First test HC to make sure they're correct
samples <- rownames(metaProxC[  metaProxC$Context1 == "none" & metaProxC$Subgroup != "Unk" & metaProxC$outliers == "in" & metaProxC$FOS == "N" & metaProxC$Mouse_condition == "HC" ,])
dat <- na.exclude(tpmProxC[, samples])
met <- metaProxC[match(samples,metaProxC$Sample_ID),]
colnames(dat) <- paste(met$Subgroup,1:ncol(dat),sep = ".")
#Random Forest
all.genes <- unique(c(marker.genes,ca1.genes,dg.genes,vip.genes,ca3.genes,in.genes,neg.genes))
dat2 <- t(dat[all.genes, ])
dat2 <- data.frame(dat2)
dat2$Subgroup.1 <- as.character(met$Subgroup)
dat2$Subgroup.1 <- as.factor(dat2$Subgroup.1)
rf.model <- randomForest(Subgroup.1 ~ . , dat2)
gini <- rf.model$importance
gini <- gini[order(gini,decreasing=TRUE),]
############
#Run again to refine the gene list
############
next.genes <- names(gini[gini > quantile(gini, 0.7)])
dat2 <- t(dat[next.genes, ])
dat2 <- data.frame(dat2)
p <- apply(dat2, 2,rawExp)
dat2 <- dat2[,p > 0]
dat2$Subgroup.1 <- as.character(met$Subgroup)
dat2$Subgroup.1 <- as.factor(dat2$Subgroup.1)
rf.model <- randomForest(Subgroup.1 ~ . , dat2)
gini <- rf.model$importance
gini <- gini[order(gini,decreasing=TRUE),]

# Now predict with adding in Protein Stain
samples <- rownames(metaProxC[  metaProxC$Context1 == "none" & metaProxC$outliers == "in" ,])
dat <- na.exclude(tpmProxC[, samples])
met <- metaProxC[match(samples,metaProxC$Sample_ID),]

#Random Forest
dat2 <- t(dat[next.genes, ])
dat2 <- data.frame(dat2)
pred.all <- predict(rf.model,newdata = dat2)
table(pred.all,met$Brain_Region)
#Check HC again
samples <- rownames(metaProxC[  metaProxC$Context1 == "none" & metaProxC$Subgroup != "Unk" & metaProxC$outliers == "in" & metaProxC$FOS == "N" & metaProxC$Mouse_condition == "HC" ,])
met <- metaProxC[match(samples,metaProxC$Sample_ID),]
table(data.frame(met$Brain_Region, pred.all[samples]))
#Note: HC is Perfect!!!!
#metaProxC[names(pred.all),"Subgroup2"] <- as.character(pred.all)
############
# Now predict
samples <- rownames(metaProxC[  metaProxC$Context1 == "none" & metaProxC$outliers == "in" ,])
dat <- na.exclude(tpmProxC[, samples])
met <- metaProxC[match(samples,metaProxC$Sample_ID),]
dat2 <- t(dat[next.genes, ])

p <- apply(dat2, 2,rawExp)
dat2 <- dat2[,p > 0]
#dat2 <- dat2[,-c(grep("Rik",all.genes))]
dat2 <- as.data.frame(dat2)
dat2$Subgroup.1 <- as.character(met$Subgroup2)
dat2$Subgroup.1 <- as.factor(dat2$Subgroup.1)
rf.model <- randomForest(Subgroup.1 ~ . , dat2)
gini <- rf.model$importance
gini <- gini[order(gini,decreasing=TRUE),]


#Random Forest
samples <- names(pred.all[pred.all == "CA1" | pred.all == "CA3" ])
met <- metaProxC[samples,]
samples <- rownames(met[met$FOS == "N" & met$Subgroup != "Unk" & met$Mouse_condition == "HC" ,])
#samples.1 <- rownames(met[met$Subgroup == "CA1",])
#samples.2 <- rownames(met[met$Subgroup == "CA3",])
#samples.1 <- sample(x = samples.1,size = length(samples.2))
#samples <- c(samples.1,samples.2)
dat <- na.exclude(countProxC[, samples])
met <- metaProxC[match(samples,metaProxC$Sample_ID),]

all.genes <- unique(c(ca3.genes))
dat2 <- t(dat[all.genes, ])
dat2 <- data.frame(dat2)
dat2$Subgroup.1 <- as.character(met$Subgroup)
dat2$Subgroup.1 <- as.factor(dat2$Subgroup.1)
rf.model2 <- randomForest(Subgroup.1 ~ . , dat2, replace = FALSE, ntree = 10000)
rf.model2
gini <- rf.model2$importance
gini <- gini[order(gini,decreasing=TRUE),]

###
samples <- names(pred.all[pred.all == "CA1" | pred.all == "CA3" ])
met <- metaProxC[match(samples,metaProxC$Sample_ID),]
#samples <- rownames(met[met$Brain_Region == "CA3_other_negs",])
dat <- na.exclude(countProxC[, samples])
met <- metaProxC[match(samples,metaProxC$Sample_ID),]

all.genes <- unique(c(ca3.genes))
dat2 <- t(dat[all.genes, ])
dat2 <- data.frame(dat2)
pred.second <- predict(rf.model2,newdata = dat2)
table(data.frame(pred.second,met$Brain_Region))


#metaProxC[names(pred.second),"Subgroup2"] <- as.character(pred.second)

