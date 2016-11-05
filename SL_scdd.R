#library(devtools)
#devtools::install_github("kdkorthauer/scDD")
library(scDD)
library(Biobase)

#construct starting object
samples <- rownames(metaProxC[metaProxC$FOS != "L"  & metaProxC$Arc_2.5 != "greater" & metaProxC$Context1 == "none" & metaProxC$cluster_outlier == "in" & metaProxC$outliers == "in",])
dat <- na.exclude(countProxC[, samples])
dat <- dat[rowSums(dat) > 0,]
met <- metaProxC[samples,]
countTable2 <- t(apply(round(dat[rowSums(dat) > 0,],0),1, as.integer))
colnames(countTable2) <- colnames(dat)

cts <- as.matrix(countTable2[,samples])
pheno <- data.frame(condition = as.factor(met$FOS == "F"),row.names = samples)

eset.scdd <- ExpressionSet(assayData = cts,  
                           phenoData = as(pheno, "AnnotatedDataFrame"))

rm(cts)

#Create list of prior parameters for model fitting, why are these the defaults? Should play around with this
prior_param=list(alpha=0.01, mu0=0, s0=0.01, a0=0.01, b0=0.01)

#Detect and classify differentially distributed genes
set.seed(6767) # set random seed for reprodicbility
dd.results <- scDD(eset.scdd, prior_param=prior_param, permutations=0, testZeroes=FALSE)

res <- dd.results$Genes
res <- res[order(res$nonzero.pvalue.adj),]
rownames(res) <- res$gene

### Plot per gene
g <- "Baiap2"
sideViolin(exprs(eset.scdd)[rownames(eset.scdd) ==  g,], phenoData(eset.scdd)$condition,
           title.gene=g)

#res.vip_In <- res
#eset.scdd.vip_In <- eset.scdd
