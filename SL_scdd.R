#library(devtools)
#devtools::install_github("kdkorthauer/scDD")
library(scDD)
library(Biobase)

#construct starting object
samples <- rownames(t[t$k == 5  | t$k == 8,])
cts <- as.matrix(countProxC[,samples])
pheno <- data.frame(condition = as.factor(ifelse(as.numeric(tpmProxC["Vip",samples])  >7, yes = "high","low")),row.names = samples)

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
