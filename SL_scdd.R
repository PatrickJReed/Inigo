#library(devtools)
#devtools::install_github("kdkorthauer/scDD")
library(scDD)
library(Biobase)
#save(list = c("metaProxC", "countProxC"),file = "~/Documents/SalkProjects/ME/ShortLongSingature/test.rda",compress = TRUE)
load("~/Documents/SalkProjects/ME/ShortLongSingature/test.rda")

#construct starting object
samples <- rownames(metaProxC[metaProxC$Subgroup2 == "DG" & metaProxC$FOS != "L"  & metaProxC$Arc_2.5 != "greater" & metaProxC$Context1 == "none" & metaProxC$cluster_outlier == "in" & metaProxC$outliers == "in",])
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
RES.sc <- list()
i <- 0
#for (alpha in c(0.001, 0.01, 0.05 )){ running through these different alphas didn't change the results at all...why?
# original: prior_param=list(alpha=0.01, mu0=0, s0=0.01, a0=0.01, b0=0.01)
  prior_param=list(alpha, mu0=0, s0=0.01, a0=0.01, b0=0.01)
  #Detect and classify differentially distributed genes
  set.seed(6767) # set random seed for reprodicbility
  dd.results <- scDD(eset.scdd, prior_param=prior_param, permutations=0, testZeroes=TRUE)
  res <- dd.results$Genes
  res <- res[order(res$nonzero.pvalue.adj),]
  rownames(res) <- res$gene
  nm <- paste("sc",alpha, sep = ".")
  i <- i + 1
#  RES.sc[[i]] <- res
#  names(RES.sc[[i]]) <- nm
#}


### Plot per gene
g <- "Wwp1"
sideViolin(exprs(eset.scdd)[rownames(eset.scdd) ==  g,], phenoData(eset.scdd)$condition,
           title.gene=g)

#res.vip_In <- res
#eset.scdd.vip_In <- eset.scdd
