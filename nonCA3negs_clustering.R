#save(list = c("p", "df"), file = "~/Documents/SalkProjects/ME/ShortLongSingature/SLSig_R/nonCA3_negs_clust.rda")
#load( "~/Documents/SalkProjects/ME/ShortLongSingature/SLSig_R/nonCA3_negs_clust.rda")
library(scDD)
library(Biobase)
library(scde)
prior_param=list(alpha=0.01, mu0=0, s0=0.01, a0=0.01, b0=0.01)
set.seed(6767) # set random seed for reprodicbility
n.cores <- 2 # for scde
############# 
# Identify sub-groups by iterative t-SNE
############# 
samples <- rownames(t[t$k == 3 & t$FOS == "N",])
dat <- na.exclude(tpmProxC[, samples])
met <- metaProxC[samples,]

df <- vector()
#for ( i in 1:1000){
  i <- 4
  TSNE <- Rtsne(as.matrix(t(na.exclude(dat))),initial_dims=3,perplexity=i,theta=0,check_duplicates=FALSE,dims = 2,max_iter = 500)
  df <- rbind(df, TSNE$Y[,1], TSNE$Y[,2])
}
p <- pvclust(df)

############# 
# Identify genes associated with the subgroups
############# 
# get samples for difExp analysis
dat <- na.exclude(countProxC[,samples])
dat <- dat[rowSums(dat) > 0,]
met <- metaProxC[samples,]

# Cut the bootstrapped tree from above
h <- as.factor(cutree(p$hclust,k = 5)) 
group <- h == 2 
Pair <- levels(as.factor(as.character(group)))
#edgeR
res.edge <- exact(dat, group, Pair)

#scde
dat2 <- apply(dat,2,function(x) {storage.mode(x) <- 'integer'; x}) 
o.ifm <- scde.error.models(counts=dat2,groups=group,n.cores=n.cores,threshold.segmentation=T,save.crossfit.plots=F,save.model.plots=F,verbose=1)
valid.cells <- o.ifm$corr.a >0
o.prior <- scde.expression.prior(models=o.ifm,counts=dat,length.out=400,show.plot=F)
groups <- as.factor(group)
names(group) <- row.names(o.ifm);
# run differential expression tests on all genes.
eset.scde <- scde.expression.difference(o.ifm,dat,o.prior,groups=group,n.randomizations=100,n.cores=n.cores,verbose=1)
eset.scde <- eset.scde[order(eset.scde$Z,decreasing=T),]

#scdd
pheno <- data.frame(condition = group,row.names = samples)
eset.scdd <- ExpressionSet(assayData = dat,phenoData = as(pheno, "AnnotatedDataFrame"))
dd.results <- scDD(eset.scdd, prior_param=prior_param, permutations=0, testZeroes=FALSE)
res.scde <- dd.results$Genes
res.scde <- res[order(res.scde$nonzero.pvalue.adj),]
rownames(res.scde) <- res.scde$gene


############# 
# Plot
############# 
group <- as.factor(cutree(p$hclust,k = 5))

dat <- na.exclude(tpmProxC[, samples])
met <- metaProxC[samples,]

g <- "Wdr26"
tmp <- data.frame(tpm = as.numeric(dat[g,]),met)
tmp[samples,"h"] <-  group

ggplot(tmp, aes(h, tpm, fill = h, alpha = 0.3))+
  geom_violin()+
  geom_jitter()

