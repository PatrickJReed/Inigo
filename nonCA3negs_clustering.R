#save(list = c("p.nonca3_all","p.nonca3_fosN", "df.nonca3_all","df.nonca3_fosN"), file = "~/Documents/SalkProjects/ME/ShortLongSingature/SLSig_R/nonCA3_negs_clust.rda")
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
samples <- rownames(t[t$a == 1 | t$a == 2,])
dat <- na.exclude(tpmProxC[, samples])
met <- metaProxC[samples,]

df <- vector()
for ( i in 1:100){
  i <- 8
  TSNE <- Rtsne(as.matrix(t(na.exclude(dat))),initial_dims=3,perplexity=i,theta=0,check_duplicates=FALSE,dims = 2,max_iter = 500)
  df <- rbind(df, TSNE$Y[,1], TSNE$Y[,2])
}
p <- pvclust(df)
plot(p)
a <- cutree(p$hclust,k=2)


###
i <- 10
TSNE <- Rtsne(as.matrix(t(na.exclude(df))),initial_dims=10,perplexity=i,theta=0,check_duplicates=FALSE,dims = 2,max_iter = 500)
t <- as.data.frame(TSNE$Y)
colnames(t) <- c("T1","T2")
t <- cbind(t,met)
#t$group <- paste(t$Brain_Region, t$Mouse3, sep =".")
t$gene <- as.numeric(tpmProxC["Dcn",rownames(t)])
#k <- kmeans(t[,c(1:2)],centers =3,nstart = 2000)
#t$k <- as.factor(k$cluster)
#tiff("~/Documents/SalkProjects/ME/ShortLongSingature/MolecDissec_Figs_Tables/Figures_vD/tsne_all_a.tiff",width = 9,height = 6.5,units = 'in',res = 600,compression = 'lzw')
t$a <- a
ggplot(t, aes(T1,T2, color = as.factor(a)  , shape = FOS))+
  geom_point(size = 5)+
  theme_bw()+
  xlab("TSNE1")+
  ylab("TSNE2")+
  theme(text=element_text(size=20))+
  theme(panel.border = element_rect(colour=c("black"),size=2),
        axis.ticks = element_line(size=1.5),
        panel.grid.major = element_line(size = 1))
  
############# 
# Identify genes associated with the subgroups
############# 
# get samples for difExp analysis
dat <- na.exclude(countProxC[,samples])
dat <- dat[rowSums(dat) > 0,]
met <- metaProxC[samples,]

# Cut the bootstrapped tree from above
h <- as.factor(cutree(p$hclust,k = 2)) 
group <- h == 1
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
group <- as.factor(cutree(p$hclust,k = 2))

dat <- na.exclude(tpmProxC[, samples])
met <- metaProxC[samples,]

g <- "Il34"
tmp <- data.frame(tpm = as.numeric(dat[g,]),met)
tmp[samples,"h"] <-  group

ggplot(tmp, aes(h, tpm, fill = h, alpha = 0.3))+
  geom_violin()+
  geom_jitter()

####
samples <- c(rownames(t[t$k == 4,]),
             rownames(tmp[tmp$h == 2,]))
group <- c(rep("CA3",sum(t$k == 4)),
           rep("unk",sum(tmp$h == 2))
)
