############
## Cummulative Trait Detection
############
library(ggplot2)
library(parallel)
PCA.HT <- function(dat,nPC = 10){
  require(pcaMethods)
  p <- pca(t(dat),nPcs = nPC)
  return(p@scores)
} 

ICA.HT <- function(){
  require(fastICA)
  ica <- fastICA(t(dat),n.comp = 2)
}

TSNE.HT <- function(dat,perp = 8, K = 8){
  require(Rtsne)
  TSNE <- Rtsne(as.matrix(t(na.exclude(dat))),initial_dims=10,perplexity=perp,theta=0.1,check_duplicates=FALSE)
  t <- as.data.frame(TSNE$Y)
  k <- kmeans(t,K)
  return(as.character(k$cluster))
}

Distance <- function(reduced){
  require(cluster)
  d <- daisy(reduced, metric = "gower")
  return(d)
}

ModelMe <- function(g,tmp){
 model <- anova(lm(value ~ as.factor(k), tmp[tmp$X2 == g,]))
 return(as.vector(model$coefficients))
}
GeneDescriptor <- function(dat,k, genes){
  cl <- makeCluster(getOption("cl.cores", 3))
  tmp <- melt(t(dat))
  tmp$k <- k
  clusterExport(cl=cl, varlist=c("tmp"))
  tmp2 <- do.call("rbind",parLapply(cl=cl,X=genes, fun=ModelMe, tmp=tmp))
  stopCluster(cl)
  colnames(tmp2) <- paste(c("int",levels(as.factor(k))[-c(1)]),rep(c("Est","StdErr","tValue","pValue"),each= length(unique(k))),sep = ".")
  rownames(tmp2) <- genes
  return(tmp2)
}

############
# Find Hidden Traits
############
# Subset dataframe
samples <- metaProxC[metaProxC$Brain_Region == "DG" & metaProxC$alignable >  100000 &  metaProxC$Smartseq2_RT_enzyme_used == "ProtoscriptII" ,"Sample_ID"]#
dat <- na.exclude(tpmProxC[, samples])
met <- metaProxC[match(samples,metaProxC$Sample_ID),]

##########
p <- PCA.HT(dat)
t <- TSNE.HT(dat,perp = 10)
reduced <- data.frame(p[,1:3],
                      t,
                      as.numeric(apply(dat,2,propExp)), met$alignable)
k <- as.factor(cutree(hclust(Distance(reduced)),k = 6))


tmp <- data.frame(p,k,fos = met$FOS,prox = met$PROX1, cond = met$Mouse_condition, align = met$alignable,
                  gene = as.numeric(dat["Arc",]))
ggplot(tmp, aes(PC1,PC2, colour = k, shape = fos))+
  geom_point(size = 3,position = position_jitter(width = 0.8, height = 0.5))






tmp <- melt(t(dat))
tmp$k <- as.factor(k)
tmp$FOS <- met$FOS
tmp$cond <- met$Mouse_condition

gene <- "Fos"

ggplot(tmp[tmp$X2 == gene,], aes(k, value,colour = FOS))+
  geom_violin()+
  labs(title = gene)



geneTest <- GeneDescriptor(dat, k, genes[1:100] )
