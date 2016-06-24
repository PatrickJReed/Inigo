############
## Negative Subtypes
############
#load(c("~/Documents/SalkProjects/ME/ShortLongSingature/SLSig_R/tsne.rda"))
# Above saved in the Visualizations R script
############
library(pcaMethods)
library(edgeR)
library(scatterplot3d)
library(fdrtool)
library(ggplot2)
### other functions
GLM <- function(dat, variable1, variable2 = NULL, prefit=FALSE){
  if(is.null(variable2)){
    design <- model.matrix(~variable1)
    Coef = 2
  }else{
    design <- model.matrix(~variable2*variable1)
    Coef = 4
  }
  cds <- DGEList(dat)
  if (prefit == FALSE){
    cds <- calcNormFactors(cds)
    cds <- estimateGLMCommonDisp( cds )
    cds <- estimateGLMTrendedDisp(cds)
    fit <- glmQLFit(y=cds,design)
  }
  fit <- glmQLFit(y=cds,design)
  lrt <- glmLRT(fit, coef=Coef)#interaction term = 4
  edg <- data.frame(lrt$table)
  edg <- edg[order(edg$PValue),]
  f <- p.adjust(edg$PValue)
  edg$f <- f
  return(edg)
}
exact <- function(dat, variable, Pair){
  ########################
  ### Pre-process Data
  ########################
  cds <- DGEList(dat,group=group)
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
tsneplot <- function(t,gene = NULL){
  
  if (!is.null(gene)){
    t$gene <- as.numeric(tpmProxC[gene,rownames(t)])
    p <-  ggplot(t, aes(T1,T2, colour = gene))
  }else{
      p <-  ggplot(t, aes(T1,T2, colour = k))
  }
  p <- p + geom_point(size = 4, alpha = 0.7) +
     geom_point(size = 4,  shape = 1) +theme_bw()+
    xlab("TSNE1")+
    ylab("TSNE2")+
    labs(title="Homecage Negs\nt-SNE")+
    theme(text=element_text(size=20))+
    theme(panel.border = element_rect(colour=c("black"),size=2),
          axis.ticks = element_line(size=1.5),
          panel.grid.major = element_line(size = 1))
  
print(p)
}
############
# Assign subgroups using k-means
#############
#k <- kmeans(t.hc.neg[,c(1,2)],centers = 5,100)
#t.hc.neg$k <- as.factor(k$cluster)
#plot the hc plot
#############
# Plot t-SNE
#############
tsneplot(t.hc.neg,gene = "Sst")
#############
## EDGER
#############
t.hc.neg$k <- as.numeric(as.character(t.hc.neg$k))
samples <- rownames(t.hc.neg)#[t.hc.neg$k == 1|t.hc.neg$k == 2 ,])
dat <- na.exclude(countProxC[, samples])
dat <- dat[rowSums(dat) > 0,]
met <- metaProxC[match(samples,metaProxC$Sample_ID),]
#
group <- t.hc.neg$k == 4
Pair <- levels(as.factor(as.character(group)))
res <- exact(dat, group, Pair)
##########
# Combine Results
#RES.neg <- list()
i <- 2
RES.neg[[i]] <- res
names(RES.neg)[i] <- "hc_FOSN_neg_subpop_1_notvs5_4"

###########
# PCA for linear separation of putative CA3 cells
###########
t.hc.neg$k <- as.numeric(as.character(t.hc.neg$k))
samples <- rownames(t.hc.neg)#rownames(t.hc.neg[t.hc.neg$k <4 ,])
dat <- na.exclude(tpmProxC[, samples])
dat <- dat[rowSums(dat) > 0,]
met <- metaProxC[match(samples,metaProxC$Sample_ID),]
p <- pca(t(dat[,samples]),nPcs = 10)
scores <- as.data.frame(p@scores)
scores <- cbind(scores, met)
loading <- as.data.frame(p@loadings)
#tiff(filename = "~/Documents/SalkProjects/ME/ShortLongSingature/SLSig_tiff/PCA_HC_N.tiff",width = 6.5,height = 5,units = 'in',res = 300)
ggplot(scores,aes(PC1, PC2, colour =AMP_Date))+# as.factor(t.hc.neg[samples,"k"])))+
  geom_point(size = 4)

#############
## SUbgroup t-SNE
#############
samples <- rownames(t.hc.neg[t.hc.neg$k <4 ,])
dat <- na.exclude(tpmProxC[, samples])
dat <- dat[rowSums(dat) > 0,]
met <- metaProxC[match(samples,metaProxC$Sample_ID),]


i <- 2#17
TSNE <- Rtsne(as.matrix(t(na.exclude(dat))),initial_dims=5,perplexity=i,theta=0,check_duplicates=FALSE,dims = 2)
t <- as.data.frame(TSNE$Y)
colnames(t) <- c("T1","T2")#,"T3")
t <- cbind(t,met)
k <- kmeans(t[,c(1,2)],centers = 5,100)
t$k <- as.factor(k$cluster)
#plot the hc plot
tsneplot(t)



