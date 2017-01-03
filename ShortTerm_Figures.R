###################
## Figures for Short-term
###################
library(ggplot2)
library(reshape)
library(plyr)
library(pcaMethods)
library(fastICA)

############################################
## boxplots of normalized tpm values
############################################

tmp <- data.frame(counts = colSums(countProx),sample = labelsProx$well, labels = labelsProx$fos, date = labelsProx$date)

ggplot(tmp, aes(reorder(sample,as.numeric(as.factor(labels))), counts,fill= labels))+
  geom_bar(stat="identity")+
  xlab("Sample")+
  ylab("Count")+
  labs(title="Raw Counts")+
  geom_hline(yint = mean(tmp$counts),linetype="dashed")+
  theme_bw(base_size=20)+
  theme(panel.border = element_rect(colour=c("black"),size=2),
        axis.ticks = element_line(size=1.5))+
  theme(axis.text.x = element_text(angle = 70, hjust = 1))
  
  
############################################
## Mean Expression by Variance
############################################
tmp <- na.exclude(data.frame(mean = apply(X=tpmProx,MARGIN=1,FUN=meanNoZero), var = apply(X=tpmProx,MARGIN=1,FUN=varNoZero)))

ggplot(na.exclude(tmp), aes(mean,var))+
  geom_point()+
  xlab("Mean Expression")+
  ylab("Variance")+
  labs(title="Normalized Prox1+ Expression")+
  theme_bw(base_size=20)+
  theme(panel.border = element_rect(colour=c("black"),size=2),
        axis.ticks = element_line(size=1.5))

############################################
## Mean Expression by P(detection)
############################################
tmp <- na.exclude(data.frame(mean = apply(X=tpmProx,MARGIN=1,FUN=meanNoZero), prop = apply(X=tpmProx,MARGIN=1,FUN=propExp)))

ggplot(na.exclude(tmp), aes(mean,prop))+
  geom_point()+
  xlab("Mean Expression")+
  ylab("P(detection)")+
  labs(title="Normalized Prox1+ Expression")+
  theme_bw(base_size=20)+
  theme(panel.border = element_rect(colour=c("black"),size=2),
        axis.ticks = element_line(size=1.5))

############################################
## Hierarchical Cluster
############################################
dat <- tpmProx
colnames(dat) <- labelsProx$fos
h <- hclust(dist(t(dat)))
plot(h)
Hgroups <- cutree(tree=h,k=2)
############################################
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
## SECTION 1: Prox1+ Homecage
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
############################################
############################################
## PCA of Prox1+ HC
############################################
dat <- tpmProx[apply(X=tpmProx,MARGIN=1,FUN=meanNoZero) > 0,]
p <- pca(t(dat))
l <- as.data.frame(p@loadings)
l <- l[order(l$PC1,decreasing=TRUE),]
g <- "Arc"
tmp <- data.frame(PC1 = p@scores[,1], PC2 = p@scores[,2],
                  fos = labelsProx$fos,
                  gene = as.numeric(tpmProx[g,]))
#pdf("~/Documents/SalkProjects/BenLacar/ManuscriptFigures/pca_nooutliers.pdf",width=7,height=7)
ggplot(tmp, aes(PC1, PC2,shape= as.factor(Hgroups),colour = conc))+
  geom_point(size=5, alpha=0.7)+
  #geom_point(size=5, shape=1, colour="black")+
  theme_bw()+
  labs(title="Variability Between Samples")+
  xlab(paste("PC1: ", 100 * round(p@R2[1],2),"%",sep=""))+
  ylab(paste("PC2: ", 100 * round(p@R2[2],2),"%",sep=""))+
  scale_colour_gradient2(high="red",mid="grey",low="blue",midpoint=0.5)
#dev.off()




############################################
## Prox1 HC Individual Genes (tpm)
############################################
gene <- "Per1"
tmp <- data.frame(exp = as.numeric(tpmProx[gene,]),
                  Fos = as.factor(labelsProx$fos),
                  date = as.factor(labelsProx$date))

#pdf("~/Documents/SalkProjects/BenLacar/ManuscriptFigures/Camk4.pdf",width=6,height=5)
ggplot(tmp, aes(Fos,exp,colour=conc))+
  geom_boxplot(outlier.shape=NA)+
  geom_point(position=position_jitter(width=0.1,height=0))+
  theme_bw(base_size=20)+
  ylim(c(0,max(tmp$exp) + 0.2*(max(tmp$exp))))+
  theme(panel.border = element_rect(colour=c("black"),size=2),
        axis.ticks = element_line(size=1.5))+
  labs(title=paste(gene))
#dev.off()


############################################
## Prox1 HC Individual Genes by Monocle Grouping (tpm)
############################################
gene <- "Htr1a"
tmp <- data.frame(exp = as.numeric(tpmProx[gene,]),
                  Fos = as.factor(labelsProx$fos),
                  date = as.factor(labelsProx$date),
                  state  = ICA.states[colnames(tpmProx),"State"])

#pdf("~/Documents/SalkProjects/BenLacar/ManuscriptFigures/Camk4.pdf",width=6,height=5)
ggplot(tmp, aes(state,exp,colour=Fos))+
  geom_boxplot(outlier.shape=NA)+
  geom_point(position=position_jitter(width=0.1,height=0))+
  theme_bw(base_size=20)+
  ylim(c(0,max(tmp$exp) + 0.2*(max(tmp$exp))))+
  theme(panel.border = element_rect(colour=c("black"),size=2),
        axis.ticks = element_line(size=1.5))+
  labs(title=paste(gene))
#dev.off()


############################################
## Heatmap
############################################
dat <- tpmQC
pgenes <- apply(X=dat,MARGIN=1,FUN=propExp)
genes <- rownames(dat[pgenes > 0.2 & pgenes < 0.9,])
tmp <- dat[genes,]
tmp[tmp == 0] <- NA
tmp <- apply(X=tmp,MARGIN=2,FUN=scaleME)
tmp2 <- melt(t(tmp))

#sample order 
p <- apply(X=dat,MARGIN=2,FUN=propExp)
sampleorder <- rep(p,times=nrow(tmp))
#gene order 
p <- apply(X=dat[genes,],MARGIN=1,FUN=propExp)
geneorder <- rep(p,each=ncol(tmp))
#pdf("~/Documents/SalkProjects/BenLacar/ManuscriptFigures/Camk4.pdf",width=6,height=5)
ggplot(tmp2, aes(reorder(X1,sampleorder),reorder(X2,geneorder),fill=value))+
  geom_tile()+
  scale_fill_gradient(high="red",low="blue")#+
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),legend.position="none",
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank())

############################################
## DESeq Plot
############################################

difdat <- as.data.frame(resOrdered)
difdat$sig <- difdat$padj < 0.05

ggplot(difdat, aes(log2FoldChange, -log(pvalue),colour = sig))+
  geom_point()+
  xlim(-4.5,4.5)+
  labs(title="Grouped by hcluster k = 2 All Genes")+
  xlab("Log2FoldChange group1 < group2")+
  theme_bw(base_size=20)
