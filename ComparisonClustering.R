#####################
###Comparing Clustering
#####################
#save(list = c("score.df.all","score.df.all"),file = "~/Documents/SalkProjects/ME/ShortLongSingature/SLSig_R/compareclustering.rda",compress = TRUE)
#load(c("~/Documents/SalkProjects/ME/ShortLongSingature/SLSig_R/compareclustering.rda"))
###
library(ggplot2)
library(reshape)
library(pcaMethods)
library(Rtsne)
library(WGCNA)
library(cluster)
#####################
## Step One: Choose Data
#####################

samples <- metaProxC[ metaProxC$FOS == "N" & metaProxC$Mouse_condition == "HC" & metaProxC$Context != "A" & metaProxC$alignable >  100000 & metaProxC$Smartseq2_RT_enzyme_used == "ProtoscriptII" ,"Sample_ID"]
dat <- tpmProxC[, samples]
met <- metaProxC[match(samples,metaProxC$Sample_ID),]
#####################
## Step Two: Calculate Clusters by Clustering or DimRed then Clustering
#####################

group <- list()
j <- 0
dat <- na.exclude(dat)
  # Linear: PCA
for (n in seq(2,22,5)){  
    for (K in seq(3,15,3)){
    j <- j+1
    p <- pca(t((dat)),nPcs = n)
    p.k <- as.numeric(kmeans(p@scores, centers = K)$cluster)
    #Assign
    nm <- paste("p",K,n,sep = ".")
    group[[j]] <- p.k
    names(group)[j] <- nm
  }
  # Non-Linear: t-SNE
  for (per in seq(3,12,3)){
    for (K in seq(3,15,3)){
      TSNE <- Rtsne(as.matrix(t((dat))),initial_dims=50,perplexity=per,theta=0.1,dims = 3,check_duplicates=FALSE)
      t <- as.data.frame(TSNE$Y)
      t.k <- as.numeric(kmeans(t,K)$cluster)
      #Assign
      j <- j+1
      nm <- paste("t",K,per,sep = ".")
      group[[j]] <- t.k
      names(group)[j] <- nm
    }
  }
  # Merged Linear/Non-linear 
}
#####################
## Step Three: Test ability of the clusters to separate out Protein Labels and Arc, Fos, Inhba
#####################

score <- list()

genes <- c("Prox1","Gad1","Sema3e","Meg3","Inhba")
proteins <- c("FOS")
for (method in names(group)){
  if (is.null(score[[method]])){
    score[[method]] <- 0
  }
  for(gene in genes){
    model <- anova(lm(as.numeric(dat[gene,]) ~ as.factor(group[[method]])))
    score[[method]] <- score[[method]] + model[1,4]
  }
 # for(protein in proteins){
#    model <- anova(lm(as.numeric(met[,protein]) ~ as.factor(group[[method]])))
#    score[[method]] <- score[[method]] + model[1,4]
#  }
}

#
score.df <- melt(as.data.frame(score))
score.df[grep("p",score.df$variable),"method"] <- "pca"
score.df[grep("t",score.df$variable),"method"] <- "tsne"
score.df[grep("d",score.df$variable),"method"] <- "joint"
score.df$k <- seq(3,15,3)
score.df$i <- rep(c(seq(2,22,5),seq(3,12,3)),each = length(seq(3,15,3)))

#tiff("~/Documents/SalkProjects/ME/ShortLongSingature/SLSig_tiff/DimRed.tiff",width = 8,height = 7,res = 300,units = 'in')
ggplot(score.df, aes(k, value, colour = as.factor(i),shape = as.factor(method)))+
  geom_point(size = 4)+
  geom_line()+
  labs(title= "Iterated Dimensionality Reduction\nHomecage FOSN")+
  ylab("Supervised Clustering Score")
#dev.off()

####
p <- pca(t(na.exclude(dat)),nPcs = 5)
TSNE <- Rtsne(as.matrix(t(na.exclude(dat))),initial_dims=50,perplexity=6,theta=0.1,check_duplicates=FALSE)
t <- as.data.frame(TSNE$Y)
t.k <- as.numeric(kmeans(t,K)$cluster)

reduced <- data.frame(p@scores,as.factor(t.k))
d <- daisy(reduced, metric = "gower")

s.2 <- vector()
for(K in 3:20){
d.k <- as.numeric(cutree(hclust(d),k = K))

s <- 0
for(gene in genes){
  model <- anova(lm(as.numeric(dat[gene,]) ~ as.factor(d.k)))
  s <- s + model[1,4]
}
s.2 <- c(s.2,s)
}


##
tmp <- melt(t(dat))
tmp$k <- as.factor(d.k)
tmp$FOS <- met$FOS
tmp$cond <- met$Mouse_condition
tmp$BR <- met$Brain_Region
gene <- "Arc"

ggplot(tmp[tmp$X2 == gene,], aes(k, value))+
  geom_violin()+
  labs(title = gene)

scores <- as.data.frame(p@scores)
scores <- cbind(scores,met)
ggplot(scores, aes(PC1,PC2, colour = FOS))+#colour = as.factor(group[["p.3.5"]])))+
  geom_point()
