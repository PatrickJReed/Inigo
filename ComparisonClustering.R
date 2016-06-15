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
library(pvclust)
#####################
## Step One: Choose Data
#####################

samples <- metaProxC[ metaProxC$FOS == "N" & metaProxC$Mouse_condition == "HC" & metaProxC$Context != "A" & metaProxC$alignable >  100000 & metaProxC$Smartseq2_RT_enzyme_used == "ProtoscriptII" ,"Sample_ID"]
samples <- samples[-c(match(outliers,samples))]
dat <- tpmProxC[, samples]
met <- metaProxC[match(samples,metaProxC$Sample_ID),]
#####################
## Step Two: Calculate the number of useful components
#####################

pc.it <- vector()
dat <- na.exclude(dat)
  # Linear: PCA
for (n in seq(2,25,2)){  
    j <- j+1
    p <- pca(t((dat)),nPcs = 30)
    pv <- pvclust(p@scores[,c(1:n)],nboot = 500)
    pv <- pv$edges
    pv$nPC <- n
    #Assign
    pc.it <- rbind(pc.it,pv )
}


ggplot(pc.it, aes(v,c, colour = as.factor(nPC)))+
  geom_point()





  # Non-Linear: t-SNE
  for (per in seq(3,12,3)){
    for (K in seq(3,15,3)){
      TSNE <- Rtsne(as.matrix(t((dat))),initial_dims=7,perplexity=per,theta=0,dims = 3,check_duplicates=FALSE)
      t <- as.data.frame(TSNE$Y)
      t.k <- as.numeric(kmeans(t,K)$cluster)
      #Assign
      j <- j+1
      nm <- paste("t",K,per,sep = ".")
      group[[j]] <- t.k
      names(group)[j] <- nm
    }
  }
    for (K in seq(3,15,3)){
      k.k <- as.numeric(kmeans(t(dat),K)$cluster)
      #Assign
      j <- j+1
      nm <- paste("k",K,per,sep = ".")
      group[[j]] <- k.k
      names(group)[j] <- nm
    }
  for (K in seq(3,15,3)){
    h.k <- as.numeric(cutree(hclust(dist(t(dat))),K))
    #Assign
    j <- j+1
    nm <- paste("h",K,per,sep = ".")
    group[[j]] <- h.k
    names(group)[j] <- nm
  }
  # Merged Linear/Non-linear 

#####################
## Step Three: Test ability of the clusters to separate out Protein Labels and Arc, Fos, Inhba
#####################

score <- list()

genes <- c("Prox1","Vip","Gad1","Meg3")
proteins <- c("PROX1","CTIP2")
for (method in names(group)){
  if (is.null(score[[method]])){
    score[[method]] <- 0
  }
  for(gene in genes){
    model <- anova(lm(as.numeric(dat[gene,]) ~ as.factor(group[[method]])))
    score[[method]] <- score[[method]] + model[1,4]
  }
  for(protein in proteins){
    model <- anova(lm(as.numeric(met[,protein]) ~ as.factor(group[[method]])))
    score[[method]] <- score[[method]] + model[1,4]
  }
}

#
score.df <- melt(as.data.frame(score))
score.df[grep("p",score.df$variable),"method"] <- "pca"
score.df[grep("t",score.df$variable),"method"] <- "tsne"
score.df[grep("d",score.df$variable),"method"] <- "joint"
score.df[grep("k",score.df$variable),"method"] <- "kmeans"
score.df[grep("h",score.df$variable),"method"] <- "hclust"
score.df$k <- seq(3,15,3)
score.df$i <- do.call("rbind",strsplit(as.character(score.df$variable), split = ".",fixed = TRUE))[,3]

#tiff("~/Documents/SalkProjects/ME/ShortLongSingature/SLSig_tiff/DimRed.tiff",width = 8,height = 7,res = 300,units = 'in')
ggplot(score.df[score.df$method == "pca",], aes(k, value, colour = as.factor(i)))+
  geom_point(size = 4)+
  geom_line()+
  labs(title= "Iterated Dimensionality Reduction\nPCA")+
  ylab("Supervised Clustering Score")
#dev.off()

#################