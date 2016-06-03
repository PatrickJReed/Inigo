#####################
###Comparing Clustering
#####################
library(ggplot2)
library(pcaMethods)
library(Rtsne)
library(WGCNA)
#####################
## Step One: Choose Data
#####################
samples <- metaProxC[ metaProxC$Brain_Region == "DG" & metaProxC$alignable >  100000 & metaProxC$Smartseq2_RT_enzyme_used == "ProtoscriptII" ,"Sample_ID"]
dat <- tpmProxC[, samples]
met <- metaProxC[match(samples,metaProxC$Sample_ID),]
#####################
## Step Two: Calculate Clusters by Clustering or DimRed then Clustering
#####################
K <- 4
# Linear: PCA
p <- pca(t(dat),nPcs = 4)
p.k <- as.numeric(kmeans(p@scores, centers = K)$cluster)
# Non-Linear: t-SNE
TSNE <- Rtsne(as.matrix(t(na.exclude(dat))),initial_dims=50,perplexity=10,theta=0.1,check_duplicates=FALSE)
t <- as.data.frame(TSNE$Y)
t.k <- as.numeric(kmeans(t,K)$cluster)
# Merged Linear/Non-linear 
reduced <- data.frame(p@scores[,1:3],t)
d <- daisy(reduced, metric = "gower")
d.k <- as.numeric(cutree(hclust(d),k = K))
# WGCNA
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",
     main = paste("Scale independence"));
#choose the power for the experiment
POWER = 20
#Calculate the network
net = blockwiseModules(datExpr,maxBlockSize=60, power = POWER,
                       TOMType = "unsigned", minModuleSize = 2,
                       reassignThreshold = 0, numericLabels = TRUE, pamRespectsDendro = FALSE,
                       TOMDenom="min")
# Convert labels to colors for plotting
w.k = as.numeric(as.factor(labels2colors(net$colors)))

#####################
## Step Three: Test ability of the clusters to separate out Protein Labels and Arc, Fos, Inhba
#####################
