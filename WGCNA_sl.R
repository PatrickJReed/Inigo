#WGCNA
#Note: WGCNA is not the right tool for 
#      analyzing gene networks OVER TIME
library("WGCNA")
library("DESeq2")
library("Rsamtools")

samples <- rownames(t[t$k == 3,])
dat <- tpmProxC[, samples]
met <- metaProxC[samples,]

v <- apply(dat, 1, var)
m <- apply(dat, 1, mean)
m <- rowSums(dat)
vstMat2 <- dat[names(na.exclude(v[v > 7 & m > 5])),]
##################################
#### Setup the WGCNA-Style Objects
##################################
#Define data set dimensions
datExpr <- as.data.frame(na.exclude(t(vstMat2)))
datTraits <- met[rownames(datExpr),]

##################################
#### Run WGCNA
##################################

options(stringsAsFactors = FALSE);
disableWGCNAThreads()
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 1, to=25, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
#plot the different power levels
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",
     main = paste("Scale independence"));
abline(h=0.90,col="red")
#choose the power for the experiment
POWER = 5
#Calculate the network
l <- vector()
for(x in seq(2,200,5)){
net = blockwiseModules(datExpr, power = POWER,maxBlockSize=5000,
                       TOMType = "unsigned", minModuleSize = 22,
                       reassignThreshold = 0, #mergeCutHeight = 0.99,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       TOMDenom="min",
                       #saveTOMs = TRUE,
                       #saveTOMFileBase = "hnp",
                       verbose = 3)
l <- c(l, length(unique(labels2colors(net$colors))))
}
#plot the network
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module Colors",
                    dendroLabels = FALSE,hang =0.03,
                    addGuide = TRUE, guideHang = 0.1)

#######
# Identify Association with trait
#######
moduleColors = labels2colors(net$colors)
MEsO = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs <- data.frame(orderMEs(MEsO))

heatmap(as.matrix(MEs))

colors <- c()
for (i in 1:nrow(MEs)){
  j <- as.numeric(MEs[i,(-c(which(colnames(MEs) == "MEgrey")))])
  j.max <- max(j)
  colors <- c(colors, colnames(MEs)[which(j == j.max)])
}

sample.colors <- data.frame(samples = colnames(t(datExpr)), colors = colors,row.names = colnames(t(datExpr)))


#identify labels/colors
moduleLabels = net$colors
table(moduleLabels)
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
modE <- moduleEigengenes(t(vstMat2),moduleColors)

# ##########################################
# Plot eigengenes
modE2 <- melt(t(modE[[1]]))
sample <- 1

ggplot(modE2, aes(X1, value,colour = X2))+
  geom_point()
# ##########################################

# Recalculate topological overlap if needed
enableWGCNAThreads()
TOM = TOMsimilarityFromExpr((datExpr), power = POWER);
# Read in the annotation file
#annot = read.csv(file = "GeneAnnotation.csv");
# Select modules
head(moduleLabels[which(moduleColors == "lightgreen")])
modules = c(18)#, "red");
# Select module probes
probes = rownames(vstMat2)
inModule = is.finite(match(moduleLabels, modules));
modProbes = probes[inModule];
modGenes <- modProbes
write.table(x=modGenes, file="~/Documents/test.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)
#modGenes = annot$gene_symbol[match(modProbes, annot$substanceBXH)];
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) <- list(modProbes)
colnames(TOM) <- (probes)
rownames(TOM) <- (probes)
dimnames(TOM) <- list(probes)

# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
                               weighted = TRUE,
                               threshold = 0.01,
                               nodeNames = modProbes,
                               nodeAttr = moduleColors[inModule]);
write.table(na.exclude(cyt$edgeData),"~/Documents/test.txt",quote=FALSE,row.names=FALSE)

##########
tmp <- as.data.frame(cyt$edgeData)
ggplot(tmp , aes(reorder(fromNode,weight), reorder(toNode,weight),fill=weight))+
  geom_tile()+
  scale_fill_continuous(high="red",low="blue")+
  theme(axis.text.x=element_text(angle=50,hjust=1,size=3),
        axis.text.y=element_text(angle=50,hjust=1,size=3))

###################################

#####
tmp2 <- melt(scale(t(na.exclude(datExpr))))
tmp2$module <- moduleColors
tmp2$br <- rep(met$Brain_Region, each = ncol(datExpr)) 
#line
m <- "red"
ggplot(tmp2[tmp2$module == m ,], aes(br, value))+
  geom_violin()+
  geom_jitter()+
  labs(title = m)+
  theme_bw(base_size = 18)

unique(as.character(tmp2[tmp2$module == m, "X1"]))

