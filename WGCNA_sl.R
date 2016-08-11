#WGCNA
#Note: WGCNA is not the right tool for 
#      analyzing gene networks OVER TIME
library("WGCNA")
library("DESeq2")
library("Rsamtools")

samples <- rownames(metaProxC[ metaProxC$Subgroup2 == "DG" & metaProxC$FOS != "L"  & metaProxC$Context1 == "none" & metaProxC$Subgroup2!= "HDG" & metaProxC$EE_ArcGroup != "Unk" & metaProxC$outliers == "in",])
dat <- tpmProxC[, samples]
rownames(dat) <- toupper(rownames(dat))
met <- metaProxC[match(samples,metaProxC$Sample_ID),]

m <- apply(X = dat,1,meanNoZero)
vstMat2 <- dat[m > 3,]
vstMat2 <- vstMat2[-c(grep("(^GM)",rownames(vstMat2), perl = TRUE)),]
##################################
#### Setup the WGCNA-Style Objects
##################################
#Define data set dimensions
datExpr <- as.data.frame(na.exclude((vstMat2)))
datTraits <- met[rownames(datExpr),]

##################################
#### Run WGCNA
##################################

options(stringsAsFactors = FALSE);
disableWGCNAThreads()
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
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
net = blockwiseModules(datExpr,maxBlockSize=500, power = POWER,
                       TOMType = "unsigned", minModuleSize = 2,
                       reassignThreshold = 0, #mergeCutHeight = 0.99,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       TOMDenom="min",
                       #saveTOMs = TRUE,
                       #saveTOMFileBase = "hnp",
                       verbose = 3)

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
MEs$Mouse_condition <- factor(datTraits$Mouse_condition,levels = c("HC","EE"))
MEs$FOS <- factor(datTraits$FOS, levels = c("N","F"))
df <- vector()
for(i in 1:(ncol(MEs)-2)){
  tmp <- data.frame(ME = MEs[,i], Mouse_condition = MEs$Mouse_condition,FOS = MEs$FOS)
  model <- summary(lm(ME ~ Mouse_condition*FOS, tmp))$coefficients
  df <- rbind(df, c(as.vector(as.matrix(model)), colnames(MEs)[i]))
}
colnames(df) <- c(paste(rep(c("Est","StdErr","t","p"), each = 3), rep(c("Int","Mouse_condition","FOS"),4),sep= "."), "ME")
df <- data.frame(df)
df2 <- as.data.frame(t(data.frame(apply(df[,c(1:12)], 1, as.numeric))))
df2$ME <- as.vector(df$ME)
colnames(df2) <- c(paste(rep(c("Est","StdErr","t","p"), each = 3), rep(c("Int","Mouse_condition","FOS"),4),sep= "."), "ME")

rownames(df2) <- df2$ME
#view this
sizeGrWindow(10,6)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
              # textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))


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
