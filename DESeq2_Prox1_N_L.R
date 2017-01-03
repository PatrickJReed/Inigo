####################################################
## Title: DESeq2
## Created: 01/26/2015, Slinker
## Purpose: A script to run the basic functions of Deseq, specially written for Kumar :D
## Required Input:
##   1) countTable: 
##        class= data.frame of whole number counts (Deseq does not take fractions, e.g., RSEM output)
##        columns= samples
##        rows= genes
##   2) label:
##        class= vector of factors
##        info= case/control labels for each sample. Ex: label = "case", "case", "control" ...
##   3) normal: (bad name I know)
##        class= character vector of length 1
##        info= the label that equations will be referenced to. Ex: normal = "control"
##
####################################################
load(file="~/Documents/SalkProjects/ME/ShortLongSingature/SLSig_R/DESeq2_Prox1_N_L.rda")
####################################################
## Step 1) Load Libraries
####################################################

library("GenomicRanges")
library("Rsamtools")
library("GenomicFeatures")
library("rtracklayer")
library("DESeq2")

####################################################
## Step 2) Define User Specified Variables
####################################################
condition <- label <- as.factor(metaProx$fos)

normal <- as.character("N")
countTable2 <- countProx
countTable3 <- countTable2[rowSums(countTable2) > 0,]

####################################################
## Step 3) Setup the data.frames required by DESeq2
##    Note: You shouldn't have to change anything by hand from this point on
####################################################
#Calculate the replicate order for deseq
rep.order <- vector()
factor1 <- 0
factor2 <- 0
l <- levels(label)
for (i in as.factor(label)){
  if (i == l[1]){
    factor1 <- factor1 + 1
    rep.order <- c(rep.order,factor1)
  }else{
    factor2 <- factor2 + 1
    rep.order <- c(rep.order,factor2)
  }
}

#prepare the dataframe for DESeq2
colData <- data.frame(condition,c(rep("paired-end",length(rep.order))))
colnames(colData) <- c("condition","type")
dds <- DESeqDataSetFromMatrix(countData = countTable3,
                              colData = colData,
                              design = ~ condition)

dds$condition <- relevel(dds$condition, normal)

####################################################
## Step 4) Run DESeq2 
####################################################
dds <- tryCatch(DESeq(dds),error=function(e) NULL)

####################################################
## Step 5) Get the Results!! 
####################################################
res <- results(dds)
#put them in order of signficance
resOrdered <- as.data.frame(res[order(res$padj),])
head(resOrdered)

# Most recently calc on hcluster k =2 all genes
save(list=c("resOrdered"),file="~/Documents/SalkProjects/ME/ShortLongSingature/SLSig_R/DESeq2_Prox1_N_L.rda",compress=TRUE)

####################################################
## Extra Goodies) Get normalized data for downstream analyses
####################################################
#Rlog transformation with no prior assumption of groups
rld <- rlogTransformation(dds,blind=TRUE)
rlogMat <- assay(rld)
#rlogMat <- rld
#variance stabalized data
vsd <- varianceStabilizingTransformation(dds)
vstMat <- assay(vsd)
