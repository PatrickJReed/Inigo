##########
## Error Rate
##########
library(pcaMethods)
library(edgeR)
library(scatterplot3d)
library(fdrtool)
library(ggplot2)
### other functions
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

#Get Samples
samples <- metaProxC[ metaProxC$Mouse_condition == "HC" & metaProxC$Subgroup != "Unk" & metaProxC$FOS == "N" & metaProxC$Context1 == "none" & metaProxC$outliers == "in" ,"Sample_ID"]#
dat <- na.exclude(tpmProxC[, samples])
met <- metaProxC[match(samples,metaProxC$Sample_ID),]
