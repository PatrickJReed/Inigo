## EDGER
library(pcaMethods)
library(edgeR)
library(scatterplot3d)
library(fdrtool)
library(ggplot2)
### other function
GLM <- function(dat, variable1, variable2 = NULL, prefit=FALSE){
  if(is.null(variable2)){
    design <- model.matrix(~variable1)
    Coef = 2
  }else{
    design <- model.matrix(~variable2+variable1)
    Coef = 4
  }
  cds <- DGEList(dat)
  if (prefit == FALSE){
    cds <- calcNormFactors(cds)
    cds <- estimateGLMCommonDisp( cds )
    cds <- estimateGLMTrendedDisp(cds)
    fit <- glmQLFit(y=cds,design)
  }
  #fit <- glmQLFit(y=cds,design)
  lrt <- glmLRT(fit, coef=Coef)#interaction term = 4
  edg <- data.frame(lrt$table)
  edg <- edg[order(edg$PValue),]
  f <- p.adjust(edg$PValue)
  edg$f <- f
  a <- apply(tpmProxC[,colnames(dat)], 1, propExp)
  edg$propExp <- as.vector(a[rownames(edg)])
  return(edg)
}
exact <- function(dat, group, Pair){
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
  a <- apply(tpmProxC[,colnames(dat)[group == TRUE]],1,rawExp,1)
  b <- apply(tpmProxC[,colnames(dat)[group == FALSE]],1,rawExp,1)
  res$a <- a[rownames(res)] / sum(group == TRUE)
  res$b <- b[rownames(res)] / sum(group == FALSE)
  return(res)
}

########################
###Load and Format Data
########################
#save(list = c("celltypeorder.dg", "celltypeorder.pin", "celltypeorder.ca1", "celltypeorder.neg","celltypeorder","activitygenes","celltypegenes","celltypegenes.hdg", "celltypegenes.dg", "celltypegenes.ca1", "celltypegenes.neg","celltypegenes.ca23","celltypegenes.in","activitygenes.ca1","activitygenes.dg","activitygenes.hdg","activitygenes.neg","RES","RES2"),file = "~/Documents/SalkProjects/ME/ShortLongSingature/SLSig_R/edgeR_slsig.rda",compress = TRUE)
#load("~/Documents/SalkProjects/ME/ShortLongSingature/SLSig_R/edgeR_slsig.rda")
###
samples <- rownames(metaProxC[metaProxC$Mouse_condition == "HC" & metaProxC$FOS == "N"  & metaProxC$Context1 == "none" & metaProxC$Subgroup2!= "HDG" & metaProxC$outliers == "in",])
dat <- na.exclude(countProxC[, samples])
dat <- dat[rowSums(dat) > 0,]
met <- metaProxC[match(samples,metaProxC$Sample_ID),]
###################
#Assign groups
###################
group <- met$Subgroup2 == "VIP"
Pair <- levels(as.factor(as.character(group)))
###################
# Test genes
###################
#!!!Run exact test
res <- exact(dat, group, Pair)
############################
### Keep track of gene lists 
############################
variable1 <- pData(my.data5)$Pseudotime
res <- GLM(dat, variable1)
#######
S_matrix <- reducedDimS(my.data5)
lib_info_with_pseudo <- pData(my.data5)
ica_space_df <- data.frame(t(S_matrix[c(x, y), ]))
colnames(ica_space_df) <- c("ICA_dim_1", "ICA_dim_2")
ica_space_df$sample_name <- row.names(ica_space_df)
df <- as.data.frame(merge(ica_space_df, lib_info_with_pseudo, 
                          by.x = "sample_name", by.y = "row.names"))
variable1 <- df$ICA_dim_2   
res <- GLM(dat, variable1)
