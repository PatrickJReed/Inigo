## EDGER
library(pcaMethods)
library(edgeR)
library(scatterplot3d)
library(fdrtool)
library(ggplot2)
library(randomForest)
### other function
maximum <- function(i,g){
    return(sum(tpmProxC[names(Max[i]),colnames(dat)[group == g]] > Max[i]) / sum(group == g))
}
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
#NEWER (After CA3 tighter)
#save(list = c("RES.celltype","RES.activity","RES.FosF.HCN","RES.FosN.HCN","SigGenes"), file = "~/Documents/SalkProjects/ME/ShortLongSingature/MolecDissec_Figs_Tables/Figures_vC/edger.res",compress = TRUE)
#load(file = "~/Documents/SalkProjects/ME/ShortLongSingature/MolecDissec_Figs_Tables/Figures_vC/edger.res")
###
for (g2 in c("DG","CA1","VIP","CA3","Neg")){
  samples <- rownames(metaProxC[metaProxC$Mouse_condition == "EE" & metaProxC$Subgroup2 == g2 & metaProxC$FOS == "F" & metaProxC$Arc_2.5 != "greater" & metaProxC$outliers == "in"|
                                  metaProxC$Mouse_condition == "EE" &  metaProxC$Subgroup2 == g2 & metaProxC$FOS == "N" & metaProxC$Arc_2.5 != "greater"  & metaProxC$outliers == "in",])
  dat <- na.exclude(countProxC[, samples])
  dat <- dat[rowSums(dat) > 0,]
  met <- metaProxC[match(samples,metaProxC$Sample_ID),]
  ###################
  #Assign groups
  ###################
  group <-  met$FOS == "F"
  Pair <- levels(as.factor(as.character(group)))
  ###################
  # Test genes
  ###################
  #RES <- list()
  res <- exact(dat, group, Pair)
  i <- i + 1
  RES[[i]] <- res
  names(RES)[i] <- g2
}

dg <- RES[["DG"]]
ca1 <- RES[["CA1"]]
ca3 <- RES[["CA3"]]
vip <- RES[["VIP"]]
neg <- RES[["Neg"]]

#all genes
genes <- unique(c(rownames(dg[dg$f < 0.05 & dg$a > 0.2 | dg$f < 0.05 & dg$b > 0.2, ]),
                rownames(ca1[ca1$f < 0.05 & ca1$a > 0.2| ca1$f < 0.05 & ca1$b > 0.2,]),
                rownames(ca3[ca3$f < 0.05 & ca3$a > 0.2| ca3$f < 0.05 & ca3$b > 0.2,]),
                rownames(vip[vip$f < 0.05 & vip$a > 0.2| vip$f < 0.05 & vip$b > 0.2,]),
                rownames(neg[neg$f < 0.05 & neg$a > 0.2| neg$f < 0.05 & neg$b > 0.2,])))
#all shared genes
genes <- table(c(rownames(dg[dg$f < 0.01 & dg$a > 0.2 | dg$f < 0.01 & dg$b > 0.2, ]),
                  rownames(ca1[ca1$f < 0.01 & ca1$a > 0.2| ca1$f < 0.01 & ca1$b > 0.2,]),
                  rownames(ca3[ca3$f < 0.01 & ca3$a > 0.2| ca3$f < 0.01 & ca3$b > 0.2,]),
                  rownames(vip[vip$f < 0.01 & vip$a > 0.2| vip$f < 0.01 & vip$b > 0.2,]),
                  rownames(neg[neg$f < 0.01 & neg$a > 0.2| neg$f < 0.01 & neg$b > 0.2,])))


#############
# Random Forest
#############
genes2 <- genes[-c(grep("Rik",genes))]
genes2 <- genes2[-c(grep("-",genes2))]
dat <- na.exclude(tpmProxC[, samples])
met <- metaProxC[match(samples,metaProxC$Sample_ID),]
#

dat2 <- as.data.frame(t(dat[genes2, ]),strings.as.factors = FALSE)
dat2$Subgroup <- met$Subgroup2
dat2$Subgroup <- as.factor(dat2$Subgroup)
rf.model <- randomForest(Subgroup ~ . , dat2)
gini <- rf.model$importance
gini <- gini[order(gini,decreasing=TRUE),]

