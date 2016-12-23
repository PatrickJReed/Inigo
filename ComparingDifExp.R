library(scde)
library(scDD)
library(Biobase)
library(edgeR)
library(SummarizedExperiment)
require(biomaRt)
require(plyr)

#save(list = c('RES.activity_meta'),file =  "~/Documents/SalkProjects/ME/ShortLongSingature/SLSig_R/meta.rda",compress = TRUE)
#load("~/Documents/SalkProjects/ME/ShortLongSingature/SLSig_R/meta.rda")
#only need mart if ranges = TRUE, and it's not functional right now as such
#mart <- useEnsembl(dataset = "mmusculus_gene_ensembl", biomart = "ensembl", version = 83)

GLM <- function(dat, variable1, variable2 = NULL, prefit=FALSE){
  if(is.null(variable2)){
    design <- model.matrix(~variable1)
    Coef = 2
  }else{
    design <- model.matrix(~variable2*variable1)
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
  Max_a <- apply(tpmProxC[rownames(res),colnames(dat)[group == TRUE]],1,max)
  Max_b <- apply(tpmProxC[rownames(res),colnames(dat)[group == FALSE]],1,max)
  tmp_a <- cbind(dat[rownames(res),colnames(dat)[group == TRUE]], Max_b[rownames(res)])
  Col <<- sum(group == TRUE)
  res$max_a <- as.numeric(apply(tmp_a,1,maximum2))
  tmp_b <- cbind(dat[rownames(res),colnames(dat)[group == FALSE]], Max_a[rownames(res)])
  Col <<- sum(group == FALSE)
  res$max_b <- as.numeric(apply(tmp_b,1,maximum2))
  
  return(res)
}
SCDE <- function(dat, group){
  ##############################################
  ## Step 1: Load in Count Table
  ##############################################
  countTable2 <- t(apply(round(dat[rowSums(dat) > 0,],0),1, as.integer))
  colnames(countTable2) <- colnames(dat)
  ##############################################
  ## Step 2: Fit error model
  ##############################################
  
  # number of local process cores to use in processing
  # For goodness sake don't go above 2 or your computer will die...blegh
  n.cores <- 2
  # calculate models
  o.ifm <- scde.error.models(counts=countTable2,
                             groups=group,
                             n.cores=n.cores,
                             threshold.segmentation=T,
                             save.crossfit.plots=F,
                             save.model.plots=F,
                             verbose=1)
  # Filter out cells with really bad fits
  valid.cells <- o.ifm$corr.a >0
  table(valid.cells)
  # estimate gene expression prior
  o.prior <- scde.expression.prior(models=o.ifm,
                                   counts=countTable2,
                                   length.out=400,
                                   show.plot=F
                                   )
  
  ##############################################
  ## Step 3: Test Differential Expression
  ##############################################
  # define two groups of cells
  groups <- as.factor(group)
  names(groups) <- row.names(o.ifm);
  # run differential expression tests on all genes.
  ediff <- scde.expression.difference(o.ifm,
                                      countTable2,
                                      o.prior,
                                      groups=groups,
                                      n.randomizations=100,
                                      n.cores=n.cores,
                                      verbose=1)
  ediff2 <- ediff[order(abs(ediff$Z),decreasing=T),]
  ediff2$f <- 2*pnorm(-abs(ediff2$cZ))
  a <- apply(tpmProxC[,colnames(dat)[group == TRUE]],1,rawExp,1)
  b <- apply(tpmProxC[,colnames(dat)[group == FALSE]],1,rawExp,1)
  ediff2$a <- a[rownames(ediff2)] / sum(group == TRUE)
  ediff2$b <- b[rownames(ediff2)] / sum(group == FALSE)
 return(ediff2)
}
SCDD <- function(dat, group, samples, ranged = FALSE){
  countTable2 <- t(apply(round(dat[rowSums(dat) > 0,],0),1, as.integer))
  colnames(countTable2) <- colnames(dat)
  cts <- as.matrix(countTable2[,samples])
  pheno <- data.frame(condition = group,row.names = samples)
  if (ranged == TRUE){# Add row meta data into the eset object
    genes <- as.character(rownames(countTable2))
    convert <- getBM(
      filters = "external_gene_name",
      attributes = c("external_gene_name","chromosome_name","start_position","end_position","strand"),#,"external_gene_name"),
      values = c(genes),
      mart = mart,
      verbose = FALSE
    )
    
    convert2 <- convert[match(genes, convert$external_gene_name),]
    rowRanges <- GRanges(convert2$chromosome_name,
                        IRanges(convert2$start_position, width=(convert2$end_position - convert2$start_position)),
                        strand=convert2$strand,
                        feature_id=convert2$external_gene_name)
    eset.scdd <- SummarizedExperiment(assays = list(cts = cts), 
                                      rowRanges = rowRanges,
                               colData = DataFrame(pheno))
  }else{
    eset.scdd <- SummarizedExperiment(assays = list(cts = cts), 
                                      colData = DataFrame(pheno))
  }
  
  #Create list of prior parameters for model fitting, why are these the defaults? Should play around with this
  prior_param=list(alpha=0.01, mu0=0, s0=0.01, a0=0.01, b0=0.01)
  
  #Detect and classify differentially distributed genes
  set.seed(6767) # set random seed for reprodicbility
  dd.results <- scDD(eset.scdd, prior_param=prior_param, permutations=0, testZeroes=TRUE)
  
  res <- as.data.frame(dd.results@metadata)[,c(1:9)]
  res <- res[order(res$Genes.nonzero.pvalue.adj),]
  rownames(res) <- res$Genes.gene
  a <- apply(tpmProxC[,colnames(dat)[group == TRUE]],1,rawExp,1)
  b <- apply(tpmProxC[,colnames(dat)[group == FALSE]],1,rawExp,1)
  res$a <- a[rownames(res)] / sum(group == TRUE)
  res$b <- b[rownames(res)] / sum(group == FALSE)
  
  return(res)
  
}
maximum <- function(i,g){
  return(sum(tpmProxC[names(Max[i]),colnames(dat)[group == g]] > Max[i]) / sum(group == g))
}
maximum2 <- function(x){
  return(sum(x[1:Col] > x[length(x)]) / Col)
}
fisher.meta <- function(Ps){
  Ps <- Ps[!is.na(Ps)]
  pchisq(-2 * sum(log(Ps)), df = 2 * length(Ps),lower.tail = FALSE)
}
rawExp <- function(x,i=2){
  sum(na.exclude(x) > i)
}
###############
## Run loop over celltypes
###############
RES.activity_meta <- list()
i <- 0
#for (g2 in c("DG","CA1","VIP")){
for (g2 in c("DG","CA1","VIP")){
  i <- i  + 1
  samples <- rownames(metaProxC[metaProxC$predicted == g2 & metaProxC$Mouse_condition == "EE" & metaProxC$FOS == "N" & metaProxC$cluster_outlier == "in" & metaProxC$outliers == "in" & metaProxC$Arc_2.5 != "greater"  |
                                  metaProxC$predicted == g2 & metaProxC$Mouse_condition == "EE" & metaProxC$FOS == "F" & metaProxC$cluster_outlier == "in" & metaProxC$outliers == "in" & metaProxC$Arc_2.5 != "greater"  ,])#
  dat <- na.exclude(countProxC[, samples])
  dat <- dat[rowSums(dat) > 0,]
  met <- metaProxC[samples,]
  group = as.factor(met$FOS == "F")
  Pair = levels(group)
  ###############
  ## EDGER EXACT
  ###############
  res.exact <- exact(dat, group, Pair)
  ###############
  ## SCDD
  ###############
  res.scdd <- SCDD(dat, group, samples)
  ###############
  ## SCDE
  ###############
  res.scde <- SCDE(dat, group)
  ###############
  # Meta analysis
  ###############
 
  df2 <- data.frame(p1 = res.exact[rownames(dat),"f"], 
               p2 = res.scdd[rownames(dat), "Genes.nonzero.pvalue.adj"],
               p3 = res.scdd[rownames(dat), "Genes.zero.pvalue.adj"],
               p4 = res.scde[rownames(dat),"f"],
               row.names = rownames(dat)
               )
  
  
  all  <- apply(X = df2, MARGIN = 1, FUN = fisher.meta)
  res.comb <- data.frame(res.exact[rownames(dat),], res.scdd[rownames(dat),], res.scde[rownames(dat),],meta_padj = all)
  res.comb <- res.comb[order(res.comb$meta_padj),]
  RES.activity_meta[[i]] <- res.comb
  names(RES.activity_meta)[i] <- g2
}


#################
## Number of genes
#################

#A <- unique(c(rownames(res.exact[(res.exact$f < 0.05 & res.exact$logFC > 0 & res.exact$a > 0.2),]) ,  rownames(res.exact[(res.exact$f < 0.05 & res.exact$logFC < 0 & res.exact$b > 0.2),])))
A <- unique(c(rownames(res.exact[(res.exact$f < 0.05 & res.exact$logFC > 0 ),]) ,  rownames(res.exact[(res.exact$f < 0.05 & res.exact$logFC < 0 ),])))
B <- rownames(res.scde[(res.scde$f < 0.05),])
C.1 <- res.scdd[!is.na(res.scdd$Genes.nonzero.pvalue.adj),]
C.2 <- res.scdd[!is.na(res.scdd$Genes.zero.pvalue.adj),]
C <- unique(c(rownames(C.1[C.1$Genes.nonzero.pvalue.adj < 0.05 ,]), rownames(C.2[C.2$Genes.zero.pvalue.adj < 0.05 ,])))

nm <- table(c(A,A,A,A,B,B,C))
table(nm)
genes <-names(nm[nm == 2])
a <- res.scde[genes,]
a <- a[order(a$f),]
tail(a)
#4 = exact only
#5 = exact and scdd
#6 = exact and scde
#3 = scdd and scde
#2 = scde only
#1 = scdd only