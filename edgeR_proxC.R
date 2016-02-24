library(edgeR)
GLM <- function(dat, variable1, variable2 = NULL){
  require(fdrtool)
  if(is.null(variable2)){
    design <- model.matrix(~variable1)
    Coef = 2
  }else{
    design <- model.matrix(~variable2*variable1)
    Coef = 4
  }
    cds <- DGEList(dat, group = group)
    cds <- calcNormFactors(cds)
    cds <- estimateGLMCommonDisp( cds )
    cds <- estimateGLMTrendedDisp(cds)
    fit <- glmQLFit(y=cds,design)
  
  all <- vector()
  for(i in 2:ncol(fit$coefficients)){
    lrt <- glmLRT(fit, coef=i)
    edg <- data.frame(lrt$table)
    f <- fdrtool(x=edg$PValue,statistic="pvalue",plot=FALSE)
    edg$f <- f$qval
    colnames(edg) <- paste(colnames(edg),
                           strsplit(colnames(fit$coefficients)[i],split = "variable1",fixed = TRUE)[[1]][2],sep="_")
    all <- cbind(all,as.matrix(edg))
  }
  return(all)
}
getFit.glm <- function(dat, variable1, variable2){
  if(is.null(variable2)){
    design <- model.matrix(~variable1)
    Coef = 2
  }else{
    design <- model.matrix(~variable2*variable1)
    Coef = 4
  }
  cds <- DGEList(dat)
  cds <- calcNormFactors(cds)
  cds <- estimateGLMCommonDisp( cds )
  cds <- estimateGLMTrendedDisp(cds)
  cds <- estimateGLMTagwiseDisp( cds )
  fit <- glmQLFit(y=cds,design)
  return(fit$fitted.values)
}
EXACT <- function(dat, group){
  cds <- DGEList(dat, group=group)
  cds <- calcNormFactors( cds )
  cds <- estimateCommonDisp( cds )
  cds <- estimateTagwiseDisp( cds )
  
  de.cmn <- exactTest( cds , pair = levels(as.factor(group)) )
  res <- as.data.frame(de.cmn$table)
  f <- fdrtool(x=res$PValue,statistic="pvalue",plot=FALSE)
  res$f <- f$qval
  res <- res[order(res$PValue),]
  return(res)
}
groupname <- function(x){
  return(paste(x[1],x[2],sep="."))
}

metaProxC <- metaProxC[colnames(tpmProxC),]
dat <- tpmProxC[,metaProxC$Mouse_condition == "HC" ]
meta <- metaProxC[metaProxC$Mouse_condition == "HC" ,]
group <- as.factor(apply(cbind(as.character(meta$PROX1),as.character(meta$CTIP2)),1,groupname))

#All cell types
hc_nuclei_edger <- as.data.frame(GLM(dat, group))
hc_nuclei_edger <- hc_nuclei_edger[order(hc_nuclei_edger$f_N.C),]

#DG versus hDGC
dat <- tpmProxC[,metaProxC$Mouse_condition == "HC" & metaProxC$PROX1 == "P" ]
meta <- metaProxC[metaProxC$Mouse_condition == "HC" & metaProxC$PROX1 == "P" ,]
group <- meta$CTIP2
hDG <- as.data.frame(EXACT(dat, group))

#CA1 specific
dat <- tpmProxC[,metaProxC$Mouse_condition == "HC" ]
meta <- metaProxC[metaProxC$Mouse_condition == "HC",]
group <- meta$Brain_Region == "CA1"
CA1 <- as.data.frame(EXACT(dat, group))
