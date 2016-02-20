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
  cds <- DGEList(dat)
  
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
groupname <- function(x){
  return(paste(x[1],x[2],sep="."))
}

metaProxC <- metaProxC[colnames(tpmProxC),]
dat <- tpmProxC[,metaProxC$Mouse_condition == "HC" & metaProxC$Smartseq2_RT_enzyme_used == "Protoscript_II"]
meta <- metaProxC[metaProxC$Mouse_condition == "HC" & metaProxC$Smartseq2_RT_enzyme_used == "Protoscript_II",]
group <- as.factor(apply(cbind(as.character(meta$PROX1),as.character(meta$CTIP2)),1,groupname))


hc_nuclei_edger <- as.data.frame(GLM(dat, group))
hc_nuclei_edger <- hc_nuclei_edger[order(hc_nuclei_edger$f_N.N),]
