## EDGER
library(edgeR)
library(fdrtool)
########################
###Load and Format Data
########################

GLM <- function(dat, variable1, variable2 = NULL, prefit=TRUE){
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
  fit <- glmQLFit(y=cds,design)
  lrt <- glmLRT(fit, coef=Coef)#interaction term = 4
  edg <- data.frame(lrt$table)
  edg <- edg[order(edg$PValue),]
  f <- fdrtool(x=edg$PValue,statistic="pvalue",plot=FALSE)
  edg$f <- f$qval
  return(edg)
}


dat <- na.exclude(countProxC[,metaProxC$prox == "P" & metaProxC$cond == "EE"])
meta <- metaProxC[metaProxC$prox == "P" & metaProxC$cond == "EE",]
group <- as.factor(meta$fos)
cds <- DGEList(dat, group= group)

########################
### Pre-process Data
########################
# Filter out really lowly detected reads (*from EdgeR tutorial)
#cds <- cds[rowSums(1e+06 * cds$counts/expandAsMatrix(cds$samples$lib.size, dim(cds)) > 1) >= 3, ]
cds <- calcNormFactors( cds )
########################
### Estimate parameters from data
########################
cds <- estimateCommonDisp( cds )
cds <- estimateTagwiseDisp( cds )

de.cmn <- exactTest( cds , pair = c( "N" , "F" ) )
res <- as.data.frame(de.cmn$table)
f <- fdrtool(x=res$PValue,statistic="pvalue",plot=FALSE)
res$f <- f$qval
res <- res[order(res$PValue),]
#edgeR_P_N_fospos <- res
#edgeR_P_N_fosn <- res
edgeR_F_N_prox <- res


a <- table(c(rownames(edgeR_P_N_fospos[edgeR_P_N_fospos$f < 0.05,]),
              rownames(edgeR_P_N_fospos[edgeR_P_N_fospos$f < 0.05,]),
              rownames(edgeR_P_N_fosn[edgeR_P_N_fosn$f < 0.05,])))


edgeR_P_N_fospos_min <- edgeR_P_N_fospos[names(a[a==2]),]
edgeR_P_N_fospos_min <- edgeR_P_N_fospos_min[order(edgeR_P_N_fospos_min$PValue),]

edgeR_P_N_fospos_min_prox <- na.exclude(edgeR_P_N_fospos_min[rownames(edgeR_F_N_prox[edgeR_F_N_prox$f < 0.05,]),])
edgeR_P_N_fospos_min_prox <- edgeR_P_N_fospos_min_prox[order(edgeR_P_N_fospos_min_prox$PValue),]
edgeR_P_N_fospos_min_prox$proxp <- edgeR_P_N_fospos_min_prox[order(edgeR_P_N_fospos_min_prox$PValue),"PValue"]
edgeR_P_N_fospos_min_prox <- edgeR_P_N_fospos_min_prox[order(edgeR_P_N_fospos_min_prox$proxp),]
########################
### Interaction Term between Prox and Fos
########################

dat <- na.exclude(countProxC[,metaProxC$fos != "L"])
meta <- metaProxC[metaProxC$fos != "L",]
variable1 <- as.factor(meta$prox)
variable2 <- as.factor(as.character(meta$fos))

modelME <- function(x){
  require(lme4)
  m <- summary(lmer(value ~ prox + (1|fos), data=tmp[tmp$X2 == x,]))
  return(m$coefficients[2,])
}

tmp <- melt(t(dat))
tmp <- cbind(tmp, meta)
save(list=c("dat","tmp"),file="~/Documents/test.rda")

p <- apply(X=dat,MARGIN=1, FUN=propExp)
#
cl = makeCluster(getOption("cl.cores",3))
clusterExport(cl=cl,varlist=c("tmp"))
bleg <- parLapply(cl=cl,X=rownames(dat),fun=modelME)
res_p <- as.data.frame(do.call("rbind",bleg))
rownames(res_p) <- rownames(dat)
colnames(res_p) <- c("Estimate","stderr","t")
stopCluster(cl=cl)