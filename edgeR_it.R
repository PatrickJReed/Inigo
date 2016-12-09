#################
## Iterative edgeR
#################
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
maximum2 <- function(x){
  return(sum(x[1:Col] > x[length(x)]) / Col)
}
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
N <- 15
it <- 20

df <- matrix(0,nrow = nrow(countProxC),ncol = 3)
rownames(df) <- rownames(countProxC)
j <- 0
for(g2 in c("DG","VIP","CA1")){
  j <- j+1
  genes <- vector()
  for (i in 1:it){
    samples.f <- rownames(metaProxC[metaProxC$Subgroup == g2 & metaProxC$Mouse_condition == "EE" & metaProxC$Context1 == "none" & metaProxC$FOS == "F" & metaProxC$cluster_outlier == "in" & metaProxC$outliers == "in" & metaProxC$Arc_2.5 != "greater",])
    samples.n <- rownames(metaProxC[metaProxC$Subgroup == g2 & metaProxC$Mouse_condition == "EE" & metaProxC$Context1 == "none" & metaProxC$FOS == "N" & metaProxC$cluster_outlier == "in" & metaProxC$outliers == "in"& metaProxC$Arc_2.5 != "greater",])
    samples <- c(sample(x = samples.f, size = N, replace = FALSE),
                 sample(x = samples.n, size = N, replace = FALSE)
                 )
    dat <- na.exclude(countProxC[, samples])
    dat <- dat[rowSums(dat) > 0,]
    met <- metaProxC[samples,]
    group <- factor(met$FOS)
    Pair <- levels(group)
    res <- exact(dat, group, Pair)
    genes <- c(genes,rownames(res[res$f < 0.05,]))
  }
  a <- table(genes)
  df[names(a),j] <- as.numeric(a)
}

df <- data.frame(df)
colnames(df) <- c("DG","VIP","CA1")
df2 <- melt(t(df))
df2$X1 <- factor(df2$X1, c("DG","CA1","VIP"))
ggplot(df2[df2$value > 5,], aes(value > 5, fill = X1))+
  geom_bar(position = "dodge")

sum(df2[df2$X1 == "DG","value"] > 0)
sum(df2[df2$X1 == "CA1","value"]> 0)
sum(df2[df2$X1 == "VIP","value"] > 0)

