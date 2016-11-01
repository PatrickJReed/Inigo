#save(list = c("sample.colors"),file = "~/Documents/SalkProjects/ME/ShortLongSingature/neg.rda",compress = TRUE)
#load("~/Documents/SalkProjects/ME/ShortLongSingature/neg.rda")

dat <- na.exclude(tpmProxC[, samples])
met <- metaProxC[samples,]
group <-  group#as.factor(ifelse(as.numeric(tpmProxC["Vip",samples])  >7, yes = "high","low")) 


#####
g <- "Htr2c"
tmp <- data.frame(tpm = as.numeric(dat[g,]),
                  group = group)

ggplot(tmp, aes(group, tpm, fill = group, alpha = 0.3))+
  geom_violin()+
  geom_jitter()+
  labs(title = g)

###
Cor <- cor(t(tpmProxC[,samples]))
a <- Cor["Cntn2",]
a <- a[order(a)]

#res.vip_vip <- res
#cor.vip <- Cor

#########
samples <- rownames(t[t$glut_gaba == "GABA" & t$Brain_Region != "P+C-",])
dat <- tpmProxC[,samples]
met <- metaProxC[samples,]
colnames(dat) <- paste(met$FOS,as.factor(ifelse(as.numeric(tpmProxC["Vip",samples])  >7, yes = "high","low")) , variable1 , sep = ".")
h <- hclust(dist(t(dat[genes,])))
variable1 <- k <- as.numeric(cutree(h, k = 3))
plot(h)

###GLM
if(is.null(variable2)){
  design <- model.matrix(~factor(variable1))
  dat <- na.exclude(countProxC[, samples])
  dat <- dat[rowSums(dat) > 0,]
  cds <- DGEList(dat)
  
    cds <- calcNormFactors(cds)
    cds <- estimateGLMCommonDisp( cds )
    cds <- estimateGLMTrendedDisp(cds)
    fit <- glmQLFit(y=cds,design)

  lrt <- glmLRT(fit, coef=2:3)#interaction term = 4
  edg <- data.frame(lrt$table)
  edg <- edg[order(edg$PValue),]
  f <- p.adjust(edg$PValue)
  edg$f <- f
  a <- apply(tpmProxC[,colnames(dat)], 1, propExp)
  edg$propExp <- as.vector(a[rownames(edg)])
  return(edg)
}

