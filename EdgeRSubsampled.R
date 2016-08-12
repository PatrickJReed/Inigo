#save(list = c("SubgroupNeg","SubgroupCA2"),file = "~/Documents/SalkProjects/ME/ShortLongSingature/SLSig_R/subsampledifexp.rda",compress = TRUE)
#load("~/Documents/SalkProjects/ME/ShortLongSingature/SLSig_R/subsampledifexp.rda")
####
## Reduced sample size difexp
####
SubsampleEdgeR <- function(samples.2, condition = "DG" , met){
  genes <- vector()
  for (i in 1:20){
    a <- is.na(match(rownames(met), samples.2))
    samples.1 <- rownames(met)[a]
    samples <- c(samples.2, sample(x = samples.1,size = length(samples.2),replace = FALSE))
    dat <- na.exclude(countProxC[, samples])
    dat <- dat[rowSums(dat) > 0,]
    met <- metaProxC[match(samples,metaProxC$Sample_ID),]
    ###################
    #Assign groups
    ###################
    group <- met$Subgroup == condition
    Pair <- levels(as.factor(as.character(group)))
    ###################
    # Test genes
    ###################
    #!!!Run exact test
    res <- exact(dat, group, Pair)
    genes <- c(genes, rownames(res[res$logFC > 0 & res$f < 0.05,]))
  }
  
  genes2 <- as.data.frame(table(genes))
  genes2 <- genes2[order(genes2$Freq,decreasing = TRUE),]
  return(genes2)
}
samples.2 <- rownames(df[df$group == "close",])

genes <- vector()
for (i in 1:20){
samples.1 <- rownames(df[df$group == "far",])
samples <- c(samples.2, sample(x = samples.1,size = length(samples.2),replace = FALSE))
dat <- na.exclude(countProxC[, samples])
dat <- dat[rowSums(dat) > 0,]
met <- metaProxC[match(samples,metaProxC$Sample_ID),]
###################
#Assign groups
###################
group <- rep(c("close","far"), each = length(samples.2))# & as.numeric(log(dat["Gad2",])) > 4
Pair <- levels(as.factor(as.character(group)))
###################
# Test genes
###################
#!!!Run exact test
res <- exact(dat, group, Pair)
genes <- c(genes, rownames(res[res$logFC > 0 & res$f < 0.05,]))
}

genes2 <- as.data.frame(table(genes))
genes2 <- genes2[order(genes2$Freq,decreasing = TRUE),]
#SubgroupNeg <- genes2
write.table(x = as.character(genes2[genes2$Freq > 12,"genes"]),file = "~/Documents/test.txt",quote = FALSE,row.names = FALSE,col.names = FALSE)

#dg.genes2 <- genes2