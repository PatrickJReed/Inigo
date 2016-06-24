#save(list = c("SubgroupNeg"),file = "~/Documents/SalkProjects/ME/ShortLongSingature/SLSig_R/subsampledifexp.rda",compress = TRUE)
#load("~/Documents/SalkProjects/ME/ShortLongSingature/SLSig_R/subsampledifexp.rda")
####
## Reduced sample size difexp
####

samples.2 <- metaProxC[metaProxC$Subgroup == "Neg" & metaProxC$outlier == "in" & metaProxC$Mouse_condition == "HC" & metaProxC$FOS == "N" &metaProxC$Context1 == "none" & metaProxC$alignable >  100000 & metaProxC$Smartseq2_RT_enzyme_used == "ProtoscriptII" ,"Sample_ID"]

genes <- vector()
for (i in 1:100){
samples.1 <- metaProxC[metaProxC$Subgroup != "VIP" & metaProxC$Subgroup != "Erbb4+" & metaProxC$Subgroup != "Neg" & metaProxC$outlier == "in" & metaProxC$Mouse_condition == "HC" & metaProxC$FOS == "N" &metaProxC$Context1 == "none" & metaProxC$alignable >  100000 & metaProxC$Smartseq2_RT_enzyme_used == "ProtoscriptII" ,"Sample_ID"]
samples <- c(samples.2, sample(x = samples.1,size = length(samples.2),replace = FALSE))
dat <- na.exclude(countProxC[, samples])
dat <- dat[rowSums(dat) > 0,]
met <- metaProxC[match(samples,metaProxC$Sample_ID),]
###################
#Assign groups
###################
group <- met$Subgroup == "Neg"# & as.numeric(log(dat["Gad2",])) > 4
Pair <- levels(as.factor(as.character(group)))
###################
# Test genes
###################
#!!!Run exact test
res <- exact(dat, group, Pair)
genes <- c(genes, rownames(res[res$logFC > 0 & res$f < 0.05,]))
}

genes2 <- as.data.frame(table(genes))
#SubgroupNeg <- genes2