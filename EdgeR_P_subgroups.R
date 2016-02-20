## EDGER

library(edgeR)
library(fdrtool)

#save(list="res.HC_P_sub",file="~/Documents/SalkProjects/ME/ShortLongSingature/SLSig_resTables/edgeR_P_group.res",compress=TRUE)
#load("~/Documents/SalkProjects/ME/ShortLongSingature/SLSig_resTables/edgeR_P_group.res")
########################
###Load and Format Data
########################
dat <- countProx[,metaProx$prox == "N"]
group <- metaProx[metaProx$prox == "N","fos"]#tempGroupingProx1HC
group[group == "TRUE"] <- "A"
group[group == "FALSE"] <- "B"
cds <- DGEList(dat, group= group)

########################
### Pre-process Data
########################
# Filter out really lowly detected reads (*from EdgeR tutorial)
cds <- cds[rowSums(1e+06 * cds$counts/expandAsMatrix(cds$samples$lib.size, dim(cds)) > 1) >= 3, ]
cds <- calcNormFactors( cds )
########################
### Estimate parameters from data
########################
cds <- estimateCommonDisp( cds )
cds <- estimateTagwiseDisp( cds )

de.cmn <- exactTest( cds , pair = c( "N" , "F" ) )
res <- as.data.frame(de.cmn$table)
res <- res[order(res$PValue),]
f <- fdrtool(x=res$PValue,statistic="pvalue")
res$f <- f$qval
#res.HC_P_sub <- res

