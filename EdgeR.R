## EDGER

library(edgeR)
library(fdrtool)

#save(list=c("res.HC_N_P_1","res.HC_P_sub","res_CA1_Fos"),file="~/Documents/SalkProjects/ME/ShortLongSingature/SLSig_resTables/edgeR.res",compress=TRUE)
#load("~/Documents/SalkProjects/ME/ShortLongSingature/SLSig_resTables/edgeR.res")
########################
###Load and Format Data
########################
dat <- countProx[,metaProx$cond == "HC" & metaProx$date != "150629"]
group <- metaProx[metaProx$cond == "HC" & metaProx$date != "150629", "prox"]
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

de.cmn <- exactTest( cds , pair = c( "N" , "P" ) )
res <- as.data.frame(de.cmn$table)
res <- res[order(res$PValue),]
f <- fdrtool(x=res$PValue,statistic="pvalue")
res$f <- f$qval
res.HC_N_P_1 <- res

write.table(x=rownames(res.HC_N_P_1[res.HC_N_P_1$f < 0.05 & res.HC_N_P_1$logFC > 0,]),file="~/Documents/test.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)
