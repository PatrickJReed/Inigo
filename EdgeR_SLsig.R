## EDGER
library(pcaMethods)
library(edgeR)
library(scatterplot3d)
library(fdrtool)
library(ggplot2)
########################
###Load and Format Data
########################
#save(list = c("res.HC_P_NvC_N","res.HC_PvsN_NvsC_N","res.HC_PvsN_N_N","res.HC_PvN_C_N","res.HC_PvsN_CvsN_N","res.HC_N_CvsN_N","res.NE_P_N_FvN","res.NE_P_C_FvN","res.NE_N_C_FvN"),file = "~/Documents/SalkProjects/ME/ShortLongSingature/SLSig_R/edgeR_slsig.rda",compress = TRUE)
#load("~/Documents/SalkProjects/ME/ShortLongSingature/SLSig_R/edgeR_slsig.rda")
###
samples <- metaProxC[metaProxC$FOS != "L" & metaProxC$PROX1 == "N" & metaProxC$CTIP2 == "C" &  metaProxC$Mouse_condition == "EE" & metaProxC$alignable >  500000 & metaProxC$Smartseq2_RT_enzyme_used == "ProtoscriptII" ,"Sample_ID"]#
                       #metaProxC$CTIP2 == "N" & metaProxC$PROX1 == "N" & metaProxC$FOS == "N" & metaProxC$Mouse_condition == "HC" & metaProxC$alignable >  500000 & metaProxC$Smartseq2_RT_enzyme_used == "ProtoscriptII"  ,"Sample_ID"]
dat <- tpmProxC[, samples]
met <- metaProxC[match(samples,metaProxC$Sample_ID),]
#Take a quick look through met to make sure everything is kosher
met 
#Assign groups
group <- met$FOS
Pair <- levels(as.factor(as.character(group)))

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
f <- fdrtool(x=res$PValue,statistic="pvalue",plot=FALSE)
res$f <- f$qval
res.NE_N_C_FvN <- res

####
sum(res$f < 0.05)
sum(res$f < 0.05 & res$logFC > 0)
sum(res$f < 0.05 & res$logFC < 0)
####
## Get overlapping genes
group1 <- res.HC_PvsN_N_N
group2 <- res.HC_PvsN_CvsN_N
group3 <- res.HC_N_CvsN_N

group1 <- group1[group1$f < 0.05,]
group2 <- group2[group2$f < 0.05,]
group3 <- group3[group3$f < 0.05,]

tmp <- c(rownames(group1[group1$logFC < 0,]),
         rownames(group2[group2$logFC < 0,]),
         rownames(group3[group3$logFC > 0,]))
tmp2 <- table(tmp)
sum(tmp2==3)
sum(tmp2==2)
sum(tmp2==1)

genes <- c(genes, names(tmp2[tmp2 == 3]))


res.fosN_fosP <- res
res.fosN_fosP[res.fosN_fosP$logFC > 0 & res.fosN_fosP$f < 0.05,"sig"] <- "P High"
res.fosN_fosP[res.fosN_fosP$logFC < 0 & res.fosN_fosP$f < 0.05,"sig"] <- "N High"
res.fosN_fosP[res.fosN_fosP$f > 0.05 & res.fosN_fosP$PValue < 0.05,"sig"] <- "p < 0.05"
res.fosN_fosP[res.fosN_fosP$f > 0.05 & res.fosN_fosP$PValue > 0.05,"sig"] <- "ns"

write.table(x=rownames(res.fosN_fosP[res.fosN_fosP$f < 0.05 & res.fosN_fosP$logFC > 0,]),file="~/Documents/test.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(x=res.fosN_fosP,file="~/Documents/SalkProjects/ME/WCellActivation_Manuscript/SummaryExcelTables/fosN_fosP_res.txt",quote=FALSE)
####### Volcano
res.fosN_fosP[,"spotME"] <- FALSE
res.fosN_fosP["Egr1","spotME"] <- TRUE
#pdf("~/Documents/SalkProjects/ME/WCellActivation_Manuscript/WCellMan_tiff/fosN_fosP_edgeR_tmp.pdf",width=9,height=7)
ggplot(res.fosN_fosP, aes(logFC, -log(PValue),colour = spotME))+
  geom_point(shape = 1)+
  #geom_point(size=2, shape = 1, colour = "black")+
  xlim(c(-15,15))+
  #scale_colour_manual(values =c("blue","grey","black","red"))+
  labs(title=c("EE Nuclei"))+
  theme_bw(base_size=15)+
  theme(panel.border = element_rect(colour=c("black"),size=2),
        axis.ticks = element_line(size=1.5))+
  geom_vline(xintercept = 0, linetype = "dashed")+
  geom_hline(yintercept = -log(0.05), linetype = "dashed")+
  geom_hline(yintercept = -log(0.004038490), linetype = "dashed")+
  theme(axis.text.x = element_text(angle = 70, hjust = 1))
#dev.off()
