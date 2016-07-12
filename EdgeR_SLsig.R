## EDGER
library(pcaMethods)
library(edgeR)
library(scatterplot3d)
library(fdrtool)
library(ggplot2)
### other function
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
  fit <- glmQLFit(y=cds,design)
  lrt <- glmLRT(fit, coef=Coef)#interaction term = 4
  edg <- data.frame(lrt$table)
  edg <- edg[order(edg$PValue),]
  f <- p.adjust(edg$PValue)
  edg$f <- f
  return(edg)
}
exact <- function(dat, variable, Pair){
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
  return(res)
}

########################
###Load and Format Data
########################
#save(list = c("celltypeorder.dg", "celltypeorder.pin", "celltypeorder.ca1", "celltypeorder.neg","celltypeorder","activitygenes","celltypegenes","celltypegenes.hdg", "celltypegenes.dg", "celltypegenes.ca1", "celltypegenes.neg","celltypegenes.ca23","celltypegenes.in","activitygenes.ca1","activitygenes.dg","activitygenes.hdg","activitygenes.neg","RES","RES2"),file = "~/Documents/SalkProjects/ME/ShortLongSingature/SLSig_R/edgeR_slsig.rda",compress = TRUE)
#load("~/Documents/SalkProjects/ME/ShortLongSingature/SLSig_R/edgeR_slsig.rda")
###
samples <- rownames(metaProxC[ metaProxC$Brain_Region == "DG" & metaProxC$Context2 != "C" & metaProxC$outliers == "in",])
dat <- na.exclude(countProxC[, samples])
dat <- dat[rowSums(dat) > 0,]
met <- metaProxC[match(samples,metaProxC$Sample_ID),]
###################
#Assign groups
###################
group <- met$Mouse_condition == "5hpA"
Pair <- levels(as.factor(as.character(group)))
###################
# Test genes
###################
#!!!Run exact test
res <- exact(dat, group, Pair)
############################
### Keep track of gene lists 
############################

##RES = using the 500000 aligned reads cutoff
#Combine all results into a single list
#RES <- list(res.HC_P_NvC_N,res.HC_PvsN_NvsC_N ,res.HC_PvsN_N_N ,res.HC_PvN_C_N ,res.HC_PvsN_CvsN_N ,res.HC_N_CvsN_N ,res.NE_N_N_FvN, res.NE_N_C_FvN, res.NE_P_C_FvN, res.NE_P_N_FvN)
#names(RES) <- c("res.HC_P_NvC_N","res.HC_PvsN_NvsC_N" ,"res.HC_PvsN_N_N" ,"res.HC_PvN_C_N" ,"res.HC_PvsN_CvsN_N" ,"res.HC_N_CvsN_N" ,"res.NE_N_N_FvN", "res.NE_N_C_FvN", "res.NE_P_C_FvN", "res.NE_P_N_FvN")
#RES[[(length(RES) + 1)]] <- res.HCNE_P_N_FN
#names(RES)[[(length(RES))]] <- "res.HCNE_P_N_FN"

##RES2 = using the 100000 aligned reads cutoff
##RES3 = groups from 'subgroup2'
i <- 3
RES3[[i]] <- res.HC_IN
names(RES3)[i] <- "res.HC_IN"
#### Find significant genes in each group
i <- 1
sum(RES[[i]]$f < 0.05)
sum(RES[[i]]$f < 0.05 & RES[[i]]$logFC > 0)
sum(RES[[i]]$f < 0.05 & RES[[i]]$logFC < 0)

#################
## Get celltype-specific overlapping genes
#################
#pick groups with the pertinent comparisons
# DG = 1(<0),4(>0),5(>0)
# pIN =  1 (>0), 2 (>0), 3 (>0)
# CA1 = 2 (<0), 4 (<0), 6 (<0)
# NEG = 3 (<0), 5 (<0), 6 (>0) 
group1 <- RES2[[3]] 
group2 <- RES2[[5]]
group3 <- RES2[[6]]

nm <- table(c(rownames(group1),
              rownames(group2),rownames(group3)))
nm <- names(nm[nm == 3])
group1 <- group1[nm,]
group2 <- group2[nm,]
group3 <- group3[nm,]

###Gene names
a <-rownames(group2[group1$logFC < 0 & group1$f < 0.05 &
                                     group2$logFC < 0 & group2$f < 0.05 &
                                     group3$logFC > 0 & group3$f < 0.05  
                                    ,])

#celltypegenes.neg  <- a 
#celltypegenes <- c(celltypegenes.pin, celltypegenes.dg, celltypegenes.ca1, celltypegenes.neg)
###Gen Order
all <- rbind(group1[a,], group2[a,], group3[a,])
all$gene <- a
all <- all[order(abs(all$logFC), decreasing=TRUE),]
all$dup <- duplicated(all$gene)
celltypeorder.neg <- abs(all[all$dup == FALSE,"logFC"])


celltypeorder <- c(celltypeorder.pin, celltypeorder.dg, celltypeorder.ca1,celltypeorder.neg)
celltypeorder2 <- c(20 * celltypeorder.pin, celltypeorder.dg, 1000*celltypeorder.ca1,80*celltypeorder.neg)


#### For activity genes
hdg <- RES[[10]] 
dg <- RES[[9]]
ca1 <- RES[[8]]
neg <- RES[[7]]

nm <- table(c(rownames(hdg),
  rownames(dg),rownames(ca1),rownames(neg)))
nm <- names(nm[nm == 4])
hdg <- hdg[nm,]
dg <- dg[nm,]
ca1 <- ca1[nm,]
neg <- neg[nm,]


foshigenes.hdg <-rownames(hdg[hdg$logFC < 0 & hdg$f < 0.05 & dg$PValue > 0.05 & ca1$PValue > 0.05 & neg$PValue > 0.05,])
foshigenes.dg <- rownames(dg[dg$logFC < 0 & dg$f < 0.05 & hdg$PValue > 0.05 & ca1$PValue > 0.05 & neg$PValue > 0.05,])
foshigenes.ca1 <- rownames(ca1[ca1$logFC < 0 & ca1$f < 0.05 & hdg$PValue > 0.05 & dg$PValue > 0.05 & neg$PValue > 0.05,])
foshigenes.neg <- rownames(neg[neg$logFC < 0 & neg$f < 0.05 & hdg$PValue > 0.05 & ca1$PValue > 0.05 & dg$PValue > 0.05,])


foshigenes <- na.exclude(unique(c(foshigenes.dg,foshigenes.hdg,foshigenes.ca1,foshigenes.neg)))
####NEg genes
samples <- metaProxC[metaProxC$PROX1 == "N"  & metaProxC$CTIP2 == "N" &  metaProxC$Mouse_condition == "HC" & metaProxC$alignable >  500000 & metaProxC$Smartseq2_RT_enzyme_used == "ProtoscriptII" ,"Sample_ID"]
p.neg <- apply(X = dat[,samples],MARGIN = 1, FUN = propExp)
cvar <- apply(X = dat[,samples],MARGIN = 1, FUN = cvarNoZero)
m.neg <- apply(X = dat[,samples],MARGIN = 1, FUN = meanNoZero)
genes <- names(p.neg[cvar > 0.5 & cvar < 0.6 & m.neg > 2])

plot(apply(X = dat[,samples],MARGIN = 1, FUN = cvarNoZero),p.neg)


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
