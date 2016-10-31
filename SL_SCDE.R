##################
# SCDE for Prox1+/-
# Slinker 12/2014
##################
library(ggplot2)
library(reshape)
library(scde)
##################
#save(list = c("ca1.ediff2","dg.ediff2","ca3.ediff2","in.ediff2","vip.ediff2","ca1_NE.ediff2","dg_NE.ediff2","in_NE.ediff2","vip_NE.ediff2"),file = "~/Documents/SalkProjects/ME/ShortLongSingature/SLSig_R/scde_difexp.rda",compress = TRUE)
##############################################
## Step 1: Load in Count Table
##############################################
samples <- rownames(t[t$T1 < 0 & t$T2 > -10,])
t.1 <- t[t$T1 < 0 & t$T2 > -10,]
#samples <- rownames(metaProxC[metaProxC$FOS != "L"  & metaProxC$Arc_2.5 != "greater" & metaProxC$Context1 == "none" & metaProxC$cluster_outlier == "in" & metaProxC$outliers == "in",])
dat <- na.exclude(countProxC[, samples])
dat <- dat[rowSums(dat) > 0,]
met <- metaProxC[samples,]
countTable2 <- t(apply(round(dat[rowSums(dat) > 0,],0),1, as.integer))
colnames(countTable2) <- colnames(dat)

###################
#Assign groups
###################
labels <- ifelse(t[samples,"a"] == 2 ,yes = "F",no = "N")

##############################################
## Step 2: Fit error model
##############################################

# number of local process cores to use in processing
# For goodness sake don't go above 2 or your computer will die...blegh
n.cores <- 2
# calculate models
o.ifm <- scde.error.models(counts=countTable2,
                           groups=labels,
                           n.cores=n.cores,
                           threshold.segmentation=T,
                           save.crossfit.plots=F,
                           save.model.plots=F,
                           verbose=1)
# Filter out cells with really bad fits
valid.cells <- o.ifm$corr.a >0
table(valid.cells)
# estimate gene expression prior
o.prior <- scde.expression.prior(models=o.ifm,
                                 counts=countTable2,
                                 length.out=400,
                                 show.plot=F)

##############################################
## Step 3: Test Differential Expression
##############################################
# define two groups of cells
groups <- as.factor(labels)
names(groups) <- row.names(o.ifm);
# run differential expression tests on all genes.
ediff <- scde.expression.difference(o.ifm,
                                    countTable2,
                                    o.prior,
                                    groups=groups,
                                    n.randomizations=100,
                                    n.cores=n.cores,
                                    verbose=1)
ediff2 <- ediff[order(ediff$Z,decreasing=T),]
ediff2$rank_high <- c(1:nrow(ediff2))
ediff2 <- ediff[order(ediff$Z,decreasing=F),]
ediff2$rank_low <- c(1:nrow(ediff2))

tail(ediff2)
# Plot posterior for a single gene

gene <- "Arc"
scde.test.gene.expression.difference(gene,
                                       models=o.ifm,
                                       counts=countTable2,
                                     n.randomizations=100
                                       ,prior=o.prior,
                                     groups=groups)

#dev.off()

