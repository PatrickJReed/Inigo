##################
# SCDE for Prox1+/-
# Slinker 12/2014
##################
library(ggplot2)
library(reshape)
library(scde)
##################
##############################################
## Step 1: Load in Count Table
##############################################

countTable2 <- countProx[rowSums(countProx) > 0,]
labels <- labelsProx$fos
##############################################
## Step 2: Fit error model
##############################################

# number of local process cores to use in processing
n.cores <- 2;
# calculate models
o.ifm <- scde.error.models(counts=countTable2,
                           groups=labels,
                           n.cores=n.cores,
                           threshold.segmentation=T,
                           save.crossfit.plots=F,
                           save.model.plots=F,
                           verbose=1);
# Filter out cells with really bad fits
valid.cells <- o.ifm$corr.a >0;
table(valid.cells)
# estimate gene expression prior
o.prior <- scde.expression.prior(models=o.ifm,
                                 counts=countTable2,
                                 length.out=400,
                                 show.plot=F)

# get expression magntiude estimates
o.fpm <- as.data.frame(scde.expression.magnitude(o.ifm,counts=countTable2))

