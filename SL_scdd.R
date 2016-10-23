#library(devtools)
#devtools::install_github("kdkorthauer/scDD")
library(scDD)
library(Biobase)

#construct starting object
samples <- rownames(t)
cts <- as.matrix(countProxC[,samples])
pheno <- data.frame(condition = as.factor(t$glut_gaba),row.names = samples)

eset.scdd <- ExpressionSet(assayData = cts,  
                           phenoData = as(pheno, "AnnotatedDataFrame"))

rm(cts)

#Create list of prior parameters for model fitting, why are these the defaults? Should play around with this
prior_param=list(alpha=0.01, mu0=0, s0=0.01, a0=0.01, b0=0.01)

#Detect and classify differentially distributed genes
