###########################
## ShortTerm DataProcessing
###########################
## libraries
library(monocle)
library(ggplot2)

###########################
#### Separation Into States by ICA
###########################
#### Step 1) 
exprs <- tpmProx
lab <- as.character(labelsProx$fos)

df <- data.frame(cells = colnames(exprs),labels = lab,row.names=colnames(exprs))
phenoData <- new("AnnotatedDataFrame",data=df)
df <- data.frame(gene_short_name=rownames(exprs))
rownames(df) <- rownames(exprs)
featureData <- new("AnnotatedDataFrame",data=df)

my.data <- newCellDataSet(exprs,phenoData = phenoData, featureData = featureData)


## Step2) Look at gene categories

my.data2 <- detectGenes(my.data,min_expr=1)
p <- apply(X=tpmProx,MARGIN=1,FUN=propExp)
genes <- c(rownames(tpmProx[p < 0.8 & p > 0.5,]))

marker_genes <- row.names(subset(fData(my.data2),gene_short_name %in% genes))
my.data3 <- setOrderingFilter(my.data2, marker_genes)
my.data4 <- reduceDimension(my.data3, use_irlba=FALSE)
my.data5 <- orderCells(my.data4,num_paths=2)
plot_spanning_tree(my.data5,color_by="State")

#Export state information
ICA.states <- pData(my.data5)

###########################
#### Monocle Differential Expression based on states
###########################

diff_test_res <- differentialGeneTest(my.data5[marker_genes,],fullModelFormulaStr="expression ~VGAM::s(Pseudotime)")
diff_test_res <- na.exclude(diff_test_res[order(diff_test_res$qval),])
diff_test_res <- diff_test_res[order(diff_test_res$pval),]


