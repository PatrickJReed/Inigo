### monocle permutation
###########################
## Monocle
###########################
library(monocle)
library(ggplot2)
#save(list=c(),file="",compress=TRUE)

###########################
## Functions
###########################
pseudoPlot <- function(gene){
  tmp <- pheno
  tmp[,"TPM"] <- as.numeric(dat[gene,])
  
  plt <- ggplot(tmp, aes(Pseudotime, TPM,colour = brain_regions))+
    geom_point()+
    labs(title = gene)+
    geom_smooth(size = 1, colour = "black",group = "1")+
    theme_bw(base_size = 15)+
    xlab("Pseudotime")+
    theme(panel.border = element_rect(colour=c("black"),size=2),
          axis.ticks = element_line(size=1.5))
  return(plt)
}
###########################
## Step 1) Create CellDataSet object
###########################

###Monocle requires normalized counts
samples <- metaProxC[metaProxC$FOS == "N" &  metaProxC$Mouse_condition == "HC" & metaProxC$alignable >  500000 & metaProxC$Smartseq2_RT_enzyme_used == "ProtoscriptII" ,"Sample_ID"]
exprs <- tpmProxC[, samples]
met <- metaProxC[match(samples,metaProxC$Sample_ID),]

lab <- as.character(met$Brain_Region)


tmp <- as.numeric(scale(colSums(exprs))[,1])

colnames(exprs) <- nm2 <- paste(met$Brain_Region, c(1:nrow(met)),sep = "_")
df <- data.frame(cells = nm2,labels = lab)
rownames(df) <- nm2
phenoData <- new("AnnotatedDataFrame",data=df)
df <- data.frame(gene_short_name=rownames(exprs))
rownames(df) <- rownames(exprs)
featureData <- new("AnnotatedDataFrame",data=df)

my.data <- newCellDataSet(exprs,
                          phenoData = phenoData,
                          featureData = featureData)


###########################
## Step2) Identify gene set
###########################
my.data2 <- detectGenes(my.data,min_expr=1)
expressed_genes <- row.names(subset(fData(my.data2),num_cells_expressed >=10))
##
genes <- celltypegenes


marker_genes <- row.names(subset(fData(my.data2),gene_short_name %in% genes))
my.data3 <- setOrderingFilter(my.data2, marker_genes)

###########################
## Step3) Dimernsionality reduction
###########################
my.data4 <- reduceDimension(my.data3, use_irlba=FALSE)
my.data5 <- orderCells(my.data4,num_paths=4)
pheno <- pData(my.data5)
pheno$brain_region <- met$Brain_Region

###########################
## Step4) Gene models
###########################

#full_model_fits <- fitModel(my.data5[marker_genes,], modelFormulaStr = "expression ~VGAM::s(Pseudotime)")
#expression_curve_matrix <- responseMatrix(full_model_fits)
#clusters <- clusterGenes(expression_curve_matrix, k = 5)
#plot_clusters(my.data5,clusters)

###########################
## Step5) Plot results
###########################
pData(my.data5)$brain_regions  <- as.character(met$Brain_Region)
pData(my.data5)$brain_regions[met$subgroup == "IN"] <- "IN"
plot_spanning_tree(my.data5,color_by="brain_regions")
#pData(my.data5)$Gabra1  <- as.numeric(tpmQC["Atm",colnames(my.data5)])
#plot_spanning_tree2(my.data5,color_by="Gabra1", cellsize=3,tit="EE DGC Nuclei")
###
pseudoPlot("Synpr")
