### monocle permutation
###########################
## Monocle
###########################
library(monocle)
library(ggplot2)
library(parallel)
#save(list=c("linear_monocle.dg","linear_monocle.ca1","res.n.ee.left_v_right","res.n.ee.right","res.n.ee.left","my.data5.dg","my.data5.ca1","pheno.dg","pheno.ca1"),file="~/Documents/SalkProjects/ME/ShortLongSingature/SLSig_R/sl_monocle.rda",compress=TRUE)
load("~/Documents/SalkProjects/ME/ShortLongSingature/SLSig_R/sl_monocle.rda")
###########################
## Functions
###########################
pseudoPlot <- function(gene){
  tmp <- pheno
  tmp[,"TPM"] <- as.numeric(dat[gene,])
  
  plt <- ggplot(tmp, aes(Pseudotime, TPM,colour = FOS))+
    geom_point()+
    labs(title = gene)+
    geom_smooth(size = 1, colour = "black",group = "1")+
    theme_bw(base_size = 15)+
    xlab("Pseudotime")+
    theme(panel.border = element_rect(colour=c("black"),size=2),
          axis.ticks = element_line(size=1.5))
  return(plt)
}
modelMe <- function(g){
  model <- summary(lm(value ~ ps, tmp[tmp$X2 == g,]))
  model <- model$coefficients
  return(c(model[1,],model[2,]))
}

###########################
## Step 1) Create CellDataSet object
###########################

###Monocle requires normalized counts
samples <- metaProxC[ metaProxC$Brain_Region == "DG" &  metaProxC$alignable >  500000 & metaProxC$Smartseq2_RT_enzyme_used == "ProtoscriptII" ,"Sample_ID"]
outlier <- rownames(pheno[which(pheno$Pseudotime > 12 & pheno$FOS == "F.EE"),])
#samples <- samples[-c(match(outlier, samples))]
exprs <- dat <- na.exclude(tpmProxC[, samples])
met <- metaProxC[match(samples,metaProxC$Sample_ID),]

lab <- as.character(met$FOS)


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
#a <- RES[[8]]
genes <- c(activitygenes)


marker_genes <- row.names(subset(fData(my.data2),gene_short_name %in% genes))
my.data3 <- setOrderingFilter(my.data2, marker_genes)

###########################
## Step3) Dimernsionality reduction
###########################
my.data4 <- reduceDimension(my.data3, use_irlba=FALSE)
my.data5 <- orderCells(my.data4,num_paths=4)
pheno <- pData(my.data5)
pheno$FOS <- paste(as.character(met$FOS), as.character(met$Mouse_condition),sep = ".")
samples.1 <- as.character(pheno[ ,"cells"])
samples <- colnames(dat)[match(samples.1, colnames(exprs))]
rownames(pheno) <- samples
###########################
## Step4) Gene models
###########################

#full_model_fits <- fitModel(my.data5[marker_genes,], modelFormulaStr = "expression ~VGAM::s(Pseudotime)")
#expression_curve_matrix <- responseMatrix(full_model_fits)
#clusters <- clusterGenes(expression_curve_matrix, k = 5)
#plot_clusters(my.data5,clusters)

metaProxC[colnames(dat),"State_DG"] <- pData(my.data5)$State

###########################
## Step5) Plot results
###########################
pData(my.data5)$FOS  <- paste(as.character(met$FOS), as.character(met$Mouse_condition),sep = ".")
g <- "Pcdha5"
pData(my.data5)$gene  <-as.numeric(dat[g,])
plot_spanning_tree2(my.data5,color_by="gene",tit =g )

plot_spanning_tree(my.data5,color_by = "State" )

###
pseudoPlot("Acvr1c")

###########################
## Extract samples of interest
###########################

samples.1 <- as.character(pheno[pheno$State == 4 & pheno$FOS == "N.EE" | pheno$State == 5  & pheno$FOS == "N.EE"
                   #pheno$State == 6  
                   ,"cells"])
samples <- colnames(dat)[match(samples.1, colnames(exprs))]
group <- pheno[match(samples.1,pheno$cells),"State"]
group <- as.factor(group == 4)
Pair <- levels(as.factor(group))

##
tmp <- dat[,samples]
tmp <- melt(t(tmp[rowSums(tmp) > 0, ]))
tmp$ps <- pheno[samples,"Pseudotime"]

####
tmp$group <- ifelse(group == TRUE, "left","right")
gene <- "Bdnf"
ggplot(tmp[tmp$X2 == gene,], aes(group, value))+
  geom_violin()+
  geom_point()+
  labs(title = gene)

####################
## 

#Run regression
p <- apply(dat[genes,samples],MARGIN = 1,FUN = propExp )
genes <- names(p[p > 0.3])
cl <- makeCluster(getOption("cl.cores", 3))
clusterExport(cl=cl, varlist=c("tmp"))
tmp3 <- as.data.frame(do.call("rbind",parLapply(cl=cl,X=genes, fun=modelMe)))
stopCluster(cl)
#format dataframe
colnames(tmp3) <- as.vector(t(outer(c("int","ps"),c("est","err","t","p"),paste,sep=".")))
tmp3$f <- p.adjust(tmp3$ps.p)
rownames(tmp3) <- genes
tmp3 <- tmp3[order(tmp3$ps.p), ]
tmp3$originalFC <- RES[[8]][rownames(tmp3),"logFC"]
tmp3$originalF<- RES[[8]][rownames(tmp3),"f"]
linear_monocle.ca1_allAct <- tmp3
#save results on main page

#######################
##
res <- cbind(res.n.ee.left, res.n.ee.right[rownames(res.n.ee.left),])
colnames(res) <- c(paste(colnames(res),c(1:ncol(res)),sep = "."))
ggplot(res, aes(logFC.1, logFC.5))+
  geom_point(aes(logFC.1, logFC.5), data = res[res$f.4 > 0.05,], colour = "grey")+
  geom_point(aes(logFC.1, logFC.5), data = res[res$f.4 < 0.05 & res$PValue.7 > 0.05,], colour = "red")+
  geom_point(aes(logFC.1, logFC.5), data = res[res$f.8 < 0.05 & res$PValue.3 > 0.05,], colour = "green")+
  theme_bw(base_size = 15)+
  xlab("Fold Change\nleft FOS- NE vs else")+
  ylab("Fold Change\nright FOS- NE vs else")+
  labs(title =c("Distinguishing Left and Right\nFOS- NE populations"))

a <- res[res$logFC.1 < 0 & res$f.4 < 0.05 & res$PValue.7 > 0.05,]
a <- a[order(a$PValue.3),]


