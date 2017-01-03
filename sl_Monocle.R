### monocle permutation
###########################
## Monocle
###########################
library(monocle)
library(ggplot2)
library(parallel)
#save(list=c("linear_monocle.dg","linear_monocle.ca1","res.n.ee.left_v_right","res.n.ee.right","res.n.ee.left","my.data5.dg","my.data5.ca1","pheno.dg","pheno.ca1"),file="~/Documents/SalkProjects/ME/ShortLongSingature/SLSig_R/sl_monocle.rda",compress=TRUE)
#load("~/Documents/SalkProjects/ME/ShortLongSingature/SLSig_R/sl_monocle.rda")
##### Below is newer data as of 7/20/2016
#save(list=c("allnuclei","inhib","vip","IN"),file="~/Documents/SalkProjects/ME/ShortLongSingature/SLSig_R/sl_monocle2.rda",compress=TRUE)
#load("~/Documents/SalkProjects/ME/ShortLongSingature/SLSig_R/sl_monocle2.rda")
#save(list=c("my.data5_ca3"),file="~/Documents/SalkProjects/ME/ShortLongSingature/SLSig_R/sl_monocle3.rda",compress=TRUE)
#load("~/Documents/SalkProjects/ME/ShortLongSingature/SLSig_R/sl_monocle3.rda")
###########################
## Functions
###########################
pseudoPlot <- function(gene,color_by = "FOS", shape_by = NULL , COLORS = NULL, reverse = FALSE, Max = NULL){
  tmp <- pData(my.data5)
  tmp[,"TPM"] <- as.numeric(dat[gene,])
  tmp[,"color_by"] <- tmp[,color_by]
  if(reverse == TRUE){
    tmp$Pseudotime <- -tmp$Pseudotime
  }
  if(!is.null(shape_by)){
    tmp[,"shape_by"] <- tmp[,shape_by]
    plt <- ggplot(tmp, aes(Pseudotime, TPM,colour = color_by, shape= shape_by))+
      geom_point()
  }else{
  plt <- ggplot(tmp, aes(Pseudotime, TPM,colour = color_by))+
    geom_point()
  }
   plt <- plt + labs(title = gene)+
    geom_smooth(size = 1, colour = "black",group = 1)+
    theme_bw(base_size = 15)+
    xlab("Pseudotime")+
    theme(panel.border = element_rect(colour=c("black"),size=2),
          axis.ticks = element_line(size=1.5))
     if(!is.null(COLORS)){
       plt <- plt + scale_colour_manual(values = COLORS)
     }
   if(!is.null(Max)){
     plt <- plt + ylim(c((min(tmp$TPM) - sd(tmp[tmp$TPM <=0,"Pseudotime"])), Max))
   }

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
samples <- rownames(metaProxC[metaProxC$Brain_Region == "DG" & metaProxC$FOS != "L" & metaProxC$cluster_outlier == "in" & metaProxC$outliers == "in",])
exprs <- dat <- na.exclude(tpmProxC[, samples])
met <- metaProxC[samples,]

colnames(exprs) <- nm2 <- paste(met$Subgroup2, c(1:nrow(met)),sep = "_")
df <- data.frame(cells = nm2,labels = as.character(met$FOS))
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
#expressed_genes <- row.names(subset(fData(my.data2),num_cells_expressed >=2))

marker_genes <- row.names(subset(fData(my.data2),gene_short_name %in% genes))
my.data3 <- setOrderingFilter(my.data2, marker_genes)

###########################
## Step3) Dimernsionality reduction
###########################
my.data4 <- reduceDimension(my.data3, use_irlba=FALSE)
my.data5 <- orderCells(my.data4,num_paths=1)
pheno <- pData(my.data5)
pheno$FOS <- paste(as.character(met$FOS), as.character(met$Mouse_condition),sep = ".")
samples.1 <- as.character(pheno[ ,"cells"])
samples <- colnames(dat)[match(samples.1, colnames(exprs))]
rownames(pheno) <- samples

###########################
## Step5) Plot results
###########################
pData(my.data5)$FOS <- paste(as.character(met$FOS),as.character(met$Mouse_condition),sep = ".")
pData(my.data5)$Subgroup2 <- factor(met$Subgroup2,levels = c("CA1","CA1b","CA3","DG","VIP","IN","CA2"))
pData(my.data5)$samples <- colnames(exprs)
pData(my.data5)$Mouse_condition <- paste(as.character(met$Mouse_condition),sep = ".")
pData(my.data5)$gene <- as.numeric(dat["Arc",])
i <- "N.EE"
pData(my.data5)$group <- ifelse(pData(my.data5)$FOS == i, i,"other")

plot_spanning_tree2(my.data5,
                    color_by = "gene", tit = "Activity in the DG")# , COLORS = c("red","grey"))#,
                    #COLORS = c("red","#f98e04","#e1ba04","blue","skyblue"))#,
                    #SHAPES = c(15,17))#"#f98e04","#e1ba04","blue","skyblue"))#, Subset = "EE", Subset_col = "Mouse_condition")

###
#tiff(filename = "~/Documents/SalkProjects/ME/ShortLongSingature/MolecDissec_Figs_Tables/Figures_vE/pseudotime/pseudo_CA3_Neat1.tiff",width = 10,height = 4,units = 'in',res = 500)
pseudoPlot("Dgat2l6",color_by = "FOS",reverse = TRUE )
#dev.off()
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



############
# 
S_matrix <- reducedDimS(my.data5)
lib_info_with_pseudo <- pData(my.data5)
ica_space_df <- data.frame(t(S_matrix[c(x, y), ]))
colnames(ica_space_df) <- c("ICA_dim_1", "ICA_dim_2")
ica_space_df$sample_name <- row.names(ica_space_df)
df <- as.data.frame(merge(ica_space_df, lib_info_with_pseudo, 
                                 by.x = "sample_name", by.y = "row.names"))
df$labels2 <- df$labels
df[df$labels2 == "L","labels2"] <- "F"
ggplot(df[df$Subgroup2 != "CA2",], aes(ICA_dim_2, fill = Subgroup2, alpha = 0.5))+
  geom_density()
