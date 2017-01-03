library(monocle)
library(ggplot2)
library(parallel)
pseudoPlot <- function(gene,color_by = "FOS", shape_by = NULL , COLORS = NULL){
  tmp <- pData(my.data5)
  tmp[,"TPM"] <- as.numeric(exprs[gene,])
  tmp[,"color_by"] <- tmp[,color_by]
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
  
  return(plt)
}
modelMe <- function(g){
  model <- summary(lm(value ~ ps, tmp[tmp$X2 == g,]))
  model <- model$coefficients
  return(c(model[1,],model[2,]))
}
#save(list = c("my.data5_subset","my.data5_all"),file =  "~/Documents/SalkProjects/ME/ShortLongSingature/priming_monocle.rda",compress = TRUE)
#Obtain list of genes that are associated with activity, but not obviously increasing from basal state
#   We want to try to bias for more genes that have the potential to be priming
a <- RES[["DG"]]
genes <- table(c(rownames(a[a$f < 0.05,]), rep(SigGenes[[2]],3)))
genes <- names(genes[genes ==1])
#outlier
#outlier <- "DG_48"

###Monocle requires normalized counts
samples <- metaProxC[ metaProxC$FOS != "L" & metaProxC$Context1 == "none" & metaProxC$Subgroup2 == "DG"  & metaProxC$outliers == "in","Sample_ID"]
exprs.1 <- dat <- na.exclude(tpmProxC[, samples])
met <- metaProxC[match(samples,metaProxC$Sample_ID),]

colnames(exprs.1) <-  paste(met$Subgroup2, c(1:nrow(met)),sep = "_")
Out <- table(c(outliers, colnames(exprs.1), colnames(exprs.1)))
In <- names(Out[Out==2])
exprs <- exprs.1[,In]
met <- metaProxC[match(samples,metaProxC$Sample_ID),]
rownames(met) <-  paste(met$Subgroup2, c(1:nrow(met)),sep = "_")

met <- met[In,]
nm2 <- colnames(exprs)
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
expressed_genes <- row.names(subset(fData(my.data2),num_cells_expressed >=10))
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
## Step5) Plot results
###########################
pData(my.data5)$FOS <- paste(as.character(met$FOS),as.character(met$Mouse_condition),sep = ".")
pData(my.data5)$Subgroup2 <- factor(met$Subgroup2,levels = c("CA1","CA1b","CA3","DG","VIP","IN","CA2"))
pData(my.data5)$samples <- colnames(exprs)
pData(my.data5)$Mouse_condition <- paste(as.character(met$Mouse_condition),sep = ".")
pData(my.data5)$pickMe <- met$Mouse_condition == "EE" & met$Subgroup2 == "CA1" & met$FOS == "N"
plot_spanning_tree2(my.data5,
                    color_by = "FOS",tit = "DG" ,
                    COLORS = c("red","blue","skyblue"))
                    #COLORS = c("red","#f98e04","#e1ba04","blue","skyblue"))
###
pseudoPlot("Prdm5",color_by = ,COLORS = c("red","blue","skyblue"))

 plot_spanning_tree2(my.data5)
### Get info
S_matrix <- reducedDimS(my.data5)
lib_info_with_pseudo <- pData(my.data5)
ica_space_df <- data.frame(t(S_matrix))
colnames(ica_space_df) <- c("ICA_dim_1", "ICA_dim_2")
ica_space_df$sample_name <- row.names(ica_space_df)
df <- as.data.frame(merge(ica_space_df, lib_info_with_pseudo, 
                          by.x = "sample_name", by.y = "row.names"))
df[df$FOS == "F.EE","group"] <- "F"
df[df$FOS == "N.EE","group"] <- "N"
df[df$State == 3 & df$FOS != "F.EE", "group"] <- "close"
df[is.na(df$group), "group"] <- "far"
#rownames(df) <- met[df$samples, "Sample_ID"]

#FOS+/FOS- 
a1 <- RES[["DG"]]
#FOS+/allhc
samples <- c(rownames(df[df$FOS == "N.HC",]), rownames(df[df$FOS == "F.EE",]))
dat <- na.exclude(countProxC[, samples])
dat <- dat[rowSums(dat) > 0,]
met <- metaProxC[match(samples,metaProxC$Sample_ID),]
group <- met$FOS == "F"  
Pair <- levels(as.factor(as.character(group)))
a2 <- exact(dat, group, Pair)
#FOS+/close
samples <- c(rownames(df[df$group == "close",]), rownames(df[df$FOS == "F.EE",]))
dat <- na.exclude(countProxC[, samples])
dat <- dat[rowSums(dat) > 0,]
met <- metaProxC[match(samples,metaProxC$Sample_ID),]
group <- met$FOS == "F"  
Pair <- levels(as.factor(as.character(group)))
a3 <- exact(dat, group, Pair)
#FOS+/far
samples <- c(rownames(df[df$group == "far",]), rownames(df[df$FOS == "F.EE",]))
dat <- na.exclude(countProxC[, samples])
dat <- dat[rowSums(dat) > 0,]
met <- metaProxC[match(samples,metaProxC$Sample_ID),]
group <- met$FOS == "F"  
Pair <- levels(as.factor(as.character(group)))
a4 <- exact(dat, group, Pair)
#close/far
samples <- c(rownames(df[df$group == "far",]), rownames(df[df$group == "close",]))
dat <- na.exclude(countProxC[, samples])
dat <- dat[rowSums(dat) > 0,]
group <- df[samples, "group"] == "close"  
Pair <- levels(as.factor(as.character(group)))
a5 <- exact(dat, group, Pair)

##### Assign values
i <- table(c(rownames(a1[a1$logFC < 0 & a1$f < 0.05,]),
             rownames(a2[a2$logFC < 0 & a2$f < 0.05,]),
             rownames(a3[a3$logFC < 0 & a3$f < 0.05,]),
             rownames(a4[a4$logFC < 0 & a4$f < 0.05,])
             )
)
i <- names(i[i == 4])

ii <- table(c(rownames(a1[a1$logFC < 0 & a1$f < 0.05,]),
              rownames(a3[a3$f > 0.05,]),
              rownames(a4[a4$logFC < 0 & a4$f < 0.05,])#,
              #rownames(a5[a5$logFC < 0 & a5$f < 0.05,])
)
)
ii <- names(ii[ii == 3])

iii <- table(c(rownames(a1[a1$logFC < 0 & a1$f < 0.05,]),
              rownames(a3[a3$logFC < 0 & a3$f < 0.05,]),
              rownames(a4[a4$f > 0.05,])#,
              #rownames(a5[a5$logFC < 0 & a5$f < 0.05,])
)
)
iii <- names(iii[iii == 3])
