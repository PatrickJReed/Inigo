###################
##Figures for Single Nuclei (Biological Variables)
##Created, Slinker 09/22/2015
###################
library(ggplot2)
library(reshape)
library(plyr)
library(pcaMethods)
library(Rtsne)
###################
#EdgeR Results
#load("~/Documents/SalkProjects/ME/ShortLongSingature/SLSig_resTables/edgeR.res")
#Other raw data is stored in LoadData_ShortTerm.R
#T-SNE results
#save(list = c("t.all","t.hc"),file = "~/Documents/SalkProjects/ME/ShortLongSingature/SLSig_R/tsne.rda",compress = TRUE)
#load(c("~/Documents/SalkProjects/ME/ShortLongSingature/SLSig_R/tsne.rda"))
###############################################
## FUNCTIONS THAT WILL PLOT FOR YOU
###############################################


# PCA 2D ---------------------------------------------------------------------
PC2D <- function(scores,Var, dat, met, colorby = NULL,shapeby = NULL,Colors = NULL, gene = NA){

  if (!is.na(gene)){
    tmp <- data.frame(PC1 = scores[,1], PC2 = scores[,2],
                      gene = as.numeric(dat[gene,]))
    tmp <- cbind(tmp,met)
    if(is.null(shapeby)){
        tmp$shapeby <- 1
      }else{tmp$shapeby <- tmp[,shapeby]
      }
    #pdf("~/Documents/SalkProjects/ME/ShortLongSingature/SLSig_tiff/PCA_HC_NP.pdf",width=8.5,height=7)
    plt <- ggplot(tmp, aes(PC1, PC2,colour = gene,shape = shapeby))+
      geom_point(size=5, alpha=0.7)+
      theme_bw()+
      scale_colour_gradient(high="red",low="blue")+
      labs(title="PCA")+
      xlab(paste("PC1: ", 100 * round(Var[1],2),"%",sep=""))+
      ylab(paste("PC2: ", 100 * round(Var[2],2),"%",sep=""))+
      theme(text=element_text(size=20))+
      theme(panel.border = element_rect(colour=c("black"),size=2),
            axis.ticks = element_line(size=1.5))

  }else if(!is.null(colorby)){
    tmp <- data.frame(PC1 = p@scores[,1], PC2 = p@scores[,2],
                      dusp1 = as.numeric(dat["Meg3",])
                      )
    tmp <- cbind(tmp, met)
    tmp$PC <- paste(tmp$PROX1, tmp$CTIP2, sep = ".")
    tmp$PCF <- paste(tmp$PROX1, tmp$CTIP2, tmp$FOS, sep = ".")
    tmp$colorby <- tmp[,colorby]
    
    if(is.null(shapeby)){
      plt <- ggplot(tmp, aes(PC1, PC2,colour = colorby))+
        geom_point(size=5, alpha=0.7)
    }else{
      tmp$shapeby <- tmp[,shapeby]
      plt <- ggplot(tmp, aes(PC1, PC2,colour = colorby,shape = shapeby))+
      geom_point(size=5, alpha=0.7)
    }
    

    
    #pdf("~/Documents/SalkProjects/ME/ShortLongSingature/SLSig_tiff/PCA_HC_NP.pdf",width=8.5,height=7)
    plt <- plt+
      theme_bw()+
      #scale_colour_gradient(high="red",low="blue")+
      labs(title="PCA")+
      xlab(paste("PC1: ", 100 * round(Var[1],2),"%",sep=""))+
      ylab(paste("PC2: ", 100 * round(Var[2],2),"%",sep=""))+
      theme(text=element_text(size=20))+
      theme(panel.border = element_rect(colour=c("black"),size=2),
            axis.ticks = element_line(size=1.5))
  }else{
    tmp <- data.frame(PC1 = scores[,1], PC2 = scores[,2]
                      )
    
    plt <- ggplot(tmp, aes(PC1, PC2))+
      geom_point(size=5, alpha=0.7)+
      geom_point(size=5, shape=1)+
      theme_bw()+
      labs(title="PCA")+
      xlab(paste("PC1: ", 100 * round(p@R2[1],2),"%",sep=""))+
      ylab(paste("PC2: ", 100 * round(p@R2[2],2),"%",sep=""))+
      theme(text=element_text(size=20))+
      theme(panel.border = element_rect(colour=c("black"),size=2),
            axis.ticks = element_line(size=1.5))
  }
  if (!is.null(Colors)){
    plt <- plt + scale_color_manual(values = Colors)
  }
  return(plt)
}

# SPCA 2D ---------------------------------------------------------------------
SPC2D <- function(dat, met, gene = NA){
  require(nsprcomp)
  p <- nsprcomp(t(dat))
  Loading <- p$rotation
  if (!is.na(gene)){
    tmp <- data.frame(PC1 = p$x[,1], PC2 = p$x[,2],
                      CellType = met[,"prox"],
                      fos = met[,"fos"],
                      gene = as.numeric(dat[gene,]))
    tmp$CellType <- as.character(tmp$CellType)
    tmp[tmp$CellType == "N", "CellType"] <- "CA1"
    tmp[tmp$CellType == "P", "CellType"] <- "DGC"
    
    #pdf("~/Documents/SalkProjects/ME/ShortLongSingature/SLSig_tiff/PCA_HC_NP.pdf",width=8.5,height=7)
    plt <- ggplot(tmp, aes(PC1, PC2,colour = gene,shape = CellType))+
      geom_point(size=5, alpha=0.7)+
      theme_bw()+
      scale_colour_gradient(high="red",low="blue")+
      labs(title="PCA")+
      xlab(paste("PC1: ", 100 * round(p@R2[1],2),"%",sep=""))+
      ylab(paste("PC2: ", 100 * round(p@R2[2],2),"%",sep=""))+
      theme(text=element_text(size=20))+
      theme(panel.border = element_rect(colour=c("black"),size=2),
            axis.ticks = element_line(size=1.5))
  }else{
    tmp <- data.frame(PC1 = p$x[,1], PC2 = p$x[,2],
                      CellType = met[,"prox"],
                      fos = met[,"fos"])
    tmp$CellType <- as.character(tmp$CellType)
    tmp[tmp$CellType == "N", "CellType"] <- "CA1"
    tmp[tmp$CellType == "P", "CellType"] <- "DGC"
    
    #pdf("~/Documents/SalkProjects/ME/ShortLongSingature/SLSig_tiff/PCA_HC_NP.pdf",width=8.5,height=7)
    plt <- ggplot(tmp, aes(PC1, PC2,colour = CellType))+
      geom_point(size=5, alpha=0.7)+
      geom_point(size=5, shape=1)+
      theme_bw()+
      scale_colour_manual(values = c("darkorange2","darkorchid4"))+
      labs(title="PCA")+
      xlab(paste("PC1: ", 100 * round(p@R2[1],2),"%",sep=""))+
      ylab(paste("PC2: ", 100 * round(p@R2[2],2),"%",sep=""))+
      theme(text=element_text(size=20))+
      theme(panel.border = element_rect(colour=c("black"),size=2),
            axis.ticks = element_line(size=1.5))
  }
  return(list(plt,Loading))
  #dev.off()
}
# PCA 1D ---------------------------------------------------------------------
PC1D <- function(scores,dat, met, group,component){
tmp <- data.frame(PC1 = scores[,1], PC2 = scores[,2])
tmp <- cbind(tmp,met)
tmp$group <- tmp[,group]
tmp$component <- tmp[,component]
#pdf("~/Documents/SalkProjects/ME/ShortLongSingature/SLSig_tiff/PCA1_HC_NP.pdf",width=8.5,height=7)
p <- ggplot(tmp, aes(component, fill = group))+
  geom_density(alpha=0.7)+
  theme_bw()+
  labs(title=paste("PC",component,sep = ""))+
  xlab("PC1 Score ")+
  theme(text=element_text(size=20))+
  theme(panel.border = element_rect(colour=c("black"),size=2),
        axis.ticks = element_line(size=1.5))
#dev.off()
return(p)
}
# PC Gene -----------------------------------------------------------------
PCGene <- function(dat, met,gene, component){
  p <- pca(t(dat))
  tmp <- data.frame(PC = p@scores[,component],
                    CellType = met[,"prox"],
                    fos = met[,"fos"],
                    gene = as.numeric(dat[gene,]))
  tmp$CellType <- as.character(tmp$CellType)
  tmp[tmp$CellType == "N", "CellType"] <- "CA1"
  tmp[tmp$CellType == "P", "CellType"] <- "DGC"
  l <- summary(lm(PC ~ gene, tmp))
  l <- signif(l$coefficients[2,4],2)
  #pdf("~/Documents/SalkProjects/ME/ShortLongSingature/SLSig_tiff/PCA1_HC_NP.pdf",width=8.5,height=7)
  p <- ggplot(tmp, aes(gene,PC, colour = fos))+
    geom_point(alpha=0.8)+
    geom_smooth(method="lm",fill = NA)+
    theme_bw()+
    labs(title=paste(paste("PC",component," ~ ",gene,sep=""),paste("p<", l,sep=""),collapse='\n'))+
    xlab("TPM")+
    ylab("PC Score")+
    theme(text=element_text(size=20))+
    theme(panel.border = element_rect(colour=c("black"),size=2),
          axis.ticks = element_line(size=1.5))
  return(p)
}
# tsne --------------------------------------------------------------------


getErrors <- function(x){
  return(strsplit(stdout[x]," ")[[1]][5])
}





i <- 17
TSNE <- Rtsne(as.matrix(t(na.exclude(dat))),initial_dims=5,perplexity=i,theta=0,check_duplicates=FALSE,dims = 3)
t <- as.data.frame(TSNE$Y)
colnames(t) <- c("T1","T2","T3")
t <- cbind(t,met)
t$Brain_Region <- as.character(t$Brain_Region)
t[t$Brain_Region == "CA3_other_negs","Brain_Region"] <- "Neg"
t[t$Brain_Region == "P+Neg","Brain_Region"] <- "pIN (P+C-)"
t$gene <- as.numeric(tpmProxC["Gad1",samples])

#tiff("~/Documents/SalkProjects/ME/ShortLongSingature/SLSig_tiff/tsne_hc.tiff",width = 4,height = 4,units = 'in',res = 300,compression = 'lzw')
Plot3D.TSNE(t)
#dev.off()

#tiff("~/Documents/SalkProjects/ME/ShortLongSingature/SLSig_tiff/tsne_all.tiff",width = 10,height = 8,units = 'in',res = 300,compression = 'lzw')
ggplot(t, aes(T1,T2, colour = Brain_Region))+
  geom_point(alpha = 0.7, size = 3)+
  geom_point(shape = 1, size = 3)+
  theme_bw()+
  xlab("TSNE1")+
  ylab("TSNE2")+
  labs(title="Homecage FOS- Neurons\nt-SNE")+
  theme(text=element_text(size=20))+
  theme(panel.border = element_rect(colour=c("black"),size=2),
      axis.ticks = element_line(size=1.5),
      panel.grid.major = element_line(size = 1))+
  #scale_colour_manual(values= c("#6ca425"))
  #scale_colour_manual(values = c("black","#00c7e4","#a800b3","#6ca425","#e19041"))
  scale_colour_manual(values = c("#00c7e4","#a800b3","#6ca425","#e19041"))
#dev.off()



t[t$X < 0 & t$Y < -5, "group"] <- "IN"
t[is.na(t$group),"group"] <- "EX"

# Inidividual Genes (tpm) -------------------------------------------------
Indiv <- function(gene,dat,met){
  
  tmp <- data.frame(exp = as.numeric(dat[gene,]),
                    Arc = as.numeric(dat["Arc",]))
  tmp <- cbind(tmp, met)
  tmp$FOS <- as.character(tmp$FOS)
  tmp[tmp$fos == "F","FOS"] <- "High"
  tmp[tmp$fos == "L","FOS"] <- "Low"
  tmp[tmp$fos == "N","FOS"] <- "None"
  #pdf("~/Documents/SalkProjects/BenLacar/ManuscriptFigures/Camk4.pdf",width=6,height=5)
  p <- ggplot(tmp, aes(FOS,exp))+
    geom_violin()+#outlier.shape=NA)+
    geom_point(position=position_jitter(width=0.01,height=0))+
    theme_bw(base_size=20)+
    ylab("TPM")+
    #scale_colour_gradient(high="red", low="blue")+
    ylim(c(0,max(tmp$exp) + 0.2*(max(tmp$exp))))+
    theme(panel.border = element_rect(colour=c("black"),size=2),
          axis.ticks = element_line(size=1.5))+
    labs(title=paste(gene,"\n"))+
    facet_grid( Mouse_condition  ~ Brain_Region) 
return(p)
}
IndivSubgroup <- function(gene,dat,met){
  
  tmp <- data.frame(exp = as.numeric(dat[gene,]),
                    Arc = as.numeric(dat["Arc",]))
  tmp <- cbind(tmp, met)
  tmp$FOS <- as.character(tmp$FOS)
  tmp[tmp$fos == "F","FOS"] <- "High"
  tmp[tmp$fos == "L","FOS"] <- "Low"
  tmp[tmp$fos == "N","FOS"] <- "None"
  #pdf("~/Documents/SalkProjects/BenLacar/ManuscriptFigures/Camk4.pdf",width=6,height=5)
  p <- ggplot(tmp, aes(subgroup,exp, colour = FOS))+
    geom_violin()+#outlier.shape=NA)+
    geom_point(position=position_jitter(width=0.01,height=0))+
    theme_bw(base_size=20)+
    ylab("TPM")+
    #scale_colour_gradient(high="red", low="blue")+
    ylim(c(0,max(tmp$exp) + 0.2*(max(tmp$exp))))+
    theme(panel.border = element_rect(colour=c("black"),size=2),
          axis.ticks = element_line(size=1.5))+
    labs(title=paste(gene,"\n"))+
    facet_grid( Mouse_condition  ~ Brain_Region) 
  return(p)
}
IndivProx1Grouped <- function(gene){
  
  tmp <- data.frame(exp = as.numeric(dat[gene, met$prox == "P"]),
                    Arc = as.numeric(dat["Arc", met$prox == "P"]))
  tmp <- cbind(tmp, met)
  tmp$prox <- as.character(tmp$prox)
  tmp[tmp$prox == "P","prox"] <- "DGC"
  tmp[tmp$prox == "N","prox"] <- "CA1"
  tmp$fos <- as.character(tmp$fos)
  tmp[tmp$fos == "F","fos"] <- "High"
  tmp[tmp$fos == "L","fos"] <- "Low"
  tmp[tmp$fos == "N","fos"] <- "None"
  t <- t.test(tmp[tmp$tempGroupingProx1HC == TRUE, "exp"],tmp[tmp$tempGroupingProx1HC == FALSE, "exp"])
  #pdf("~/Documents/SalkProjects/BenLacar/ManuscriptFigures/Camk4.pdf",width=6,height=5)
  p <- ggplot(tmp, aes(tempGroupingProx1HC,exp))+
    geom_boxplot(outlier.shape=NA)+
    geom_point(position=position_jitter(width=0.01,height=0))+
    theme_bw(base_size=20)+
    ylab("TPM")+
    #scale_colour_gradient(high="red", low="blue")+
    ylim(c(0,max(tmp$exp) + 0.2*(max(tmp$exp))))+
    theme(panel.border = element_rect(colour=c("black"),size=2),
          axis.ticks = element_line(size=1.5))+
    labs(title=paste(gene,":",round(t$p.value,2),"\n"))+
    facet_grid( cond  ~ prox) 
  return(p)
}
IndivByDate <- function(gene,dat,met){
  
  tmp <- data.frame(exp = as.numeric(dat[gene,]),
                    Arc = as.numeric(dat["Arc",])
                    )
  tmp <- cbind(tmp, met)
  tmp$prox <- as.character(tmp$prox)
  tmp[tmp$prox == "P","prox"] <- "DGC"
  tmp[tmp$prox == "N","prox"] <- "CA1"
  tmp$fos <- as.character(tmp$fos)
  tmp[tmp$fos == "F","fos"] <- "High"
  tmp[tmp$fos == "L","fos"] <- "Low"
  tmp[tmp$fos == "N","fos"] <- "None"
  tryCatch(
    l <- summary(lm(exp ~ date*fos, tmp))
  )
  #pdf("~/Documents/SalkProjects/BenLacar/ManuscriptFigures/Camk4.pdf",width=6,height=5)
  p <- ggplot(tmp, aes(date,exp,fill = fos))+
    geom_boxplot(outlier.shape=NA)+
    geom_point(position=position_jitter(width=0.01,height=0))+
    theme_bw(base_size=14)+
    ylab("TPM")+
    #scale_colour_gradient(high="red", low="blue")+
    ylim(c(0,max(tmp$exp) + 0.2*(max(tmp$exp))))+
    theme(panel.border = element_rect(colour=c("black"),size=2),
          axis.ticks = element_line(size=1.5))+
    labs(title=paste(gene,"\n","p, date < ",signif(l$coefficients[2,4],2),"\n","p, fos < ",signif(l$coefficients[3,4],2),"\n","p, date:fos < ",signif(l$coefficients[4,4],2)))+
    facet_grid( cond  ~ prox) 
  return(p)
}
Plot3D.TSNE <- function(t,group = "Brain_Region",COLORS = c("#00c7e4","#a800b3","#e19041","#6ca425")){
  require(scatterplot3d)
  #source('http://www.sthda.com/sthda/RDoc/functions/addgrids3d.r')
  t$COL <- t[,group]
  j <- 0
  for (i in levels(as.factor(t[,group]))){
    j <- j + 1
    t[t$COL == i, "COL"] <- COLORS[j]
  }
  t2 <- t[t$Brain_Region == "HDG",]
  scatterplot3d(x = t$T1, y = t$T2,z =  t$T3, pch = 16, 
            color = t$COL, angle = 35,#box= FALSE,
            xlab = "TSNE1",
            ylab = "TSNE2", 
            zlab = "TSNE3")
}
# facet_grid( cond  ~ prox)
#dev.off()
# Correlations (2genes) ---------------------------------------------------
Plot2Genes <- function(a,b,dat,met,group="FOS"){
  tmp <- data.frame(A = as.numeric(dat[a,]),
                    B = as.numeric(dat[b,]),
                    group = met[,group],
                    condition = met$Mouse_condition,
                    dusp1 = as.numeric(dat["Npas4",])
                    )
  
  #pdf("~/Documents/SalkProjects/BenLacar/ManuscriptFigures/Camk4.pdf",width=6,height=5)
  p <- ggplot(tmp, aes(A,B,shape = condition,colour = group))+
    geom_point(size=4)+
    theme_bw(base_size=20)+
    xlab(a)+
    ylab(b)+
    theme(panel.border = element_rect(colour=c("black"),size=2),
          axis.ticks = element_line(size=1.5))
  return(p)
  #dev.off()
}
# heatmap (needs work) ----------------------------------------------------
scaleMe <- function(x){
  scale(x)[,1]
}

heatMe <- function(dat,met,genes,k1= 10 , k2 = NULL, sampleorder = NULL,geneorder = NULL,  cutoff = NULL){
  tmp <- na.exclude(dat[genes,])
  colnames(tmp) <- paste(met$PROX1,met$CTIP2,met$FOS,met$Mouse_condition,c(1:nrow(met)),sep = ".")
  tmp.1 <- t(apply(X=tmp,MARGIN=1,FUN=scaleMe))
  tmp2 <- melt(t(na.exclude(tmp.1)))
    o1 <- sampleorder#cutree(hclust(dist(t(tmp))),k=k1)

  if(!is.null(geneorder)){
    o2 <- geneorder
  }
  if(!is.null(k2)){
    o2 <- cutree(hclust(dist(tmp)),k=k2)
  }else{
    o2 <- 1:length(genes)
  }
  tmp2$o1 <- rep(as.numeric(o1), nrow(tmp))
  tmp2$o2 <- rep(as.numeric(o2), each = ncol(tmp))
  if(!is.null(cutoff)){
    tmp2[tmp2$value > cutoff, "value"] <- cutoff + 0.1
    tmp2[tmp2$value < (-1 * cutoff), "value"] <- (-1 * cutoff) - 0.1
  }
  p1 <- ggplot(tmp2, aes(reorder(X1,o1), reorder(X2,o2) , fill = value))+
    geom_tile()+
    theme(axis.title.x=element_blank(),
          axis.title.y=element_blank()
    )+
    scale_fill_gradient2(high="red",mid = "white",low="blue")+
    #theme_bw(base_size=22)+
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  return(p1)
}
# Volcano Plot ------------------------------------------------------------

## Labels genes that are NS, nominal, or adj significant
Volcano <- function(difexp){
  direction <- res$f < 0.05
  direction[direction == FALSE] <- "NS"
  direction[direction == "NS" & res$PValue < 0.05] <- "p < 0.05"
  direction[direction == TRUE & res$logFC > 0] <- "CA1 high"
  direction[direction == TRUE & res$logFC < 0] <- "DGC high"
  direction <- as.factor(direction)
  res$direction <- direction
  
  levels(res$direction) <- c("NS","p < 0.05","CA1 high","DGC high")
  
  
  #pdf("~/Documents/SalkProjects/ME/ShortLongSingature/SLSig_tiff/difexp_HC_NP.pdf",width=7.5,height=7.5)
  p <- ggplot(res, aes(logFC, -log(PValue), colour = direction))+
    geom_point(alpha = 0.2)+
    geom_point(shape = 1)+
    geom_hline(aes(yintercept = -log(0.05)), linetype="dashed")+
    geom_hline(aes(yintercept = -log(res[res$f > 0.05,"PValue"][1])), linetype="dashed")+
    geom_vline(aes(xintercept = 0), linetype="dashed")+
    theme_bw()+
    scale_colour_manual(values = c("darkorchid4","darkorange2","grey","black"))+
    theme(text = element_text(size=18),panel.border = element_rect(colour=c("black"),size=2),
          axis.ticks = element_line(size=1.5),legend.position = "none")+
    annotate(geom="text",x=c(-5,5),y=75,label=c(paste("n= ",c(sum(direction == "CA1 high"),sum(direction == "DGC high")),sep="")))+
    labs(title="Differential Expression\nHC CA1 vs. DGC")
  #dev.off()
  return(p)
}
# Gene set by Go term


###############################################
### PLOT THESE GUYS
###############################################
#PCA 2D
samples <- metaProxC[metaProxC$Brain_Region == "CA3_other_negs" & metaProxC$alignable >  100000 & metaProxC$Smartseq2_RT_enzyme_used == "ProtoscriptII" ,"Sample_ID"]
                                             #metaProxC$Brain_Region == "CA3_other_negs" & metaProxC$Mouse_condition == "HC" & metaProxC$alignable >  100000 & metaProxC$Smartseq2_RT_enzyme_used == "ProtoscriptII","Sample_ID"]
dat <- na.exclude(tpmProxC[, samples])
met <- metaProxC[match(samples,metaProxC$Sample_ID),]
#gene <- "Meg3"
#Calculate the components
p <- pca(t(dat[,samples]),nPcs = 10)
scores <- as.data.frame(p@scores)
loading <- as.data.frame(p@loadings)
Var <- p@R2
#tiff(filename = "~/Documents/SalkProjects/ME/ShortLongSingature/SLSig_tiff/PCA_HC_N.tiff",width = 6.5,height = 5,units = 'in',res = 300)
PC2D(scores,Var,dat,met, colorby = "alignable", shapeby = "AMP_Date")#, Colors = c("red","orange","blue"))#Colors = c("#00c7e4","#6ca425","#a800b3","#e19041"))
#dev.off()
#or with out a gene
PC2D(dat,met)
#PCA 1D
component <- 2
group <- "GFAP"
PC1D(scores,dat,met,group,component)
scores2 <- cbind(scores, met)
anova(lm(PC1~alignable,scores2 ))
#PCA Score by gene
PCGene(dat,met,"FOS",component)
#T-SNE
a <- Tsne(R,dat,met,colorby="PROX1",gene="Arc")
a[1]
a[2]

# Plot Single Gene --------------------------------------------------------
samples <- metaProxC[metaProxC$alignable >  100000 &  metaProxC$Smartseq2_RT_enzyme_used == "ProtoscriptII" ,"Sample_ID"]#
#metaProxC$CTIP2 == "N" & metaProxC$PROX1 == "N" & metaProxC$FOS == "N" & metaProxC$Mouse_condition == "HC" & metaProxC$alignable >  500000 & metaProxC$Smartseq2_RT_enzyme_used == "ProtoscriptII"  ,"Sample_ID"]
dat <- tpmProxC[, samples]
met <- metaProxC[match(samples,metaProxC$Sample_ID),]
met$Mouse_condition <- as.character(met$Mouse_condition)
met[met$Mouse_condition == "EE","Mouse_condition"] <- "NE"
met$Brain_Region <- as.character(met$Brain_Region)
met[met$Brain_Region == "HDG", "Brain_Region"] <- "pIN"
met[met$Brain_Region == "CA3_other_negs", "Brain_Region"] <- "Neg"
#met[as.numeric(dat["Gad2",]) > 1 & met$Brain_Region == "Neg","Brain_Region"] <- "IN"
met$Brain_Region <- factor(met$Brain_Region, levels = c("CA1","Neg","pIN","DG"))
#tiff(filename = "~/Documents/SalkProjects/ME/ShortLongSingature/SLSig_tiff/gene.tiff",width = 6,height = 3,units = 'in',res = 300)
Indiv("Tnik",dat, met)
          #dev.off()
IndivSubgroup("Ifi203",dat, met)

IndivByDate("Prox1",dat, met)
IndivProx1Grouped("Fos")
#plot two genes
a <- "Grin3a"
b <- "Fos"
group <- "fos"
Plot2Genes(a,b, dat,met)
res <- res.HC_N_P_1
Volcano(res)

######
samples <- metaProxC[ metaProxC$Brain_Region == "DG" & metaProxC$Mouse_condition == "EE" & metaProxC$alignable >  100000 & metaProxC$Smartseq2_RT_enzyme_used == "ProtoscriptII" ,"Sample_ID"]#
#metaProxC$CTIP2 == "N" & metaProxC$PROX1 == "N" & metaProxC$FOS == "N" & metaProxC$Mouse_condition == "HC" & metaProxC$alignable >  500000 & metaProxC$Smartseq2_RT_enzyme_used == "ProtoscriptII"  ,"Sample_ID"]
tmp <- dat <- tpmProxC[, samples]
met <- metaProxC[samples,]
colnames(tmp) <- paste( met$Brain_Region,met$FOS, round(as.numeric(dat["Grin2a",])))#, met$FOS,c(1:ncol(dat)),sep = ".")
#rownames(tmp) <- toupper((rownames(tmp)))
#genes <- c("Man1a","Bcl11b","Slc6a1","Arpp21","Col15a1","Bok","Sst","Camk2a","Gad1","Gad2","Pvalb")
#genes <- c("Foxg1","Wnt5a","Dcx","Prox1","Rbfox3","Camk2a","Creb1","Neurod1","Sox11")
#upstream <- c("Creb1","Crebbp","Grin1","Grin2a","Grin2b","Gria1","Gria2","Gria3","Gria4","Gabra1","Gabra2","Gabrb","Cacna1a","Cacna1b","Cacna1c","Cacnai","Mapk3","Mapk1","Elk1","Srf","Rps6ka3")
#neg <- c("Sostdc1","Ttr","Wfs1","Pantr1","C1ql2","Pvalb","Reln","Map3k15","Sst","Gad1","Cdh24","Mpped1")
#genes <- c("Ppp1cc","Ppp1cb","Ppp1ca","Per1","Fos","Bdnf","Atf1","Creb1","Crebbp","Kcnip3","Carf")
genes <- c("Crebbp","Htr1b","Fos","Arc")
tmp2 <- tmp[genes,]
#p <- apply(X = tmp2, MARGIN = 1, FUN = rawExp)
tiff(filename = "~/Documents/SalkProjects/ME/ShortLongSingature/SLSig_tiff/test.tiff",width = 12,height = 12,units = 'in',res = 300)
h <- heatmap(as.matrix(na.exclude(tmp2)),scale = "row")
dev.off()
#genes <- celltypegenes
#a <- apply(X = dat[genes,], MARGIN = 1,FUN = rawExp, i = 4)
#genes <- table(c(genes, names(a[a>4])))
#genes <- names(genes[genes == 2])
#o1 <- data.frame(a  = met$Brain_Region, b = cutree(hclust(dist(t(dat[celltypegenes,]))),k=40))
#sampleorder <- o1$b 
#sampleorder[o1$a == "DG"] <- sampleorder[o1$a == "DG"] * 1000
#sampleorder[o1$a == "pIN"] <- sampleorder[o1$a == "pIN"] * 80
#sampleorder[o1$a == "IN"] <- sampleorder[o1$a == "IN"] * 20
#sampleorder[o1$a == "CA1"] <- sampleorder[o1$a == "CA1"] * 1
#tiff(filename = "~/Documents/SalkProjects/ME/ShortLongSingature/SLSig_tiff/celltypes_heat.tiff",width = 8,height = 12,units = 'in',res = 300)
sampleorder <- c(h$colInd)
geneorder <- c("Crebbp","Creb1","Fos","Arc","Mapk1")
heatMe(dat,met,genes,sampleorder = sampleorder,geneorder = c(1,2,3,4,5) , cutoff = 3)
#dev.off()
######
p <- apply(X=tpmProx[,metaProx$prox == "P"],MARGIN=1,FUN=propExp)
m <- apply(X=tpmProx[,metaProx$prox == "P"],MARGIN=1,FUN=meanNoZero)
plot(p,m)
rownames(tpmProx[m > 10 & p < 0.05,])



#####
samples <- metaProxC[metaProxC$Mouse_condition == "EE" & metaProxC$alignable >  100000 & metaProxC$Smartseq2_RT_enzyme_used == "ProtoscriptII" ,"Sample_ID"]
#metaProxC$Brain_Region == "CA3_other_negs" & metaProxC$Mouse_condition == "HC" & metaProxC$alignable >  100000 & metaProxC$Smartseq2_RT_enzyme_used == "ProtoscriptII","Sample_ID"]
dat.1 <- na.exclude(countProxC[, samples])
dat <- na.exclude(tpmProxC[, samples])
met <- metaProxC[match(samples,metaProxC$Sample_ID),]

df <- melt(t(dat))
df$IN <- met$IN
df$Brain_Region <- met$Brain_Region
df$Brain_Region <- as.character(df$Brain_Region)
df$FOS <- met$FOS

df[df$Brain_Region == "CA3_other_negs","Brain_Region"] <- "Neg"
df[df$Brain_Region == "HDG","Brain_Region"] <- "P+Neg"


gene <- "Serpine2"


ggplot(df[as.character(df$X2) == gene,], aes(as.factor(IN), value, colour = FOS))+
  geom_violin()+
  geom_point(data = df[as.character(df$X2) == gene,], aes(as.factor(IN), value, colour = Brain_Region))+
  labs(title  = gene)+
  theme_bw(base_size = 16)+
  facet_grid(~FOS)#+
  #scale_colour_manual(values = c("black","#00c7e4","#a800b3","#6ca425","#e19041"))
