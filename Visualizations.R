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
###############################################
## FUNCTIONS THAT WILL PLOT FOR YOU
###############################################


# PCA 2D ---------------------------------------------------------------------
PC2D <- function(scores,Var, dat, met, colorby = NULL,shapeby = NULL, gene = NA){

  if (!is.na(gene)){
    tmp <- data.frame(PC1 = scores[,1], PC2 = scores[,2],
                      gene = as.numeric(dat[gene,]))
    tmp <- cbind(tmp,met)
    tmp$shapeby <- tmp[,shapeby]

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
    tmp$colorby <- tmp[,colorby]
    tmp$shapeby <- tmp[,shapeby]

    
    #pdf("~/Documents/SalkProjects/ME/ShortLongSingature/SLSig_tiff/PCA_HC_NP.pdf",width=8.5,height=7)
    plt <- ggplot(tmp, aes(PC1, PC2,colour = colorby,shape = shapeby))+
      geom_point(size=5, alpha=0.7)+
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
  return(plt)
}
# PCA 3D ---------------------------------------------------------------------
PC2D <- function(scores,Var, dat, met, colorby = NULL,shapeby = NULL, gene = NA){
  
  if (!is.na(gene)){
    tmp <- data.frame(PC1 = scores[,1], PC2 = scores[,2],
                      gene = as.numeric(dat[gene,]))
    tmp <- cbind(tmp,met)
    tmp$shapeby <- tmp[,shapeby]
    
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
    tmp$colorby <- tmp[,colorby]
    tmp$shapeby <- tmp[,shapeby]
    
    
    #pdf("~/Documents/SalkProjects/ME/ShortLongSingature/SLSig_tiff/PCA_HC_NP.pdf",width=8.5,height=7)
    plt <- ggplot(tmp, aes(PC1, PC2,colour = colorby,shape = shapeby))+
      geom_point(size=5, alpha=0.7)+
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



#RES <- vector()
#stdout <- vector('character')
#con <- textConnection(object='stdout',open='wr')
#sink(con)
#perps <- 12#seq(5,10,1)
#for (i in perps){
#  TSNE <- Rtsne(as.matrix(t(dat)),initial_dims=2,perplexity=i,theta=0.1,check_duplicates=FALSE)
#  res <- as.data.frame(TSNE$Y)
#  res$group <- i
#  RES <- rbind(RES,res)
#}
#sink()
#close(con)

Tsne <- function(RES,dat,met,colorby,gene="Prox1"){
  RES$prox <- met$prox
  RES$cond <- met$cond
  RES$fos <- met$fos
  RES$tmpgroup <- met$tempGroupingProx1HC
  RES$Col <- RES[,colorby]
  RES$Gene <- as.numeric(dat[gene,])
  
  Errorpos <- seq(from=27,to=length(stdout),by=27) + c(0:(length(perps)-1))
  error <- data.frame(error = as.numeric(as.character(unlist(lapply(X=Errorpos,FUN=getErrors)))),
                      perplexity = perps
  )
  
  g <- error[error$error == min(error$error),"perplexity"]
  
  p1 <- ggplot(RES[RES$group == g,], aes(V1,V2,colour = Col))+
    geom_point(size=4,alpha=1 )+
    theme_bw(base_size = 13)+
    labs(title= paste("t-SNE","\nperplexity = ",g,sep=""))+
    #scale_colour_gradient(high="black",low = "red")+
    theme(panel.border = element_rect(colour=c("black"),size=2),
          axis.ticks = element_line(size=1.5))
  
  p2 <- ggplot(as.data.frame(error), aes(perplexity,error))+
          geom_point(size=4,alpha=0.5)
  
  return(list(p1,p2))
}
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


# facet_grid( cond  ~ prox)
#dev.off()
# Correlations (2genes) ---------------------------------------------------
Plot2Genes <- function(a,b,dat,met,group){
  tmp <- data.frame(A = as.numeric(dat[a,]),
                    B = as.numeric(dat[b,]),
                    group = met[,group],
                    dusp1 = as.numeric(dat["Npas4",])
                    )
  
  #pdf("~/Documents/SalkProjects/BenLacar/ManuscriptFigures/Camk4.pdf",width=6,height=5)
  p <- ggplot(tmp, aes(A,B,colour = group))+
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

heatMe <- function(dat,genes,order){
  tmp <- dat[genes,]
  tmp.1 <- t(apply(X=tmp,MARGIN=1,FUN=scaleMe))
  tmp2 <- melt(t(tmp.1))
  tmp2$order <- rep(order, ncol(dat))
  p1 <- ggplot(tmp2, aes(X1, reorder(X2,order) , fill = value))+
    geom_tile()+
    theme(axis.title.x=element_blank(),
          axis.title.y=element_blank()
    )+
    scale_fill_gradient(high="red",low="blue")+
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
samples <- metaProxC[metaProxC$Mouse_condition == "EE" & metaProxC$alignable >  500000 & metaProxC$Smartseq2_RT_enzyme_used == "ProtoscriptII" ,"Sample_ID"]
dat <- tpmProxC[, samples]
met <- metaProxC[match(samples,metaProxC$Sample_ID),]
#gene <- "Meg3"
#Calculate the components
p <- pca(t(dat[activitygenes,]),nPcs=5)
scores <- as.data.frame(p@scores)
loading <- as.data.frame(p@loadings)
Var <- p@R2
PC2D(scores,Var,dat,met,colorby = "PROX1", shapeby = "CTIP2")

#or with out a gene
PC2D(dat,met)
#PCA 1D
component <- 2
group <- "Mouse_condition"
PC1D(scores,dat,met,group,component)
scores2 <- cbind(scores, met)
anova(lm(PC2~Mouse_condition,scores2 ))
#PCA Score by gene
PCGene(dat,met,"FOS",component)
#T-SNE
a <- Tsne(RES,dat,met,colorby="fos",gene="Arc")
a[1]
a[2]

# Plot Single Gene --------------------------------------------------------
samples <- metaProxC[metaProxC$FOS != "L"  & metaProxC$Mouse_condition == "EE" &  metaProxC$alignable >  500000 & metaProxC$Smartseq2_RT_enzyme_used == "ProtoscriptII" ,"Sample_ID"]#
#metaProxC$CTIP2 == "N" & metaProxC$PROX1 == "N" & metaProxC$FOS == "N" & metaProxC$Mouse_condition == "HC" & metaProxC$alignable >  500000 & metaProxC$Smartseq2_RT_enzyme_used == "ProtoscriptII"  ,"Sample_ID"]
dat <- tpmProxC[, samples]
met <- metaProxC[match(samples,metaProxC$Sample_ID),]

Indiv("Syt4",dat, met)
IndivByDate("Uqcr11",dat, met)
IndivProx1Grouped("Nedd8")
#plot two genes
a <- "Brf1"
b <- "Bcl11b"
group <- "fos"
Plot2Genes(a,b, dat,met,group)
res <- res.HC_N_P_1
Volcano(res)

######
samples <- metaProxC[metaProxC$Brain_Region != "DG" & metaProxC$FOS != "L"  & metaProxC$Mouse_condition == "EE" &  metaProxC$alignable >  500000 & metaProxC$Smartseq2_RT_enzyme_used == "ProtoscriptII" ,"Sample_ID"]#
#metaProxC$CTIP2 == "N" & metaProxC$PROX1 == "N" & metaProxC$FOS == "N" & metaProxC$Mouse_condition == "HC" & metaProxC$alignable >  500000 & metaProxC$Smartseq2_RT_enzyme_used == "ProtoscriptII"  ,"Sample_ID"]
tmp <- dat <- tpmProxC[, samples]
met <- metaProxC[samples,]
colnames(tmp) <- paste(met$Brain_Region, met$FOS,c(1:ncol(dat)),sep = ".")
tiff(filename = "~/Documents/SalkProjects/ME/ShortLongSingature/SLSig_tiff/activity_heat.tiff",width = 20,height = 20,units = 'in',res = 300)
heatmap(as.matrix(na.exclude(tmp[activitygenes,])),scale = "col")
dev.off()
#heatMe(dat,activitygenes,c(1:length(activitygenes)))
######
p <- apply(X=tpmProx[,metaProx$prox == "P"],MARGIN=1,FUN=propExp)
m <- apply(X=tpmProx[,metaProx$prox == "P"],MARGIN=1,FUN=meanNoZero)
plot(p,m)
rownames(tpmProx[m > 10 & p < 0.05,])
