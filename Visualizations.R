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
#save(list = c("t.all","t.hc","t.hc.n","t.hc.neg"),file = "~/Documents/SalkProjects/ME/ShortLongSingature/SLSig_R/tsne.rda",compress = TRUE)
#save(list = c("t.all"),file = "~/Documents/SalkProjects/ME/ShortLongSingature/SLSig_R/tsne2.rda",compress = TRUE)
#load(c("~/Documents/SalkProjects/ME/ShortLongSingature/SLSig_R/tsne2.rda"))
#08/10/2016 changes
#save(list = c("t.hc.maingroups","t.fosN.maingroups","t.all","t.EE"),file = "~/Documents/SalkProjects/ME/ShortLongSingature/SLSig_R/tsne3.rda",compress = TRUE)
#load("~/Documents/SalkProjects/ME/ShortLongSingature/SLSig_R/tsne3.rda")
#10/23/2016: same as above, but removed outliers, both from alignment depth and cluster outliers
#save(list = c("t.all2","t.gaba","t.cge"), file = "~/Documents/SalkProjects/ME/ShortLongSingature/SLSig_R/tsne4.rda",compress = TRUE)
#load("~/Documents/SalkProjects/ME/ShortLongSingature/SLSig_R/tsne4.rda")
###############################################
## FUNCTIONS THAT WILL PLOT FOR YOU
###############################################

propExp <- function(x,i=1){
  if(sum(is.na(x)) > i){
    sum(!is.na(x))/length(x)
  }else{
    sum(x > i)/length(x)
  }
}
rawExp <- function(x,i=1){
  #if(sum(is.na(x)) > i){
  #  sum(!is.na(x))
  #}else{
  sum(na.exclude(x) > i)
  #}
}
meanNoZero <- function(x,i=0){
  a <- mean(x[x>i])
  if(sum(is.na(x)) > i){
    a <- mean(x[!is.na(x)])
  }else if(sum(x>i) == 0){
    a <- 0
  }
  return(a)
}
varNoZero <- function(x){
  sd(x[x >0])
}
MedNoZero <- function(x){
  median(x[x >0])
}
countExp <- function(x){
  sum(x > 0)
}
scaleME <- function(x){
  return(scale(x)[,1])
}
elbow <- function(x){
  x <- as.vector(scale(x[order(x)]))
  x.1 <- c(0,diff(x))
  elbow <- max(x.1)
  cutoff <- mean(c(x[which(x.1 == elbow)], x[which(x.1 == elbow) -1]))
  return(cutoff)
}
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
      scale_colour_gradient(high="red",low="grey")+
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
# Inidividual Genes (tpm) -------------------------------------------------
Indiv <- function(gene,dat,met){
  
  tmp <- data.frame(exp = as.numeric(dat[gene,]),
                    Arc = as.numeric(dat["Arc",]))
  tmp <- cbind(tmp, met)
  tmp$FOS <- as.character(tmp$FOS)
  tmp[tmp$fos == "F","FOS"] <- "High"
  tmp[tmp$fos == "L","FOS"] <- "Low"
  tmp[tmp$fos == "N","FOS"] <- "None"
  tmp$FOS <- factor(tmp$FOS, levels = c("N","L","F"))
  tmp$Brain_Region <- factor(tmp$Brain_Region, c("DG","CA1","VIP","IN"))
  #pdf("~/Documents/SalkProjects/BenLacar/ManuscriptFigures/Camk4.pdf",width=6,height=5)
  p <- ggplot(tmp, aes(FOS,exp,fill = FOS, alpha = 0.3))+
    geom_violin()+#outlier.shape=NA)+
    geom_point(position=position_jitter(width=0.01,height=0),  shape = 1, size = 0.5)+
    theme_bw(base_size=20)+
    ylab("TPM")+
    scale_fill_manual(values = c("blue","red"))+
    scale_colour_gradient(high="red", low="blue")+
    ylim(c(0,max(tmp$exp) + 0.2*(max(tmp$exp))))+
    theme(panel.border = element_rect(colour=c("black"),size=2),
          axis.ticks = element_line(size=1.5))+
    labs(title=paste(gene,"\n"))+
    facet_grid( Mouse_condition   ~  Subgroup2 ) 
return(p)
}
Indiv2 <- function(gene,dat,met,group){
  
  tmp <- data.frame(exp = as.numeric(dat[gene,]),
                    Arc = as.numeric(dat["Arc",]))
  tmp <- cbind(tmp, met)
  tmp$FOS <- as.character(tmp$FOS)
  tmp[tmp$fos == "F","FOS"] <- "High"
  tmp[tmp$fos == "L","FOS"] <- "Low"
  tmp[tmp$fos == "N","FOS"] <- "None"
  tmp$FOS <- factor(tmp$FOS, levels = c("N","L","F"))
  tmp$Brain_Region <- factor(tmp$Brain_Region, c("DG","CA1","VIP","IN"))
  tmp[,group] <- met[,group]
  colnames(tmp)[colnames(tmp) == group]  <- "group"
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
    facet_grid( Mouse_condition   ~  group  ) 
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
IndivProp <- function(gene,dat,met){
  
  tmp <- data.frame(exp = as.numeric(dat[gene,]),
                    Arc = as.numeric(dat["Arc",]))
  tmp <- cbind(tmp, met)
  tmp$FOS <- as.character(tmp$FOS)
  tmp[tmp$fos == "F","FOS"] <- "High"
  tmp[tmp$fos == "L","FOS"] <- "Low"
  tmp[tmp$fos == "N","FOS"] <- "None"
  tmp$FOS <- factor(tmp$FOS, levels = c("N","L","F"))
  tmp$Brain_Region <- factor(tmp$Brain_Region, c("DG","CA1","VIP","IN"))
  
  a <- vector()
  for (k in unique(tmp$Subgroup2)){
  for (i in unique(tmp$Mouse_condition)){
    for (j in unique(tmp$FOS)){
    a <- rbind(a, c(propExp(tmp[tmp$Mouse_condition == i & tmp$FOS == j, "exp"], i = 1), i, j,k))
    }
  }
  }
  a2 <- as.data.frame(a)
  a2$V1 <- as.numeric(as.character(a2$V1))
  colnames(a2) <- c("exp","Mouse_condition","FOS","Subgroup2")
  a2$Mouse_condition <- factor(a2$Mouse_condition, levels(tmp$Mouse_condition))
  #pdf("~/Documents/SalkProjects/BenLacar/ManuscriptFigures/Camk4.pdf",width=6,height=5)
  p <- ggplot(a2, aes(Mouse_condition, exp))+
    geom_bar(stat = "identity")+
    theme_bw(base_size=20)+
    ylab("")+
    theme(panel.border = element_rect(colour=c("black"),size=2),
          axis.ticks = element_line(size=1.5))+
    labs(title=paste(gene,"\n"))+
    facet_grid(FOS~Subgroup2)+
    ylim(c(0,1))
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
Plot3D.TSNE <- function(t,group = "Brain_Region",COLORS = c("#00c7e4","#6ca425","#a800b3","#e19041")){
  require(scatterplot3d)
  #source('http://www.sthda.com/sthda/RDoc/functions/addgrids3d.r')
  t$COL <- as.character(t[,group])
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
heatMe <- function(dat,met,genes,k1= NULL , k2 = NULL, sampleorder = NULL,geneorder = NULL,  cutoff = NULL,samplenames=NULL){
  tmp <- na.exclude(dat[genes,])
  if(!is.null(samplenames)){
    colnames(tmp) <- samplenames
  }else{
    colnames(tmp) <- paste(met$PROX1,met$CTIP2,met$FOS,met$Mouse_condition,c(1:nrow(met)),sep = ".")
  }
  tmp.1 <- t(scale(t(tmp)))
  tmp2 <- melt(t(na.exclude(tmp.1)))
  if (!is.null(k1)){
    o1 <- cutree(hclust(dist(t(na.exclude(tmp.1)))),k=k1)
  }else{
    o1 <- sampleorder
  }
  if(!is.null(geneorder)){
    o2 <- geneorder
  }
  if(!is.null(k2)){
    o2 <- cutree(hclust(dist(na.exclude(tmp.1))),k=k2)
  }else{
    o2 <- 1:nrow(na.exclude(tmp.1))
  }
  tmp2$o1 <- as.numeric(o1)
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
heatMeRaw <- function(dat,met,genes,k1= NULL , k2 = NULL, sampleorder = NULL,geneorder = NULL,  cutoff = NULL,samplenames=NULL){
  tmp <- na.exclude(dat[genes,])
  if(!is.null(samplenames)){
    colnames(tmp) <- samplenames
  }else{
    colnames(tmp) <- paste(met$PROX1,met$CTIP2,met$FOS,met$Mouse_condition,c(1:nrow(met)),sep = ".")
  }
  tmp.1 <- tmp
  tmp2 <- melt(t(na.exclude(tmp.1)))
  if (!is.null(k1)){
    o1 <- cutree(hclust(dist(t(na.exclude(tmp.1)))),k=k1)
  }else{
    o1 <- sampleorder
  }
  if(!is.null(geneorder)){
    o2 <- geneorder
  }
  else if(!is.null(k2)){
    o2 <- cutree(hclust(dist(na.exclude(tmp.1))),k=k2)
  }else{
    o2 <- 1:nrow(na.exclude(tmp.1))
  }
  tmp2$o1 <- as.numeric(o1)
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
    theme(axis.text.x = element_text(angle = 90, hjust = 1), 
          text = element_text(size = 20))
  return(p1)
}
heatMeAvg <- function(dat,met,genes,group = "Subgroup2", k1= NULL , k2 = NULL, sampleorder = NULL,geneorder = NULL,  cutoff = NULL,samplenames=NULL){
  tmp <- na.exclude(dat[genes,])
  if(!is.null(samplenames)){
    colnames(tmp) <- samplenames
  }else{
    colnames(tmp) <- paste(met$PROX1,met$CTIP2,met$FOS,met$Mouse_condition,c(1:nrow(met)),sep = ".")
  }
  tmp.1 <- na.exclude(t(scale(t(tmp))))
  tmp3.m <- vector()
  tmp3.p <- vector()
  nm <- vector()
  for(g in unique(met[,group])){
    if (!is.null(dim(tmp.1[,met[,group] == g])[1])){
      #mean
      a <- apply(tmp.1[,met[,group] == g],1,mean)
      tmp3.m <- cbind(tmp3.m, as.numeric(a))
      #prop
      a <- apply(tmp.1[,met[,group] == g],1,propExp)
      tmp3.p <- cbind(tmp3.p, as.numeric(a))
      nm <- c(nm, g)
    }
  }
  colnames(tmp3.m) <- colnames(tmp3.p) <- nm
  rownames(tmp3.m) <- rownames(tmp3.p) <- rownames(tmp.1)
  tmp3.m <- tmp3.m
  tmp3.p <- tmp3.p
  tmp4 <- rbind(melt(t(tmp3.m)), melt(t(tmp3.p)))
  tmp4[,"analysis"] <- c(rep("mean",nrow(melt(t(tmp3.m)))),rep("prop",nrow(melt(t(tmp3.m)))))
  if (!is.null(k1)){
    o1 <- cutree(hclust(dist(t(tmp3.m))),k=k1)
  }else if (!is.null(sampleorder)){
    o1 <- sampleorder
  }else{
    o1 <- 1:ncol(tmp3.m)
  }
  if(!is.null(geneorder)){
    o2 <- geneorder
  }else if(!is.null(k2)){
    o2 <- cutree(hclust(dist(na.exclude(tmp3.m))),k=k2)
  }else{
    o2 <- 1:nrow(na.exclude(tmp3.m))
  }
  tmp4$o1 <- as.numeric(o1)
  tmp4$o2 <- rep(as.numeric(o2), each = ncol(tmp3.m))
  tmp4$value <- as.numeric(as.character(tmp4$value))
  if(!is.null(cutoff)){
    tmp4[tmp4$value > cutoff, "value"] <- cutoff + 0.1
    tmp4[tmp4$value < (-1 * cutoff), "value"] <- (-1 * cutoff) - 0.1
  }
  p1 <- ggplot(tmp4[tmp4$analysis == "mean",], aes(reorder(X1,o1), reorder(X2,o2) , fill = value))+
    geom_tile()+
    theme(axis.title.x=element_blank(),
          axis.title.y=element_blank()
    )+
    scale_fill_gradient2(high="red",mid = "white",low="blue")+
    #theme_bw(base_size=22)+
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    theme(text = element_text(size = 25, colour = "black"))
    #facet_grid(~analysis,scales = "free")
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
samples <- rownames(metaProxC[ metaProxC$FOS != "L"   & metaProxC$Context1 != "none"   & metaProxC$outliers == "in"  ,])
                                             
dat <- na.exclude(tpmProxC[, samples])
met <- metaProxC[samples,]

p <- pca(t(dat[,samples]),nPcs = 2)
scores <- as.data.frame(p@scores)
scores <- cbind(scores,met)
loading <- as.data.frame(p@loadings)
loading <- loading[order(loading$PC2,decreasing=TRUE),]
Var <- p@R2
met$Brain_Region <- as.character(met$Brain_Region)
met[met$Brain_Region == "HDG", "Brain_Region"] <- "VIP"
met$group <- paste(met$FOS,  met$Mouse_condition, sep =".")
met$Vip <- as.numeric(dat["Vip",] > 7)
#tiff(filename = "~/Documents/SalkProjects/ME/ShortLongSingature/SLSig_tiff/PCA_HC_N.tiff",width = 6.5,height = 5,units = 'in',res = 300)
PC2D(scores,Var,dat,met,colorby = "alignable", shapeby = "Mouse_condition")#,Colors = c("red","blue","black"))# c("#00c7e4","#6ca425","#a800b3","#e19041"))
#dev.off()
#or with out a gene
PC2D(dat,met)
#PCA 1D
component <- 1
group <- "group"

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
samples <- rownames(metaProxC[ metaProxC$Context2 == "none" & metaProxC$Brain_Region == "DG"  &  metaProxC$FOS != "L" &  metaProxC$outliers == "in"  & metaProxC$Arc_2.5 != "greater" ,])#
#metaProxC$CTIP2 == "N" & metaProxC$PROX1 == "N" & metaProxC$FOS == "N" & metaProxC$Mouse_condition == "HC" & metaProxC$alignable >  500000 & metaProxC$Smartseq2_RT_enzyme_used == "ProtoscriptII"  ,"Sample_ID"]
dat <- na.exclude(tpmProxC[, samples])
met <- metaProxC[samples,]
met$Mouse_condition <- as.character(met$Mouse_condition)
met[met$Mouse_condition == "EE","Mouse_condition"] <- "1hr"
met[met$Mouse_condition == "5hpA","Mouse_condition"] <- "5hr"
met[met$Mouse_condition == "5hpAA","Mouse_condition"] <- "A>A"
met[met$Mouse_condition == "5hpAC","Mouse_condition"] <- "A>C"
met$Mouse_condition <- factor(x = met$Mouse_condition,levels = c("HC","1hr","5hr","A>A","A>C"))
met$Brain_Region <- as.character(met$Brain_Region)
met[met$Brain_Region == "HDG","Brain_Region"] <- "VIP"
met$Subgroup2 <- factor(met$Subgroup2, c("DG","CA3","CA1","Sub","GC","VIP","Pvalb","Lamp5"))
#tiff(filename = "~/Documents/SalkProjects/ME/ShortLongSingature/MolecDissec_Figs_Tables/Figures_vD/Gla.tiff",width = 25,height = 5,units = 'in',res = 300)#single gene = 8 x 3.5, hc and ne 8 x 5
Indiv("Smad7",dat, met)
#dev.off()
group <- "Mouse"
Indiv2("Dynlt3",dat, met,group )
IndivProp("Tmem170",dat, met)



#
a <- "Sox11"#Bap1
b <- "Dgat2l6"
group <- "Activations"
Plot2Genes(a,b, dat,met,group)
res <- res.HC_N_P_1
Volcano(res)


###########
samples <- rownames(metaProxC[ metaProxC$Mouse_condition == "HC" & metaProxC$Brain_Region != "CA3_other_negs" & metaProxC$FOS != "L"  &  metaProxC$Arc_2.5 != "greater" & metaProxC$cluster_outlier == "in" & metaProxC$Context1 == "none" & metaProxC$outliers == "in" ,
                       ])#

tmp <- na.exclude(tpmProxC[, samples])
met <- metaProxC[match(samples,metaProxC$Sample_ID),]
colnames(tmp) <- paste( met$Mouse_condition,met$FOS)
tmp2 <- tmp[genes,]
#tiff(filename = "~/Documents/SalkProjects/ME/ShortLongSingature/SLSig_tiff/test.tiff",width = 12,height = 12,units = 'in',res = 300)
h <- heatmap(as.matrix(na.exclude(tmp2)),scale = "col")
#dev.off()
tiff(filename = "~/Documents/SalkProjects/ME/ShortLongSingature/MolecDissec_Figs_Tables/Figures_vD/K_heat.tiff",width = 8,height = 12,units = 'in',res = 300)
heatMeRaw(dat,met,genes,k1 = 10,k2 = 5, samplenames = paste(met$Brain_Region,1:ncol(dat),sep="."))
dev.off()
heatMe(dat,met,genes,k1 = 10, k2 = 10 , cutoff = 3,samplenames = paste(met$FOS,met$Mouse_condition,1:ncol(dat),sep="."))

#tiff(filename = "~/Documents/SalkProjects/ME/ShortLongSingature/MolecDissec_Figs_Tables/Figures_vC/Fig6_potassium_hc.tiff",width = 8,height = 15,units = 'in',res = 300)
heatMeAvg(dat,met,genes,k2 = 8,k1 = 4 ,cutoff = 1 )
#dev.off()


###########
## T-sne
###########
samples <- rownames(metaProxC[metaProxC$Mouse_condition == "EE" & metaProxC$Brain_Region == "DG"  & metaProxC$FOS == "F"  & metaProxC$cluster_outlier == "in" & metaProxC$outliers == "in"|
                                metaProxC$Mouse_condition == "HC" & metaProxC$Brain_Region == "DG"  & metaProxC$FOS == "N"  & metaProxC$cluster_outlier == "in" & metaProxC$outliers == "in"|
                                metaProxC$Mouse_condition == "5hpA" & metaProxC$Brain_Region == "DG"  & metaProxC$FOS == "N"  & metaProxC$cluster_outlier == "in" & metaProxC$outliers == "in"
                                ,])
dat <- na.exclude(tpmProxC[, samples])
met <- metaProxC[samples,]


i <- 30#17
TSNE <- Rtsne(as.matrix(t(na.exclude(dat))),initial_dims=2,perplexity=i,theta=0,check_duplicates=FALSE,dims = 2)
#df <- rbind(df, TSNE$Y[,1], TSNE$Y[,2])

t <- as.data.frame(TSNE$Y)
colnames(t) <- c("T1","T2")
t <- cbind(t,met)
#metaProxC$Mouse3 <- 0
#for (i in unique(t$Brain_Region)){
#  a <- unique(t[t$Brain_Region == i, "Mouse2"])
#  for(j in 1:length(a)){
#    t[t$Brain_Region == i & t$Mouse2 == a[j], "Mouse3"] <- j
#  }
#}
#t$group <- paste(t$Brain_Region, t$Mouse3, sep =".")
t$gene <- as.numeric(tpmProxC["Slc4a5",rownames(t)])
#k <- kmeans(t[,c(1:2)],centers =2,nstart = 2000)
#t$k <- as.factor(k$cluster)
#tiff("~/Documents/SalkProjects/ME/ShortLongSingature/MolecDissec_Figs_Tables/Figures_vD/tsne_bcl11b.tiff",width = 9,height = 6.5,units = 'in',res = 600,compression = 'lzw')
#t$a <- a
t$Mouse_condition <- as.character(t$Mouse_condition)
t[t$Mouse_condition == "5hpA","Mouse_condition"] <- "5hr"
t[t$Mouse_condition == "EE","Mouse_condition"] <- "1hr"
t$Mouse_condition <- factor(t$Mouse_condition,c("HC","1hr","5hr"))

tiff(filename = "~/Documents/SalkProjects/ME/ShortLongSingature/Figs/5hr_EE.tiff",width = 10,height = 7,units = 'in',res = 500,compression = "lzw")
ggplot(t, aes(T1,T2, color = Mouse_condition, shape = FOS))+
  geom_point(size = 5)+
  theme_bw()+
  xlab("TSNE1")+
  ylab("TSNE2")+
  theme(text=element_text(size=20))+
  theme(panel.border = element_rect(colour=c("black"),size=2),
        axis.ticks = element_line(size=1.5),
        panel.grid.major = element_line(size = 1))+
  scale_colour_manual(values = c("black","red","blue"))
  #scale_shape_manual(values=c(8, 14, 16, 15, 17))+
#  scale_colour_manual(values = c("darkblue","orange","black"))+
   # scale_alpha_continuous(range  = c(0.6,1))+
#scale_colour_gradient2(high = "red",low = "black",mid = "grey", midpoint = 1)#+
#scale_colour_manual(values = c("#94bc68","#4b7023","#30510b","#345510","#60d6eb","#13aac5","#e9ae79","#d29258","#d97923","#7e5530","#c05cc7","#e93af5","#83108b"))
#scale_colour_manual(values = c("#6ca425","#00c7e4","#e19041","#a800b3"))
dev.off()

########
samples <- metaProxC[ metaProxC$Mouse_condition == "HC" & metaProxC$Context1 == "none" & metaProxC$FOS != "L" & metaProxC$outliers == "in" ,
                      "Sample_ID"]
dat <- na.exclude(tpmProxC[, samples])
met <- metaProxC[match(samples,metaProxC$Sample_ID),]


Cor <- cor(dat[genes,])
Cor2 <- as.data.frame(melt(t(Cor)))
Cor2$Subgroup2 <- met$Subgroup2
Cor2$FOS <- met$FOS
Cor2$group1 <- paste(met$Subgroup2, met$FOS,c(1:ncol(dat)), sep = ".")
Cor2$group2 <- rep(paste(met$Subgroup2, met$FOS,c(1:ncol(dat)) ,sep = "."), each = ncol(dat))

#tiff(filename = "~/Documents/SalkProjects/ME/ShortLongSingature/MolecDissec_Figs_Tables/Figures_vC/FOS_heat.tiff",width = 10,height = 8,units = 'in',res = 300)
ggplot(Cor2, aes(group1, group2, fill = value))+
  geom_tile()+
  scale_fill_gradient2(low = "blue" , high = "red",mid = "white",midpoint = 0.5)
#dev.off()


#####
samples <- rownames(metaProxC[ metaProxC$Context2 != "none" & metaProxC$Brain_Region == "DG" & metaProxC$FOS == "F" & metaProxC$cluster_outlier == "in" & metaProxC$outliers == "in" ,])
dat <- na.exclude(tpmProxC[, samples])
met <- metaProxC[samples,]
tmp <- data.frame(reexposure = met$Context2,sample = samples,  gene = as.numeric(dat["Dgat2l6",]))

ggplot(tmp, aes(gene, fill = reexposure))+
  geom_density(alpha = 0.5)+
  scale_fill_manual(values= c("black","blue"))
