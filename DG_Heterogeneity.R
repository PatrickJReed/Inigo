#Plot the transition through PC space
#save(list = c("t.2"),file = "~/Documents/SalkProjects/ME/ShortLongSingature/SLSig_R/dghet.rda",compress = TRUE)
#load("~/Documents/SalkProjects/ME/ShortLongSingature/SLSig_R/dghet.rda")
#t.2 = DG homecage FOS L and FOS N, k subgroups are in the ppt presentation
library(plot3D)
library(randomForest)
library(reshape)
surf3D(x = as.matrix(rnorm(100)),as.matrix(rnorm(100)),as.matrix(rnorm(100)))
###

samples <- metaProxC[metaProxC$FOS == "N" & metaProxC$Brain_Region == "DG" &  metaProxC$Mouse_condition == "HC" & metaProxC$alignable >  100000 & metaProxC$Smartseq2_RT_enzyme_used == "ProtoscriptII" ,"Sample_ID"]
dat <- na.exclude(tpmProxC[, samples])
met <- metaProxC[match(samples,metaProxC$Sample_ID),]


########
TSNE <- Rtsne(as.matrix(t(na.exclude(dat))),initial_dims=2,perplexity=i,theta=0.1,check_duplicates=FALSE)
t <- as.data.frame(TSNE$Y)
colnames(t) <- c("X","Y")
t <- cbind(t,met)
t$gene <- as.numeric(dat["Arc",])
t$Brain_Region <- as.character(t$Brain_Region)
t[t$Brain_Region == "CA3_other_negs","Brain_Region"] <- "Neg"
t[t$Brain_Region == "P+Neg","Brain_Region"] <- "pIN (P+C-)"
t$gene <- as.numeric(tpmProxC["Rbms3",samples])
k <- kmeans(t[,c(1,2)],centers = 4)
t$k <- k$cluster
ggplot(t, aes(X,Y, colour =as.factor(k), shape = FOS))+
  geom_point(size = 3)+
  theme_bw()+
  #scale_colour_gradient(high="red",low="blue")+
  #scale_colour_manual(values = c("#00c7e4","#a800b3","#6ca425","#e19041"))+
  labs(title="tSNE")+
  theme(text=element_text(size=20))+
  theme(panel.border = element_rect(colour=c("black"),size=2),
        axis.ticks = element_line(size=1.5))
  
networkgenes <- vector()
for(i in 1:4){
  for (j in 1:4){
    if(i != j){
      w <- which(k$cluster == i | k$cluster == j)
      ####
      # Round 1
      k2 <- k$cluster[w]
      tmp <- t(dat[,w])
      tmp <- as.data.frame(tmp)
      f <- randomForest(x = tmp, y = as.factor(k2))
      gini <- as.data.frame(f$importance)
      colnames(gini) <- c("MeanDecreaseGini")
      gini$mt <- c(1:nrow(gini))
      
      
      genes <- rownames(gini[gini$MeanDecreaseGini >= 0.0001,])
      
      ####
      # Round 2
  
  
      tmp <- t(dat[genes,w])
      tmp <- as.data.frame(tmp)
      f <- randomForest(x = tmp, y = as.factor(k2))
      gini <- as.data.frame(f$importance)
      colnames(gini) <- c("MeanDecreaseGini")
      gini$mt <- c(1:nrow(gini))
      genes <- rownames(gini[gini$MeanDecreaseGini >= 0.05,])
      networkgenes <- c(networkgenes, genes)
    }
  }
}
####

tmp2 <- melt(t(dat))
tmp2$k <- k$cluster

gene <- "Sgms2"
ggplot(tmp2[tmp2$X2 == gene,], aes(as.factor(k) , value))+
  geom_violin()+
  geom_point()

Cor <- cor(t(dat[networkgenes,]))
Cor2 <- melt(t(Cor))

Cor3 <- Cor2[abs(Cor2$value) > 0.6,]
write.table(x = Cor3, file = "~/Documents/test.txt",quote = FALSE)
