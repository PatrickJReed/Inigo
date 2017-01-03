## Random Forest
library(randomForest)
library(reshape)

######
qp <- 0.8

###### 
dat_P1 <- as.data.frame(t(countProx[,metaProx$cond == "HC" & metaProx$date != "150629" & metaProx$prox == "P"]))
dat_N1 <- as.data.frame(t(countProx[,metaProx$cond == "HC" & metaProx$date != "150629" & metaProx$prox == "N"]))
pp <- apply(X=dat_P,MARGIN=2,FUN=propExp)
pn <- apply(X=dat_N,MARGIN=2,FUN=propExp)
dat_P_1 <- dat_P[,pp > 0.7 & pn > 0.7 ]
dat_N_1 <- dat_N[,pp > 0.7 & pn > 0.7 ]
dat_P$group <- c("P")
dat_N$group <- c("N")
dat <- rbind(dat_P,dat_N)
colnames(dat) <- paste("X",colnames(dat),sep="")
dat$Xgroup <- as.factor(dat$Xgroup)
colnames(dat) <- gsub(pattern="-",replacement="SARA",x=colnames(dat))
dat[is.na(dat)] <- 0
#Initialize
genes <- colnames(dat)
genecount <- vector()
error <- c()
imp <- vector()
for (i in 1:2){
  rf <- randomForest(Xgroup ~ ., dat[,c(genes,"Xgroup")], ntree= 150)
  gini <- as.data.frame(rf$importance)
  #rownames(gini) <- substr(x=rownames(rf$importance),start=2,c(nchar(rownames(rf$importance))))
  colnames(gini) <- c("MeanDecreaseGini")
  gini$mt <- c(1:nrow(gini))
  
  cutoff <- quantile((gini$MeanDecreaseGini),probs=qp)
  genes <- rownames(gini[gini$MeanDecreaseGini >= cutoff,])
  if(length(which(genes == "Xgroup.1")) > 0){
    genes <- genes[-c(which(genes == "Xgroup.1"))]
  }
  error <- c(error, sum(rf$confusion[,3]))
  genecount <- c(genecount, length(genes))
  imp <- c(imp, sum(gini[genes,"MeanDecreaseGini"]))
}



#pdf("~/Documents/SalkProjects/ME/ShortLongSingature/SLSig_tiff/RandFor_HC_NP.pdf",width=7.5,height=7.5)
ggplot(gini.a, aes(resLfc, log(MeanDecreaseGini), colour = group))+
  geom_point(size=3,alpha=.5)+
  geom_point(size=3,shape=1)+
  scale_color_manual(values=c("grey","darkorange3","darkorchid3"))+
  theme_bw()+
  theme(text = element_text(size=18),panel.border = element_rect(colour=c("black"),size=2),
        axis.ticks = element_line(size=1.5),legend.position = "none")+
  labs(title="Predicitive Expression")
#dev.off()
