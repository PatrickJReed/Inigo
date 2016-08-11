


nm <- table(c(rownames(RES[[1]]),
              rownames(RES[[2]]),
              rownames(RES[[3]])
              ))

nm <- names(nm[nm==3])

dg <- RES[[1]]
ca1 <- RES[[2]]
vip <- RES[[3]]

ca1 <- ca1[ca1$logFC < 0,]
dg <- dg[dg$logFC < 0,]
vip <- vip[vip$logFC < 0,]

a <- merge(x = ca1[,c(1,4)],y =  dg[,c(1,4)],by="row.names",all = TRUE)
rownames(a) <- a$Row.names
a <- a[,-c(1)]
colnames(a) <- paste(rep(c("logFC","padj"),times = 2), rep(c("ca1","dg"),each = 2),sep = ".")

a <- merge(x = a, y = vip[,c(1,4)],by="row.names",all = TRUE)
rownames(a) <- a$Row.names
a <- a[,-c(1)]
colnames(a) <- paste(rep(c("logFC","padj"),times = 2), rep(c("ca1","dg","vip"),each = 2),sep = ".")




#a <- cbind(ifelse(RES[[1]][nm,"logFC"] >0 , RES[[1]][nm,"f"], 1) ,
#           ifelse(RES[[2]][nm,"logFC"] >0 , RES[[2]][nm,"f"], 1),
#           ifelse(RES[[3]][nm,"logFC"] >0 , RES[[3]][nm,"f"], 1),
#           ifelse(RES[[4]][nm,"logFC"] >0 , RES[[4]][nm,"f"], 1),
#           ifelse(RES[[5]][nm,"logFC"] >0 , RES[[5]][nm,"f"], 1)
#           )

#colnames(a) <- paste(names(RES), "Padj",sep = "_")
#rownames(a) <- nm
a <- as.data.frame(a)
for(i in 1:ncol(a)){
  a[,i] <- as.numeric(as.character(a[,i]))
}
a2 <- a[,seq(2,6,2)]
a2$AllCellTypes <- as.numeric(unlist(apply(a2, 1,rawExp, 0.05)))
#a2$Excitatory <- as.numeric(unlist(apply(a2[,c(1,2,3)], 1,rawExp, 0.05)))
#a2$Inhibitory <- as.numeric(unlist(apply(a2[,c(4,5)], 1,rawExp, 0.05)))
a2 <- a2[order(a2$Excitatory,decreasing=TRUE),]

a2$CellSpecific <- NA
for (i in 1:nrow(a2)){
  if(a2[i,"AllCellTypes"] == 1){
    a2[i,"CellSpecific"] <- do.call("rbind",strsplit(colnames(a2)[which(a2[i,c(1:5)] < 0.05)],split = ".",fixed = TRUE))[,2]
  }else{
    "none"
  }
}

a2 <- a2[a2$AllCellTypes > 0,]
a2 <- a2[order(a2$AllCellTypes),]

write.table(a2, "~/Documents/test.txt",quote = FALSE)


###
