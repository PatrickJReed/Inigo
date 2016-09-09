#save(list = c("SigGenes"), file = "~/Documents/SalkProjects/ME/ShortLongSingature/SLSig_R/versHC.rda",compress = TRUE)
#load("~/Documents/SalkProjects/ME/ShortLongSingature/SLSig_R/versHC.rda")

#res2 is EE FOS P versus FOS N
#res is EE FOS P versus HC FOS N
SigGenes <- list()
i <- 0
for (g2 in c("DG","CA1","CA3","VIP","Neg")){
res2 <- RES.activity[[g2]]
res <- RES.FosF.HCN[[g2]]

a <- cbind(res2, res[rownames(res2),])
a <- na.exclude(a)

#1 Fos F low genes
genes1 <- rownames(a[(a[,1] < 0  & #FOS F < FOS N
       a[,4] < 0.05 & a[,5] > 0.2 &
       a[,9] < 0 &                  #FOS F < HC N 
       a[,12] < 0.05
      ),])
#2 Fos F high genes
genes2 <- rownames(a[(a[,1] > 0  & #FOS F > FOS N
      a[,4] < 0.05 & a[,5] > 0.2 & #FOS F > HC N
      a[,9] > 0 & 
      a[,12] < 0.05
),])

#res2 is EE FOS P versus FOS N
#res is EE FOS N versus HC FOS N
res2 <- RES.activity[[g2]]
res <- RES.FosN.HCN[[g2]]
a <- cbind(res2, res[rownames(res2),])
a <- na.exclude(a)


#3 FOS L low genes
genes3 <- rownames(a[(a[,1] > 0  & #FOS N < FOS P
      a[,4] < 0.05 & a[,5] > 0.2 &
      a[,9] < 0 &                   #FOS N < HC N 
      a[,12] < 0.05
),])
#4 FOS L high genes
genes4 <- rownames(a[(a[,1] < 0  & #FOS N > FOS P
      a[,4] < 0.05 & a[,5] > 0.2 &
      a[,9] > 0 &                   #FOS N > HC N 
      a[,12] < 0.05
),])


#SigGenes <- list()
nm <- paste(g2,".FOSplow", sep = "")
i <- i + 1
SigGenes[[i]] <- genes1
names(SigGenes)[i] <- nm

nm <- paste(g2,".FOSphigh", sep = "")
i <- i + 1
SigGenes[[i]] <- genes2
names(SigGenes)[i] <- nm

nm <- paste(g2,".FOSllow", sep = "")
i <- i + 1
SigGenes[[i]] <- genes3
names(SigGenes)[i] <- nm

nm <- paste(g2,".FOSlhigh", sep = "")
i <- i + 1
SigGenes[[i]] <- genes4
names(SigGenes)[i] <- nm
}




###################
# Combine into one list
###################
FOS.F.LOW <- vector()
for (i in seq(1,length(SigGenes),4)){
  FOS.F.LOW <- c(FOS.F.LOW, SigGenes[[i]])
}
genes <- FOS.F.LOW <- unique(FOS.F.LOW)
a2 <- data.frame(row.names = genes)
for (j in c("DG","CA1","CA3","VIP","Neg")){
  nm <- paste(j,"FOSplow",sep ="." )
  a2[SigGenes[[nm]],nm] <- 1
  a2[is.na(a2[,nm]),nm] <- 0
}
a2$AllCells <- rowSums(a2)

a2$Cell <- NA
for (i in 1:nrow(a2)){
  if(a2[i,"AllCells"] == 1){
    a2[i,"Cell"] <- do.call("rbind",strsplit(colnames(a2)[which(a2[i,c(1:(ncol(a2)-2))] == 1)],split = ".",fixed = TRUE))[,1]
  }else{
    "none"
  }
}

#####
FOS.F.HIGH <- vector()
for (i in seq(2,length(SigGenes),4)){
  FOS.F.HIGH <- c(FOS.F.HIGH, SigGenes[[i]])
}
FOS.F.HIGH <- unique(FOS.F.HIGH)
genes <- FOS.F.HIGH <- unique(FOS.F.HIGH)
a2 <- data.frame(row.names = genes)
for (j in c("DG","CA1","CA3","VIP","Neg")){
  nm <- paste(j,"FOSphigh",sep ="." )
  a2[SigGenes[[nm]],nm] <- 1
  a2[is.na(a2[,nm]),nm] <- 0
}
a2$AllCells <- rowSums(a2)
a2$Cell <- NA
for (i in 1:nrow(a2)){
  if(a2[i,"AllCells"] == 1){
    a2[i,"Cell"] <- do.call("rbind",strsplit(colnames(a2)[which(a2[i,c(1:(ncol(a2)-2))] == 1)],split = ".",fixed = TRUE))[,1]
  }else{
    "none"
  }
}


FOS.N.LOW <- vector()
for (i in seq(3,length(SigGenes),4)){
  FOS.N.LOW <- c(FOS.N.LOW, SigGenes[[i]])
}
genes <- FOS.N.LOW <- unique(FOS.N.LOW)
a2 <- data.frame (row.names = genes)
for (j in c("DG","CA1","CA3","VIP","Neg")){
  nm <- paste(j,"FOSllow",sep ="." )
  a2[SigGenes[[nm]],nm] <- 1
  a2[is.na(a2[,nm]),nm] <- 0
}
a2$AllCells <- rowSums(a2)
a2$Cell <- NA
for (i in 1:nrow(a2)){
  if(a2[i,"AllCells"] == 1){
    a2[i,"Cell"] <- do.call("rbind",strsplit(colnames(a2)[which(a2[i,c(1:(ncol(a2)-2))] == 1)],split = ".",fixed = TRUE))[,1]
  }else{
    "none"
  }
}


FOS.N.HIGH <- vector()
for (i in seq(4,length(SigGenes),4)){
  FOS.N.HIGH <- c(FOS.N.HIGH, SigGenes[[i]])
}
genes <- FOS.N.HIGH <- unique(FOS.N.HIGH)
a2 <- data.frame(row.names = genes)
for (j in c("DG","CA1","CA3","VIP","Neg")){
  nm <- paste(j,"FOSlhigh",sep ="." )
  a2[SigGenes[[nm]],nm] <- 1
  a2[is.na(a2[,nm]),nm] <- 0
}
a2$AllCells <- rowSums(a2)
a2$Cell <- NA
for (i in 1:nrow(a2)){
  if(a2[i,"AllCells"] == 1){
    a2[i,"Cell"] <- do.call("rbind",strsplit(colnames(a2)[which(a2[i,c(1:(ncol(a2)-2))] == 1)],split = ".",fixed = TRUE))[,1]
  }else{
    "none"
  }
}
