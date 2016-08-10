#save(list = c("SigGenes"), file = "~/Documents/SalkProjects/ME/ShortLongSingature/SLSig_R/versHC.rda",compress = TRUE)
#load("~/Documents/SalkProjects/ME/ShortLongSingature/SLSig_R/versHC.rda")

#res2 is EE FOS P versus FOS N
#res is EE FOS P versus HC FOS N

a <- cbind(res2, res[rownames(res2),])
a <- na.exclude(a)

#1
genes1 <- rownames(a[(a[,1] < 0  & 
       a[,4] < 0.05 &
       a[,7] < 0 & 
       a[,10] < 0.05
      ),])
#2
genes2 <- rownames(a[(a[,1] > 0  & 
      a[,4] < 0.05 &
      a[,7] > 0 & 
      a[,10] < 0.05
),])


a <- cbind(res2, res[rownames(res2),])
a <- na.exclude(a)


#3
genes3 <- rownames(a[(a[,1] > 0  & 
      a[,4] < 0.05 &
      a[,7] < 0 & 
      a[,10] < 0.05
),])
#4
genes4 <- rownames(a[(a[,1] < 0  & 
      a[,4] < 0.05 &
      a[,7] > 0 & 
      a[,10] < 0.05
),])

fosp <- rownames(res2[res2$logFC > 0 & res2$f < 0.01,])
tmp <- table(c(rep(fosp,4), rep(genes2,2), genes3))
genes <- names(tmp[tmp == 4])

tmp2 <- res2[genes,]
tmp2 <- tmp2[order(tmp2$f),]

#SigGenes <- list()

#SigGenes[["VIP.FOSplow"]] <- genes1
#SigGenes[["VIP.FOSphigh"]] <- genes2
#SigGenes[["VIP.FOSllow"]] <- genes3
#SigGenes[["VIP.FOSlhigh"]] <- genes4

##
a <- scan(what = "character")
b <- scan(what = "character")
vsHomecage <- scan(what = "character")
d <- scan()

df <- data.frame(a,b,vsHomecage,d)
df2 <- df[df$vsHomecage == "Total",]
ggplot(df2[order(df2$vsHomecage, decreasing=TRUE),], aes(b, d))+#, fill = vsHomecage))+
  geom_bar(stat = "identity", alpha = 0.5)+
  geom_bar(stat = "identity", color = "black")+
  facet_grid(~a)+
  scale_fill_manual(values = c("black", "grey"))+
  xlab("FOS Stain")+
  labs(title = c("Activity-induced gene expression (Total)"))+
  ylab("Gene Count")

####

a <- scan(what = "character")
b <- scan(what = "character")
d <- scan()

df <- data.frame(a,b,d)
ggplot(df, aes(b, d))+#, fill = vsHomecage))+
  geom_bar(stat = "identity", alpha = 0.5)+
  geom_bar(stat = "identity", color = "black")+
  facet_grid(~a)+
  scale_fill_manual(values = c("black", "grey"))+
  xlab("FOS Stain")+
  labs(title = c("Activity-induced gene expression\nUndetermined Differences from Homecage"))+
  ylab("Gene Count")


df <- data.frame(b ="Total", a,d)
ggplot(df, aes(b,d))+#, fill = vsHomecage))+
  geom_bar(stat = "identity", alpha = 0.5)+
  geom_bar(stat = "identity", color = "black")+
  facet_grid(~a)+
  scale_fill_manual(values = c("black", "grey"))+
  xlab("Cell Type")+
  labs(title = c("Activity-induced gene expression\nAll Genes"))+
  ylab("Gene Count")


###################
# Combine into one list
###################
FOS.F.LOW <- vector()
for (i in seq(1,length(SigGenes),4)){
  FOS.F.LOW <- c(FOS.F.LOW, SigGenes[[i]])
}
genes <- FOS.F.LOW <- unique(FOS.F.LOW)
a2 <- data.frame(row.names = genes)
for (j in c("DG","CA1","VIP")){
  nm <- paste(j,"FOSplow",sep ="." )
  a2[SigGenes[[nm]],nm] <- 1
  a2[is.na(a2[,nm]),nm] <- 0
}
a2$AllCells <- rowSums(a2)

a2$Cell <- NA
for (i in 1:nrow(a2)){
  if(a2[i,"AllCells"] == 1){
    a2[i,"Cell"] <- do.call("rbind",strsplit(colnames(a2)[which(a2[i,c(1:3)] == 1)],split = ".",fixed = TRUE))[,1]
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
for (j in c("DG","CA1","VIP")){
  nm <- paste(j,"FOSphigh",sep ="." )
  a2[SigGenes[[nm]],nm] <- 1
  a2[is.na(a2[,nm]),nm] <- 0
}
a2$AllCells <- rowSums(a2)
a2$Cell <- NA
for (i in 1:nrow(a2)){
  if(a2[i,"AllCells"] == 1){
    a2[i,"Cell"] <- do.call("rbind",strsplit(colnames(a2)[which(a2[i,c(1:3)] == 1)],split = ".",fixed = TRUE))[,1]
  }else{
    "none"
  }
}


FOS.N.LOW <- vector()
for (i in seq(3,length(SigGenes),4)){
  FOS.N.LOW <- c(FOS.N.LOW, SigGenes[[i]])
}
FOS.N.LOW <- unique(FOS.N.LOW)
FOS.F.HIGH <- unique(FOS.F.HIGH)
genes <- FOS.F.HIGH <- unique(FOS.F.HIGH)
a2 <- data.frame(row.names = genes)
for (j in c("DG","CA1","VIP")){
  nm <- paste(j,"FOSphigh",sep ="." )
  a2[SigGenes[[nm]],nm] <- 1
  a2[is.na(a2[,nm]),nm] <- 0
}
a2$AllCells <- rowSums(a2)
a2$Cell <- NA
for (i in 1:nrow(a2)){
  if(a2[i,"AllCells"] == 1){
    a2[i,"Cell"] <- do.call("rbind",strsplit(colnames(a2)[which(a2[i,c(1:3)] == 1)],split = ".",fixed = TRUE))[,1]
  }else{
    "none"
  }
}


FOS.N.HIGH <- vector()
for (i in seq(4,length(SigGenes),4)){
  FOS.N.HIGH <- c(FOS.N.HIGH, SigGenes[[i]])
}
FOS.N.HIGH <- unique(FOS.N.HIGH)