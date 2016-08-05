a <- cbind(res2, res[rownames(res2),])
a <- na.exclude(a)
#1
genes1 <- rownames(a[(a[,1] < 0  & 
       a[,4] < 0.01 &
       a[,7] < 0 & 
       a[,10] < 0.01
      ),])
#2
genes2 <- rownames(a[(a[,1] > 0  & 
      a[,4] < 0.01 &
      a[,7] > 0 & 
      a[,10] < 0.01
),])


a <- cbind(res2, res[rownames(res2),])
a <- na.exclude(a)


#3
nrow(a[(a[,1] > 0  & 
      a[,4] < 0.01 &
      a[,7] < 0 & 
      a[,10] < 0.01
),])
#4
nrow(a[(a[,1] < 0  & 
      a[,4] < 0.01 &
      a[,7] > 0 & 
      a[,10] < 0.01
),])

fosp <- rownames(res2[res2$logFC > 0 & res2$f < 0.01,])
tmp <- table(c(rep(fosp,4), rep(genes2,2), genes3))
genes <- names(tmp[tmp == 4])

tmp2 <- res2[genes,]
tmp2 <- tmp2[order(tmp2$f),]
