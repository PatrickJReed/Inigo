#########################
## Determining fdr procedure
## Conclusion: Use p.adjust(method = "fdr") with an alpha threshold of 0.05
#########################
library(fdrtool)
library(reshape)
library(ggplot2)
##
a <- res.fosN_fosP
b <- res.NE_P_C_FvN


nm <- table(c(rownames(a),rownames(b)))
nm <- names(nm[nm == 2])

a2 <- a[nm,]
b2 <- b[nm,]

tmp <- data.frame(a2,b2)
#First just using different alphas
cosig <- data.frame()
j <- 0
for (i in seq(1,40,0.5)){
  j <- j+1
  cosig[j,"thresh"] <- i
  cosig[j,"F_old"] <- sum(tmp$logFC > 0 & -log(tmp$PValue) > i)
  cosig[j,"F_new"] <- sum(tmp$logFC.1 < 0 & -log(tmp$PValue.1) > i)
  cosig[j,"N_old"] <- sum(tmp$logFC < 0 & -log(tmp$PValue) > i)
  cosig[j,"N_new"] <- sum(tmp$logFC.1 > 0 & -log(tmp$PValue.1) > i)
  cosig[j,"F_F"] <- sum(tmp$logFC > 0 & tmp$logFC.1 < 0 & -log(tmp$PValue) > i & -log(tmp$PValue.1) > i)
  cosig[j,"N_N"] <- sum(tmp$logFC < 0 & tmp$logFC.1 > 0 & -log(tmp$PValue) > i & -log(tmp$PValue.1) > i)
  cosig[j,"F_N"] <- sum(tmp$logFC > 0 & tmp$logFC.1 > 0 & -log(tmp$PValue) > i & -log(tmp$PValue.1) > i)
  cosig[j,"N_F"] <- sum(tmp$logFC < 0 & tmp$logFC.1 < 0 & -log(tmp$PValue) > i & -log(tmp$PValue.1) > i)
  
}

## Plot raw data
cosig2 <- melt(t(cosig[,-c(1)]))
cosig2$thresh <- rep(cosig$thresh, each = 4)
ggplot(cosig2, aes(thresh, value, colour = X1 ))+
  geom_point()

## Plot as percent of total
cosig3 <- cosig
cosig3$F_F_pA <- cosig3$F_F / cosig3$F_old
cosig3$F_N_pA <- cosig3$F_N / cosig3$F_old
cosig3$N_F_pA <- cosig3$N_F / cosig3$N_old
cosig3$N_N_pA <- cosig3$N_N / cosig3$N_old
cosig3$F_F_pB <- cosig3$F_F / cosig3$F_new
cosig3$F_N_pB <- cosig3$F_N / cosig3$F_new
cosig3$N_F_pB <- cosig3$N_F / cosig3$N_new
cosig3$N_N_pB <- cosig3$N_N / cosig3$N_new
cosig3 <- cosig3[,-c(1:9)]
cosig4 <- melt(t(cosig3))
cosig4$thresh <- rep(cosig$thresh, each = 8)

ggplot(cosig4, aes(thresh, value, colour = X1))+
  geom_point()+
  geom_smooth()

###########################
## Using fdrtoools
###########################
sigF <- data.frame()
alpha <- 0.01
for(alpha in c(0.05, 0.01)){
for (j in p.adjust.methods){
  f <- p.adjust(tmp$PValue,method = j)#fdrtool(x = tmp$PValue, statistic = "pvalue",cutoff.method = "pct0",pct0 = i,plot = FALSE)
  tmp$f_A <- f
  f <- p.adjust(tmp$PValue.1,method = j)#fdrtool(x = tmp$PValue.1, statistic = "pvalue",cutoff.method = "pct0",pct0 = i, plot= FALSE)
  tmp$f_B <- f
  j <- paste(j,alpha,sep = ".")
  sigF[j,"thresh"] <- alpha
  sigF[j,"F_old"] <- sum(tmp$logFC > 0 & (tmp$f_A)  <  alpha)
  sigF[j,"F_new"] <- sum(tmp$logFC.1 < 0 & (tmp$f_B) <  alpha)
  sigF[j,"N_old"] <- sum(tmp$logFC < 0 & (tmp$f_A) <  alpha)
  sigF[j,"N_new"] <- sum(tmp$logFC.1 > 0 & (tmp$f_B) <  alpha)
  sigF[j,"F_F"] <- sum(tmp$logFC > 0 & tmp$logFC.1 < 0 & (tmp$f_A) <  alpha & (tmp$f_B) <  alpha)
  sigF[j,"N_N"] <- sum(tmp$logFC < 0 & tmp$logFC.1 > 0 & (tmp$f_A) <  alpha & (tmp$f_B) <  alpha)
  sigF[j,"F_N"] <- sum(tmp$logFC > 0 & tmp$logFC.1 > 0 & (tmp$f_A) <  alpha & (tmp$f_B) <  alpha)
  sigF[j,"N_F"] <- sum(tmp$logFC < 0 & tmp$logFC.1 < 0 & (tmp$f_A) <  alpha & (tmp$f_B) <  alpha)
}
}

sigF$F_F_pA <- sigF$F_F / sigF$F_old
sigF$F_N_pA <- sigF$F_N / sigF$F_old
sigF$N_F_pA <- sigF$N_F / sigF$N_old
sigF$N_N_pA <- sigF$N_N / sigF$N_old
sigF$F_F_pB <- sigF$F_F / sigF$F_new
sigF$F_N_pB <- sigF$F_N / sigF$F_new
sigF$N_F_pB <- sigF$N_F / sigF$N_new
sigF$N_N_pB <- sigF$N_N / sigF$N_new

## Plot raw data
sigF2 <- melt(t(sigF[,c(10:17)]))
ggplot(sigF2, aes(X2, value, fill = X1 ))+
  geom_bar(stat = 'identity', position = "dodge")

### Find genes that are being lost with the BY method
j <- "BH"
f <- p.adjust(tmp$PValue,method = j)
tmp$f_A <- f
f <- p.adjust(tmp$PValue.1,method = j)
tmp$f_B <- f
a <- rownames(tmp[tmp$logFC.1 > 0  & round(tmp$f_B,2) <=  alpha | 
                    tmp$logFC.1 < 0 & round(tmp$f_B,2) <=  alpha 
                  ,])

j <- "BY"
f <- p.adjust(tmp$PValue,method = j)
tmp$f_A <- f
f <- p.adjust(tmp$PValue.1,method = j)
tmp$f_B <- f
b <- rownames(tmp[tmp$logFC.1 > 0  & round(tmp$f_B,2) <=  alpha | 
                    tmp$logFC.1 < 0 & round(tmp$f_B,2) <=  alpha 
                  ,])

Join <- table(c(a,a,b))
sum(Join == 1)
sum(Join == 2)
sum(Join == 3)


tmp2 <- res.NE_P_C_FvN[names(Join[Join == 2]),]
tmp2 <- tmp2[order(tmp2$PValue),]
