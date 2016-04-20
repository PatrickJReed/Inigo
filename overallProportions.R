## All samples
p <- apply(X = tpmProxC, MARGIN = 1, FUN = propExp)
m <- apply(X = tpmProxC, MARGIN = 1, FUN = meanNoZero)
tmp <- data.frame(detection = p, mean = m)

#Each subgroup HC
for (y in levels(as.factor(metaProxC$Mouse_condition))){
  for (x in levels(as.factor(metaProxC$Brain_Region))){
    samples <- metaProxC[metaProxC$Brain_Region == x &  metaProxC$FOS != "L" & metaProxC$Mouse_condition == y & metaProxC$alignable >  500000 & metaProxC$Smartseq2_RT_enzyme_used == "ProtoscriptII" ,"Sample_ID"]#|
    dat <- tpmProxC[, samples]
    nm1 <- paste(x,y,"p", sep = ".")
    tmp[,nm1] <- apply(X = dat, MARGIN = 1, FUN = propExp)
    nm2 <- paste(x,y,"m", sep = ".")
    tmp[,nm2] <- apply(X = dat, MARGIN = 1, FUN = meanNoZero)
    nm3 <- paste(x,y,"100", sep = ".")
    tmp[,nm3] <- tmp[,nm1] > 0.9
  }
}

ggplot(tmp, aes(mean,detection))+
  geom_point(alpha =0.05)+
  theme_bw()+
  labs(title = "Gene Detection")
