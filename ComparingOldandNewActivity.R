a <- res.fosN_fosP
b <- RES[[9]]#res.NE_P_C_FvN

nm <- table(c(rownames(a),rownames(b)))
nm <- names(nm[nm == 2])

a2 <- a[nm,]
b2 <- b[nm,]

nm2 <- rownames(a2[a2$f < 0.05,])


a3 <- a[nm2,]
b3 <- b[nm2,]

tmp <- data.frame(a2,b2)


ggplot(tmp[-c(which(rownames(tmp) == "Arc")),], aes(logFC, -logFC.1, colour = -log(PValue.1), alpha = -log(PValue)))+
  geom_point()+
  geom_smooth(method = "lm")+
  scale_color_gradient(high = "red",low = "blue")+
  xlab("logFC\nNat Comms")+
  ylab("logFC\nMolecular Dissection")+
  labs(title = "Comparison of Results\nDG after 15min NE")

sum(tmp$f < 0.05 & tmp$logFC < 0 & tmp$f.1 < 0.05 & tmp$logFC.1 > 0)
sum(tmp$f < 0.05 & tmp$logFC > 0 & tmp$f.1 < 0.05 & tmp$logFC.1 < 0)
sum(tmp$f < 0.05 & tmp$logFC < 0 & tmp$f.1 < 0.05 & tmp$logFC.1 < 0)
sum(tmp$f < 0.05 & tmp$logFC > 0 & tmp$f.1 < 0.05 & tmp$logFC.1 > 0)
