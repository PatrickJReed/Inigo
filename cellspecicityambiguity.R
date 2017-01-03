df <- vector()
for (cell in c("CA1","CA1b","DG","IN","VIP")){
  #others F  
  samples <- metaProxC[  metaProxC$Subgroup2 != cell & metaProxC$FOS == "F" & metaProxC$Mouse_condition == "EE" & metaProxC$Context1 == "none" &  metaProxC$FOS != "L" & metaProxC$Subgroup2 != "HDG" & metaProxC$outliers == "in" ,
                           "Sample_ID"]#
    dat <- na.exclude(tpmProxC[, samples])
    met <- metaProxC[match(samples,metaProxC$Sample_ID),]
    genes <- difexp_uniq[[paste(cell,"F",sep=".")]]
    logFC <- RES.F[[cell]][genes,"logFC"]
    m.f <- as.vector(apply(dat[genes,], 1, mean))
    p.f <- as.vector(apply(dat[genes,], 1, propExp))

    samples <- metaProxC[  metaProxC$Subgroup2 != cell & metaProxC$FOS == "N" & metaProxC$Mouse_condition == "EE" & metaProxC$Context1 == "none" &  metaProxC$FOS != "L" & metaProxC$Subgroup2 != "HDG" & metaProxC$outliers == "in" ,
                           "Sample_ID"]#
    dat <- na.exclude(tpmProxC[, samples])
    met <- metaProxC[match(samples,metaProxC$Sample_ID),]
    logFC <- RES.F[[cell]][genes,"logFC"]
    m.n <- as.vector(apply(dat[genes,], 1, mean))
    p.n <- as.vector(apply(dat[genes,], 1, propExp))
    
    
    
      samples <- metaProxC[  metaProxC$Subgroup2 == cell & metaProxC$FOS == "F" & metaProxC$Mouse_condition == "EE" & metaProxC$Context1 == "none" &  metaProxC$FOS != "L" & metaProxC$Subgroup2 != "HDG" & metaProxC$outliers == "in" ,
                           "Sample_ID"]#
    dat <- na.exclude(tpmProxC[, samples])
    met <- metaProxC[match(samples,metaProxC$Sample_ID),]
    m1.f <- as.vector(apply(dat[genes,], 1, mean))
    
    samples <- metaProxC[  metaProxC$Subgroup2 == cell & metaProxC$FOS == "N" & metaProxC$Mouse_condition == "EE" & metaProxC$Context1 == "none" &  metaProxC$FOS != "L" & metaProxC$Subgroup2 != "HDG" & metaProxC$outliers == "in" ,
                           "Sample_ID"]#
    dat <- na.exclude(tpmProxC[, samples])
    met <- metaProxC[match(samples,metaProxC$Sample_ID),]
    m1.n <- as.vector(apply(dat[genes,], 1, mean))
    
    
  df <- rbind(df, cbind(genes,logFC,p.f,p.n,m.f,m.n,m1.f,m1.n,cell))
}
df <- as.data.frame(df)
df$logFC_FvsN <- -1 * as.numeric(as.character(df$logFC))
df$m1.f <- as.numeric(as.character(df$m1.f))
df$p.f <- as.numeric(as.character(df$p.f))
df$p.n <- as.numeric(as.character(df$p.n))
df$m.f <- as.numeric(as.character(df$m.f))
df$m1.n <- as.numeric(as.character(df$m1.n))
df$m.n <- as.numeric(as.character(df$m.n))
df$Remaining_FOS_N <- df$m.n
df$Remaining_FOS_F <- df$m.f

#df$Mean_Interest_FOSN


ggplot(df[df$cell != "CA1b",], aes(m1.f, m.n, colour = Remaining_FOS_F, shape = cell))+
  geom_point()+
  scale_colour_gradient(high = "red", low = "blue")+
  xlab("Mean Expression in Cell-type of Interest\nFOS+")+
  ylab("Mean Expression in Remaining Cells\nFOS-")+
  labs(title=("Cell-Specific Activity-Responsive Genes"))+
  geom_abline(slope = 1,intercept = 0)+
  theme_bw()+
  theme(text=element_text(size=16))+
  theme(panel.border = element_rect(colour=c("black"),size=2),
        axis.ticks = element_line(size=1.5),
        panel.grid.major = element_line(colour = "black"))
