#Scale expression
#Identify data
samples <- metaProxC[metaProxC$Mouse_condition == "HC" & metaProxC$outliers == "in" & metaProxC$Subgroup2 != "HDG",
                     "Sample_ID"]
dat <- na.exclude(tpmProxC[, samples])
met <- metaProxC[match(samples,metaProxC$Sample_ID),]
#genes <- c("Arc","Fos","Egr1","Fosb","Junb")
#genes <- rownames(na.exclude(dat))
#Scale data
dat.s <- na.exclude(t(scale(t(dat[genes,]))))
#For each sample, calculate the sum of scaled values
ieg <- data.frame(ieg = as.numeric(colSums(dat.s)),
                  celltype = met$Subgroup2,
                  FOS = met$FOS,
                  row.names = colnames(dat))

ieg[ieg$celltype == "CA1b","celltype"] <- "CA1"
ieg.sorted <- vector()
for (x in unique(ieg$celltype)){
  tmp <- ieg[ieg$celltype == x,]
  tmp <- tmp[order(tmp$ieg),]
  tmp$rank <- c(1:nrow(tmp))
  ieg.sorted <- rbind(ieg.sorted,tmp)
}

ieg.sorted <- ieg.sorted[ieg.sorted$FOS != "L",]
ggplot(ieg.sorted, aes(rank, ieg, colour = FOS,shape = celltype))+
  geom_point(size = 4)+
  xlab("Rank, cells")+
  ylab("All Genes Score")+
  scale_colour_manual(values=c("red","blue"))+
  facet_grid(~celltype,scales = "free_x")

####
ieg.sorted$celltype <- factor(ieg.sorted$celltype, levels = c("IN","VIP","CA1","CA3","DG"))
model <- glm(as.factor(FOS) ~ ieg, ieg.sorted[ieg.sorted$celltype == "DG",],family = "binomial")
summary(model)

ggplot(ieg.sorted, aes(FOS, ieg))+
  geom_violin()+
  geom_point()+
  facet_grid(~celltype)

####
