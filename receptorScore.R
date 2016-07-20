samples <- metaProxC[metaProxC$Context1 != "A" & metaProxC$outliers == "in" & metaProxC$Subgroup2 != "HDG",
                     "Sample_ID"]
dat <- na.exclude(tpmProxC[, samples])
met <- metaProxC[match(samples,metaProxC$Sample_ID),]
g <- read.table(as.matrix("~/Documents/SalkProjects/ME/ShortLongSingature/SLSig_R/receptors.txt"))
g$V1 <- as.character(g$V1)
#genes <- rownames(na.exclude(dat))
#Scale data
ieg.sorted2 <- vector()
for(group in as.character(unique(g$V2))){
  genes <- as.character(g[g$V2 == group,"V1"])
  dat.s <- na.exclude(t(scale(t(dat[genes,]))))
  #For each sample, calculate the sum of scaled values
  ieg <- data.frame(ieg = as.numeric(colSums(dat.s)),
                    celltype = met$Subgroup2,
                    FOS = met$FOS,
                    condition = met$Mouse_condition,
                    row.names = colnames(dat))
  
  ieg[ieg$celltype == "CA1b","celltype"] <- "CA1"
  ieg.sorted <- vector()
  for (x in unique(ieg$celltype)){
    tmp <- ieg[ieg$celltype == x,]
    tmp <- tmp[order(tmp$ieg),]
    tmp$rank <- c(1:nrow(tmp))
    ieg.sorted <- rbind(ieg.sorted,tmp)
  }
  ieg.sorted[,"group"] <- group
  ieg.sorted$ieg <- ieg.sorted$ieg / length(genes)
  ieg.sorted2 <- rbind((ieg.sorted2), as.matrix(ieg.sorted))
}
ieg.sorted2 <- as.data.frame(ieg.sorted2)
ieg.sorted2$ieg <- as.numeric(as.character(ieg.sorted2$ieg))

ieg.sorted2 <- ieg.sorted2[ieg.sorted2$FOS != "L",]
ieg.sorted2$celltype  <- as.character(ieg.sorted2$celltype)
ggplot(ieg.sorted2[ieg.sorted2$group == "Ca" &ieg.sorted2$celltype != "CA2",], aes(celltype, ieg, colour = group))+
  geom_violin()+
  facet_grid(~condition )
