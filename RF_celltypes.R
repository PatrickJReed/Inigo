library(randomForest)

samples <- rownames(metaProxC[ metaProxC$Context1 == "none" & metaProxC$FOS != "L" & metaProxC$cluster_outlier == "in" & metaProxC$outliers == "in" & metaProxC$Arc_2.5 != "greater",])
dat <- na.exclude(tpmProxC[, samples])
dat <- dat[rowSums(dat) > 0,]
met <- metaProxC[samples,]

dat2 <- data.frame(t(dat[rowSums(dat) > 10,]))
dat2$celltype <- as.factor(as.character(met$Subgroup2))

k <- round(ncol(dat2)/1000)
i <- c(seq(1,(k * 1000),1000))
j <- c(seq(1000,((k-1) * 1000),1000), (ncol(dat2)-1))

RF.gini <- vector()
for(n in 1:length(i)){
  a <- randomForest(celltype ~., dat2[,c(i[n]:j[n],25880)])
  RF.gini <- rbind(RF.gini, a$importance)
}

RF.gini <- data.frame(RF.gini,NA)
RF.gini <- RF.gini[order(RF.gini$MeanDecreaseGini,decreasing = TRUE),]

genes <- rownames(RF.gini[log(RF.gini$MeanDecreaseGini) > 0,])
genes <- genes[-c(grep("Gm",genes))]

RF2 <- randomForest(celltype ~., cbind(dat2[,genes], celltype = dat2[,25880]))
RF.gini2 <- data.frame(RF2$importance,NA)
RF.gini2 <- RF.gini2[order(RF.gini2$MeanDecreaseGini,decreasing=TRUE),]

genes2 <- rownames(RF.gini2[RF.gini2$MeanDecreaseGini > 1,])

RF3 <- randomForest(celltype ~., cbind(dat2[,genes2], celltype = dat2[,25880]))
genes2 <- genes2[-c(grep("Rik", genes2))]

tmp <- na.exclude(dat[genes2,])

met$predicted <- as.vector(RF3$predicted)

RF4 <- randomForest(predicted ~., cbind(dat2[,genes2], predicted = as.vector(RF3$predicted)))

