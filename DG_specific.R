####
#save(list = c("RES.F", "RES.Direction", "df_uniq", "difexp_uniq", "difexp","cell.Fspecific", "cell.quiet"),file = "~/Documents/SalkProjects/ME/ShortLongSingature/SLSig_R/cellspecific.rda",compress = TRUE)
#load("~/Documents/SalkProjects/ME/ShortLongSingature/SLSig_R/cellspecific.rda")
# Distance between FOS+ and FOS-
Test <- function(celltype = "DG"){
  samples <- rownames(metaProxC[metaProxC$Subgroup2 == celltype &
                                  metaProxC$Mouse_condition == "EE" & metaProxC$FOS != "L" &metaProxC$outliers == "in",])
  dat <- na.exclude(countProxC[, samples])
  dat <- dat[rowSums(dat) > 0,]
  met <- metaProxC[match(samples,metaProxC$Sample_ID),]
  res <- exact(dat, group = (met$FOS == "F"), Pair = c(TRUE,FALSE))
  return(res)
}
############################
# Step 1: Calculate Genes Differentially expressed betwen FOS+ and FOS- depending on cell types
############################
res.DG <- Test(celltype = "DG")
res.CA1 <- Test(celltype = "CA1")
res.CA1b <- Test(celltype = "CA1b")
res.CA3 <- Test(celltype = "CA3")
res.IN <- Test(celltype = "IN")
res.VIP <- Test(celltype = "VIP")

RES.F <- list(res.DG, res.CA1, res.CA1b, res.CA3, res.IN, res.VIP)
names(RES.F) <- c("DG","CA1","CA1b","CA3","IN","VIP")

############################
# Step 2: Plot a count of all differentially expressed genes by cell type
############################
p.a <- 0.3
difexp <- list()
df <- vector()
i <- 0
for (nm in names(RES.F)){
  for (j in c("F","N")){
    i <- i + 1
    r <- RES.F[[nm]]
    if (j == "F"){
      difexp[[i]] <- rownames(r[r$logFC < 0 & r$f < 0.05 & r$a > p.a,])
    }else{
      difexp[[i]] <- rownames(r[r$logFC > 0 & r$f < 0.05 & r$a < p.a ,])
    }
    names(difexp)[i] <- paste(nm, j, sep =".")
    df <- rbind(df, c(nm, j, length(difexp[[i]])))
  }
}
df <- as.data.frame(df)
df$V3 <- as.numeric(as.character(df$V3))
df$V2 <- factor(df$V2, levels = c("N","F"))

ggplot(df[df$V1 != "CA1b" & df$V1 != "CA3",], aes(V1, V3, fill = V2))+
  geom_bar(stat = 'identity', position = "dodge")+
  scale_fill_manual(values = c("blue","red"))+
  ylab("DifExp Genes, count")+
  labs(title = "Differential Expression")+
  theme_bw()+
  theme(text=element_text(size=20))+
  theme(panel.border = element_rect(colour=c("black"),size=2),
        axis.ticks = element_line(size=1.5))

############################
# Step 3: Genes that are only different in one condition
############################
df_uniq <- vector()
difexp_uniq <- list()
i <- 0
for (nm in names(difexp)){
  i <- i + 1
  other <- names(difexp)[-which(names(difexp) == nm)]
  a <- rep(difexp[[nm]],(length(difexp)+1))
  for (j in 1:length(other)){
    a <- c(a,difexp[[other[j]]])
  }
  a <- table(a)
  difexp_uniq[[i]] <- names(a[a== (length(difexp)+1)])
  names(difexp_uniq)[i] <- nm
  uniq <- sum(a == (length(difexp)+1))
  nm2 <- unlist(strsplit(nm, ".", fixed = TRUE))
  df_uniq <- rbind(df_uniq, c(nm2[1],nm2[2],uniq))
}

df_uniq <- as.data.frame(df_uniq)
df_uniq$V3 <- as.numeric(as.character(df_uniq$V3))
df_uniq$V2 <- factor(df_uniq$V2, levels = c("N","F"))

ggplot(df_uniq[df_uniq$V1 != "CA1b" & df_uniq$V1 != "CA3",], aes(V1, V3, fill = V2))+
  geom_bar(stat = 'identity', position = "dodge")+
  scale_fill_manual(values = c("blue","red"))+
  ylab("Unique DifExp Genes, count")+
  labs(title = "Differential Expression")+
  theme_bw()+
  theme(text=element_text(size=20))+
  theme(panel.border = element_rect(colour=c("black"),size=2),
        axis.ticks = element_line(size=1.5))

######
# Step 4: Generate Differential Expression matrices for
#         Subtype|FOS vs all else
#####
subs <- c("CA1","CA1b","CA3","DG","IN","VIP")

RES.Direction <- list()
i <- 0
for(s in subs){
  for(j in c("F","N")){
    i <- i + 1
    nm <- paste(s,j,sep = ".")
    if(j == "F"){
    samples <- rownames(metaProxC[metaProxC$Mouse_condition == "EE" &  metaProxC$Subgroup2 != s & metaProxC$FOS == "N" & metaProxC$Subgroup2 != "HDG" & metaProxC$outliers == "in"|
                                    metaProxC$Mouse_condition == "EE" &  metaProxC$Subgroup2 == s & metaProxC$FOS == "F" & metaProxC$outliers == "in",])
    }else{
      samples <- rownames(metaProxC[metaProxC$Mouse_condition == "EE" &  metaProxC$Subgroup2 != s & metaProxC$FOS == "F" & metaProxC$Subgroup2 != "HDG" & metaProxC$outliers == "in"|
                                      metaProxC$Mouse_condition == "EE" &  metaProxC$Subgroup2 == s & metaProxC$FOS == "N" & metaProxC$outliers == "in",])
    }
    dat <- na.exclude(countProxC[, samples])
    dat <- dat[rowSums(dat) > 0,]
    met <- metaProxC[match(samples,metaProxC$Sample_ID),]
    res <- exact(dat, group = (met$FOS == j), Pair = c(TRUE,FALSE))
    ##Save difexp results 
    RES.Direction[[i]] <- res
    names(RES.Direction)[i] <- nm 
  }
}

######
# Step 5: Save genes that are either cell-specific, or cell-quiet
######
subs <- c("CA1","CA1b","DG","IN","VIP")
cell.Fspecific <- list()
cell.quiet <- list()
i2 <- 0
for(s in unique(subs)){
  for(j in c("F","N")){
    nm <- paste(s,j,sep = ".")
    res <- RES.Direction[[nm]]
    if(j == "F"){
      i2 <- i2 + 1
      genes <- difexp_uniq[[nm]]
      res2 <- res[genes,]
      a <- res2[res2$logFC < 0 & res2$f < 0.05,]
      a <- a[order(a$PValue),]
      cell.Fspecific[[i2]] <- a
      names(cell.Fspecific)[i2] <- s
    }else{
      genes <- difexp_uniq[[nm]]
      res2 <- res[genes,]
      cell.quiet[[i2]] <- rownames(na.exclude(res2[res2$logFC > 0 & res2$f < 0.05,]))
      names(cell.quiet)[i2] <- s
    }
  }
}

######
# Step 5: Count genes that are cell-specific
######
a <- RES.F[["VIP"]][cell.Fspecific[["VIP"]],]
a <- a[order(a$PValue),]
######
# Step 6: Count genes that are cell-independent
######
excitF <- table(c(difexp[[1]], difexp[[3]]))
excitN <- table(c(difexp[[2]], difexp[[4]]))
inhibF <- table(c(difexp[[9]], difexp[[11]]))
inhibN <- table(c(difexp[[10]], difexp[[12]]))
allF <- table(c(difexp[[1]], difexp[[3]],difexp[[9]], difexp[[11]]))
allN <- table(c(difexp[[2]], difexp[[4]],difexp[[10]], difexp[[12]]))
notinhib <- table(c(rep(names(excitF[excitF == 2]),3),difexp[[9]], difexp[[11]]))

a <- head(RES.F[[1]][names(notinhib[notinhib == 3]),])
a <- a[order(a$PValue),]

sum(excitF == 2)
sum(excitN == 2)
sum(inhibF == 2)
sum(inhibN == 2)
sum(allF == 4)
sum(allN == 4)
sum(notinhib == 3)

df <- res.DGonly[genes,]
df[dg.iegspecific,"group"] <- "DG Specific"
df[dg.quiet,"group"] <- "DG"

ggplot(df, aes(group, -logFC,colour = -log(f)))+
  geom_violin()+
  geom_point()+
  scale_colour_gradient(high = "red",low = 'blue')+
  ylab("DG (F+) vs. nonDG (F-)\nlogFC")+
  theme_bw()+
  theme(text=element_text(size=16))+
  theme(panel.border = element_rect(colour=c("black"),size=2),
        axis.ticks = element_line(size=1.5),
        panel.grid.major = element_line(colour = "black"))+
  labs(title=("Genes Activated in DG FOS+ Nuclei"))


res3 <- res.DG[dg.iegspecific,]
res3 <- res3[order(res3$PValue),]
###### 
# Identify genes that are already on in other cell types
###### 
samples <- rownames(metaProxC[metaProxC$Mouse_condition == "EE" &  metaProxC$Subgroup2 != "DG" & metaProxC$FOS == "N" & metaProxC$Subgroup2 != "HDG" & metaProxC$outliers == "in"|
 metaProxC$Mouse_condition == "EE" &  metaProxC$Subgroup2 == "DG" & metaProxC$FOS == "F" & metaProxC$outliers == "in",])
dat <- na.exclude(tpmProxC[, samples])
dat <- dat[rowSums(dat) > 0,]
met <- metaProxC[match(samples,metaProxC$Sample_ID),]


prop <- vector()
for (i in unique(met$Subgroup2)){
  samples2 <- rownames(met[met$Subgroup2 == i,])
  tmp2 <- apply(X = dat[,samples2], MARGIN = 1, FUN = propExp)
  prop <- cbind(prop, as.numeric(tmp2))
}
colnames(prop) <- unique(met$Subgroup2)
rownames(prop) <- rownames(dat)

#
prop2 <- prop[genes,]
prop3 <- melt(t(prop2[,-c(3)]))
prop3$DG <- rep(prop2[,3],each = 5)

ggplot(prop3[prop3$X1 !="CA3",], aes(DG, value, colour = X1))+
  geom_point()

##
dat2 <- as.matrix(dat[genes,])
colnames(dat2) <- paste(met$Subgroup2, met$FOS)
heatmap(dat2)
