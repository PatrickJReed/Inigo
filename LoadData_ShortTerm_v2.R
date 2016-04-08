#################
## Common libraries and functions
#################

#names
Names <- function(x){
  if(length(grep(".",colnames(x)[1])) > 0){
    #colnames(x) <- as.vector(do.call("rbind",strsplit(colnames(x),fixed=TRUE,split="."))[,1])
    dotsplit <- strsplit(colnames(x),fixed=TRUE,split=".")
    #colnames(x) <- rbind.fill(lapply(x,function(y){as.data.frame(t(dotsplit),stringsAsFactors=FALSE)}))
    tmp <- vector()
    for (i in 1:length(dotsplit)){
      tmp <- c(tmp, dotsplit[[i]][1])
    }
    colnames(x) <- tmp
  }else{
    scoresplit <- strsplit(colnames(x),fixed=TRUE,split="_")
    tmp <- vector()
    for (i in 1:length(scoresplit)){
      tmp <- c(tmp, scoresplit[[i]][1])
    }
  }
    colnames(x) <- gsub(pattern="X",replacement="Sample_",x=colnames(x))
    a <- do.call("rbind",strsplit(colnames(x),split="Sample_"))[,2]
    l <- nchar(colnames(x)) - nchar(gsub(pattern="_",replacement="",x=colnames(x)))
    nam <- vector()
    for (i in 1:length(a)){
      if (l[i] == 7){
        nam <- c(nam, a[i])
      }else{
        if(add1 == TRUE){
        nam <- c(nam,paste(a[i],"_1",sep=""))
        }else{
          nam <- c(nam,paste(a[i],sep=""))
        }
      }
    }
    return(nam)
}
Meta <- function(x){
  metaProx <- as.data.frame(do.call("rbind",(strsplit(x=colnames(x),split="_"))))
  colnames(metaProx) <- c("date","well","prox","ctip","fos","cond","rep")
  rownames(metaProx) <- colnames(x)
  return(metaProx)
}

### Add QC to Meta (need to fix the code here)

propExp <- function(x,i=0){
  if(sum(is.na(x)) > i){
    sum(!is.na(x))/length(x)
  }else{
    sum(x > i)/length(x)
  }
}
rawExp <- function(x,i=2){
  if(sum(is.na(x)) > i){
    sum(!is.na(x))
  }else{
    sum(x > i)
  }
}
meanNoZero <- function(x,i=0){
  a <- mean(x[x>i])
  if(sum(is.na(x)) > i){
    a <- mean(x[!is.na(x)])
  }else if(sum(x>i) == 0){
    a <- 0
  }
  return(a)
}
varNoZero <- function(x){
  sd(x[x >0])
}
MedNoZero <- function(x){
  median(x[x >0])
}
countExp <- function(x){
  sum(x > 0)
}
scaleME <- function(x){
  return(scale(x)[,1])
}
#################
## Processed rdas
#################

#################
## Raw Data
#################
# ####alignment
# QCreadcount <- read.table("~/Documents/SalkProjects/ME/WCellActivation_Manuscript/EE_QC_alignment_stats.txt",header=TRUE)
# QCreadcount <- QCreadcount[QCreadcount$group == "mouse",]
# rownames(QCreadcount) <- do.call("rbind",strsplit(as.character(QCreadcount$sample),split="-",fixed=TRUE))[,1]
# keepQC <- rownames(QCreadcount[QCreadcount$aligned > 50000,])
# 
# ## Exclude a known experimental outlier, < 50000 aligned reads
# Expoutlier <- c("nc_ux_ti_A12_141204","nm_ui_ti_G10_141204","nm_ux_ti_F7_141204")
# #tmp <- read.table(as.matrix("~/Documents/SalkProjects/BenLacar/PC1_excludes.txt"))
# #Expoutlier <- unique(c(Expoutlier,as.character(tmp$V1)))
# ## Count
# #countQC1 <- read.table(as.matrix("~/Documents/gene_count.txt"),header=TRUE,row.names=1,fill=TRUE)
# countQC1 <- read.table(as.matrix("~/Desktop/RECOVERED_RSEM_geneSymbol_tpm_141204_allsamples.txt"),header=TRUE,row.names=1)
# countQC1 <- countQC1[,keepQC]
# countQC2 <- countQC1[,-(match(Expoutlier,colnames(countQC1)))]
# countQC <- round(countQC2)
# #countQC <- countQC[,na.exclude(match(keptnames,colnames(countQC)))]
# ## TPM
# #tpmQC1 <- read.table(as.matrix("~/Documents/SalkProjects/BenLacar/QC/QCtpm_genesymbol.txt"),header=TRUE,row.names=1)
# tpmQC1 <- read.table(as.matrix("~/Documents/SalkProjects/BenLacar/RSEM_geneSymbol_tpm_141204_allsamples.txt"),header=TRUE,row.names=1)
# tpmQC <- tpmQC1[,keepQC]
# tpmQC <- tpmQC[,-(match(Expoutlier,colnames(tpmQC)))]
# tpmQC <- log(tpmQC+1,2)
# #tpmQC <- tpmQC[,na.exclude(match(keptnames,colnames(tpmQC)))]
# ## labels
# labelsQC1 <- as.data.frame(do.call(rbind,strsplit(x=colnames(tpmQC1),split='_',fixed=TRUE)))
# rownames(labelsQC1) <- do.call("rbind",strsplit(x=colnames(tpmQC1),split=".",fixed=TRUE))[,1]
# labelsQC <- labelsQC1[keepQC,]
# labelsQC <- labelsQC[-c(match(Expoutlier,rownames(labelsQC))),]

#### ERCCs

####### 07/2015 Homecage Prox1+ Fos Low/Neg, sorted by Baptiste and Jerika, processed by Jerika
#counts
excludes <- c("150629_H4_P_C_N_HC_1") # 150629_H4_P_C_N_HC_1 has almost no genes expressed
#count.1 <- read.table(as.matrix("~/Documents/SalkProjects/ME/ShortLongSingature/raw/Prox1_48nuc_0715_geneSymb_count.txt"),header=TRUE,row.names=1)
count.1 <- read.table(as.matrix("~/Documents/SalkProjects/ME/ShortLongSingature/raw/Prox1_genenamecount.txt"),header=TRUE,row.names=1)
#countProx <- log(count.1[rowSums(count.1) > 0,]+1,2)
countProx <- round(count.1)
add1 <- FALSE
colnames(countProx) <- Names(countProx)
countProx <- countProx[,-c(which(colnames(countProx) == excludes))]


#tpm
tpm.1 <- read.table(as.matrix("~/Documents/SalkProjects/ME/ShortLongSingature/raw/Prox1_genenametpm.txt"),header=TRUE,row.names=1)
tpmProx <- log(tpm.1[rowSums(tpm.1) > 0,]+1,2)
add1 <- FALSE
colnames(tpmProx) <- Names(tpmProx)
tpmProx <- tpmProx[,-c(which(colnames(tpmProx) == excludes))]

#151207 CA1 Samples
#counts
count.1 <- read.table(as.matrix("~/Documents/SalkProjects/ME/ShortLongSingature/raw/Prox1_151207_gene_count_genename.txt"),header=TRUE,row.names=1)
colnames(count.1)[which(colnames(count.1) == "X151207_A8_N_C_N_EE_2.708.6.CAGAGAG.AGAGTAG_CAGAGAG.AGAGTAG_L006_R1_001_se")] <- "X151207_D8_N_C_N_EE_2.708.6.CAGAGAG.AGAGTAG_CAGAGAG.AGAGTAG_L006_R1_001_se"
countProx1512 <- round(count.1)
add1 <- TRUE
colnames(countProx1512) <- make.names(Names(countProx1512))
#tpm
tpm.1 <- read.table(as.matrix("~/Documents/SalkProjects/ME/ShortLongSingature/raw/Prox1_151207_gene_tpm_genename.txt"),header=TRUE,row.names=1)
colnames(tpm.1)[which(colnames(tpm.1) == "X151207_A8_N_C_N_EE_2.708.6.CAGAGAG.AGAGTAG_CAGAGAG.AGAGTAG_L006_R1_001_se")] <- "X151207_D8_N_C_N_EE_2.708.6.CAGAGAG.AGAGTAG_CAGAGAG.AGAGTAG_L006_R1_001_se"
tpmProx1512 <- log(tpm.1[rowSums(tpm.1) > 0,]+1,2)
add1 <- TRUE
colnames(tpmProx1512) <- Names(tpmProx1512)

#151214 
#counts
count.1 <- read.table(as.matrix("~/Documents/SalkProjects/ME/ShortLongSingature/raw/Prox1_151214_gene_count_genename.txt"),header=TRUE,row.names=1)
countProx151214 <- round(count.1)
b <- do.call("rbind",strsplit(colnames(count.1)[87:134],"_",fixed=TRUE))
colnames(countProx151214)[87:134] <- paste(b[,1],b[,2],sep="_")
a <- do.call("rbind",strsplit(colnames(countProx151214)[1:86],"_",fixed=TRUE))
colnames(countProx151214)[1:86] <-paste(a[,1],a[,2],a[,3],a[,4],a[,5],a[,6],a[,7],sep="_")

#tpm
tpm.1 <- read.table(as.matrix("~/Documents/SalkProjects/ME/ShortLongSingature/raw/Prox1_151214_gene_tpm_genename.txt"),header=TRUE,row.names=1)
tpmProx151214 <- log(tpm.1[rowSums(tpm.1) > 0,]+1,2)
colnames(tpmProx151214) <- colnames(countProx151214)#make.names(Names(tpmProx151214))


#160107 
#counts
count.1 <- read.table(as.matrix("~/Documents/SalkProjects/ME/ShortLongSingature/raw/Prox1_160107_gene_count.txt"),header=TRUE,row.names=1)
b <- do.call("rbind",strsplit(colnames(count.1),"_",fixed=TRUE))
colnames(count.1) <- paste(b[,1],b[,2],sep="_")
countProx160107 <- round(count.1)

#tpm
tpm.1 <- read.table(as.matrix("~/Documents/SalkProjects/ME/ShortLongSingature/raw/Prox1_160107_gene_tpm.txt"),header=TRUE,row.names=1)
tpmProx160107 <- log(tpm.1[rowSums(tpm.1) > 0,]+1,2)
colnames(tpmProx160107) <- colnames(countProx160107)

#160324 
#counts
count.1 <- read.table(as.matrix("~/Documents/SalkProjects/ME/ShortLongSingature/raw/Prox1_160324_gene_count.txt"),header=TRUE,row.names=1)
b <- do.call("rbind",strsplit(colnames(count.1),"_",fixed=TRUE))
colnames(count.1) <- paste(b[,1],b[,2],sep="_")
countProx160324 <- round(count.1)

#tpm
tpm.1 <- read.table(as.matrix("~/Documents/SalkProjects/ME/ShortLongSingature/raw/Prox1_160324_gene_tpm.txt"),header=TRUE,row.names=1)
tpmProx160324 <- log(tpm.1[rowSums(tpm.1) > 0,]+1,2)
colnames(tpmProx160324) <- colnames(countProx160324)
################################

#combine
runs <- 5
a <- table(c(rownames(tpmProx),rownames(tpmProx1512), rownames(tpmProx151214),rownames(tpmProx160107),rownames(tpmProx160324)))
tpmProxC <- cbind(tpmProx[names(a[a==runs]),],tpmProx1512[names(a[a==runs]),], tpmProx151214[names(a[a==runs]),],tpmProx160107[names(a[a==runs]),],tpmProx160324[names(a[a==runs]),])
colnames(tpmProxC) <- make.names(colnames(tpmProxC))
a <- table(c(rownames(countProx),rownames(countProx1512), rownames(countProx151214),rownames(countProx160107), rownames(countProx160324)))
countProxC <- cbind(countProx[names(a[a==runs]),],countProx1512[names(a[a==runs]),],countProx151214[names(a[a==runs]),],countProx160107[names(a[a==runs]),], countProx160324[names(a[a==runs]),])
colnames(countProxC) <- make.names(colnames(countProxC))

metaProxC <- read.table("~/Documents/SalkProjects/ME/ShortLongSingature/raw/snRNAseqSampleIDFile.txt",header=TRUE)
metaProxC$Sample_ID <- make.names(metaProxC$Sample_ID)
rownames(metaProxC) <- metaProxC$Sample_ID
