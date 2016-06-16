##################
# Compare Normalization
# - Identify which normalization method reduces variation and maintains 
#     cell specificity between known DG and CA1 genes from high, mid, and low 
#     ranges in expression
#######
# # METHODS TO TEST
# 1. Pairwise Normalization (Regev)
# 2. Scale Factors (Anders, Deseq)
# 3. ? (EdgeR) #note: not for dimensionality reduction but to ensure that difexp is ok
# 4. TPM (RSEM)
# 5. ERCC (Ding et al)
##################
prodNoZero <- function(x){
  x <- x[!is.na(x)]
  return(prod(x[x>0]))
}
MedianScale <- function(i){
  genes <- which(tpm[,i] > 0)
  medianScale <- median(tpm[genes,i] / as.numeric(geoMeanNoZero)[genes])
  return(medianScale)
}
median.scale <- function(tpm){
  geoMeanNoZero <- apply(tpm, 1, prodNoZero) ^ (1/ncol(tpm))
  scaleFactor <- unlist(lapply(1:ncol(tpm), MedianScale))
  tpm2 <- matrix(0,nrow(tpm),ncol(tpm))
  for(i in 1:ncol(tpm)){
    tpm2[,i] <- tpm[,i] / scaleFactor[i]
  }
  rownames(tpm2) <- rownames(tpm)
  colnames(tpm2) <- colnames(tpm)
  return(tpm2)
  
}
deseq.scale <- function(count){
  require("DESeq2")
    colData <- data.frame(a = rep("a",ncol(count)),c(rep("paired-end",ncol(count))))
    colnames(colData) <- c("a","type")
    dds <- DESeqDataSetFromMatrix(countData = count,
                                  colData = colData,
                                  design = ~ 1)
    norm <- DESeq(dds)
    count2 <- matrix(0,nrow = nrow(count),ncol = ncol(count))
    for(i in 1:ncol(count)){
      count2[,i] <- count[,i] / sizeFactors(norm)[i]
    }
    rownames(count2) <- rownames(count)
    colnames(count2) <- colnames(count)
    return(count2)
}
edger.scale <- function(count){
  require(edgeR)
  cds <- DGEList(count)
  cds <- cds[rowSums(1e+06 * cds$counts/expandAsMatrix(cds$samples$lib.size, dim(cds)) > 1) >= 3, ]
  cds <- calcNormFactors( cds )
  cds <- estimateCommonDisp( cds )
  cds <- estimateTagwiseDisp( cds )
  cds <- estimateTrendedDisp(cds)
  return(cds$pseudo.counts)
}
ercc.scale <- function(count, ercc){

}

##################
# Step 1. Identify Samples
##################
samples <- metaProxC[ metaProxC$Mouse_condition == "HC"  & metaProxC$alignable >  100000 & metaProxC$Smartseq2_RT_enzyme_used == "ProtoscriptII", 
                      "Sample_ID"]
met <- metaProxC[match(samples,metaProxC$Sample_ID),]
count <- na.exclude(countProxC[, samples])
tpm <- na.exclude(tpmProxC[, samples])
##################
# Step 2. Normalize Data
##################

d <- log(deseq.scale(count)+1,2)
e <- log(edger.scale(count) + 1,2)
t <- tpm
m <- median.scale(tpm)
c <- log(count+1,2)

##################
# Step 3. Test
##################

genes <- c("Prox1","Meg3","Snhg3","Rfx3")
aic <- list(0,0,0,0,0)
names(aic) <- c("d","e","t","m","c")
for (gene in genes){
  tmp <- data.frame(d = as.numeric(d[gene,]), 
                    e  = as.numeric(e[gene,]),
                    t = as.numeric(t[gene,]),
                    m = as.numeric(m[gene,]),
                    c = as.numeric(c[gene,]),
                    group = met$Brain_Region)


  model <- glm(group == "DG" ~ d , tmp, family = "binomial")
  aic[["d"]] <- aic[["d"]] + as.numeric(summary(model)$aic)
  model <- glm(group == "DG" ~ e , tmp, family = "binomial")
  aic[["e"]] <- aic[["e"]] + as.numeric(summary(model)$aic)
  model <- glm(group == "DG" ~ t , tmp, family = "binomial")
  aic[["t"]] <- aic[["t"]] + as.numeric(summary(model)$aic)
  model <- glm(group == "DG" ~ m , tmp, family = "binomial")
  aic[["m"]] <- aic[["m"]] + as.numeric(summary(model)$aic)
  model <- glm(group == "DG" ~ c , tmp, family = "binomial")
  aic[["c"]] <- aic[["c"]] + as.numeric(summary(model)$aic)
}
names(aic)[which(unlist(aic) == min(unlist(aic)))]
