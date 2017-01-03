##################
# PREDICTING SUBTYPE FROM HIPPOSEQ DATA
##################
library(biomaRt)
mart <- useEnsembl(dataset = "mmusculus_gene_ensembl", biomart = "ensembl", version = 83)
library(randomForest)
### Get Hipposeq data
diftest <- read.table("~/Documents/SalkProjects/ME/ShortLongSingature/raw/GSE74985_gene_exp.diff",header=TRUE)
mergedcounts <- read.table("~/Documents/SalkProjects/ME/ShortLongSingature/raw/GSE74985_mergedCount.txt",header=TRUE,row.names = 1)
### Identify predictive genes
difgenes <- as.character(diftest[diftest$significant == "yes","gene_id"])
difcelltype <- ifelse(test = as.character(diftest[diftest$significant == "yes","value_1"] < diftest[diftest$significant == "yes","value_2"]), yes = as.character(diftest[diftest$significant == "yes","sample_2"]), no =  as.character(diftest[diftest$significant == "yes","sample_1"]))
difgroup <- unique(data.frame(difgenes, difcelltype))
nm <- table(difgenes)
difgroup <- difgroup[match(names(nm[nm==1]), difgroup$difgenes),]


### Convert Ensembl to external gene name
#genes
genes <- as.character(difgroup$difgenes)
convert <- getBM(
  filters = "ensembl_gene_id",
  attributes = c("ensembl_gene_id","external_gene_name"),
  values = c(genes),
  mart = mart,
  verbose = FALSE
)
passed <- (convert[match(genes,convert$ensembl_gene_id),"ensembl_gene_id"])
passed <- passed[!is.na(passed)]
difgroup <- difgroup[match(passed, difgroup$difgenes),]
rownames(difgroup) <- convert[match(passed,convert$ensembl_gene_id),"external_gene_name"]

######## Build classifier using these genes
# Note: rank order the genes because the analyses were performed differently

counts <- mergedcounts[difgroup$difgenes,]
rownames(counts) <- rownames(difgroup)
counts <- t(counts)
counts <- data.frame(apply(counts, 2, rank))
a <- do.call("rbind",strsplit(rownames(counts),"_"))[,1:2]
a <- as.vector(unlist(apply(X = a,MARGIN =   1, FUN = paste, collapse = "_")))
counts$celltype <- factor(a)

Rf <- randomForest(celltype ~. , counts)
gini <- data.frame(MeanDecreaseGini = Rf$importance, NA)
gini <- gini[order(gini$MeanDecreaseGini, decreasing = TRUE),]



                      