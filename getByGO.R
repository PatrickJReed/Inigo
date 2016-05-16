library(reshape)
library(ggplot2)
library(biomaRt)
library(GO.db)

genes <- read.table(as.matrix("~/Documents/libraries/HOXgenes_panther_human.txt"))
genes <- as.character(genes$V1)

genes <- b#c(rownames(head(load,n=2)),rownames(tail(load,n=2)))
ensembl <- useMart(host="www.ensembl.org","ENSEMBL_MART_ENSEMBL",dataset="mmusculus_gene_ensembl")
#ensembl <- useMart(host="www.ensembl.org","ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl")
#filter <- c("external_gene_name") #TRY: listFilters(ensembl)#hgnc_symbol
#attrib <- c("external_gene_name","name_1006","go_id","namespace_1003")
filter = "external_gene_name"
attrib = c("external_gene_name","ensembl_gene_id")#"name_1006","go_id","namespace_1003")
res = getBM(attributes=attrib,filters=filter,values=genes,mart=ensembl)

resT <- table(res$name_1006)
resT <- resT[order(resT,decreasing = TRUE)]

resT[grep("homeobox",names(resT))]

tail(head(resT,n=80),n=20)
goterm <- "receptor binding"

b <- unique(res[grep(goterm,res$name_1006,ignore.case =TRUE),"external_gene_name"])

b2 <- RES[[8]][b,]
b2 <- b2[order(b2$f),]
