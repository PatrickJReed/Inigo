library(reshape)
library(ggplot2)
library(biomaRt)
library(GO.db)

genes <- a
ensembl <- useMart(host="www.ensembl.org","ENSEMBL_MART_ENSEMBL",dataset="mmusculus_gene_ensembl")
filter <- c("external_gene_name") #TRY: listFilters(ensembl)
attrib <- c("external_gene_name","name_1006","go_id","namespace_1003")
res = getBM(attributes=attrib,filters=filter,values=genes,mart=ensembl)

resT <- table(res$name_1006)
resT <- resT[order(resT,decreasing = TRUE)]

tail(head(resT,n=20),n=20)
goterm <- "positive regulation of transcription from RNA polymerase II promoter"
unique(res[grep(goterm,res$name_1006,ignore.case =TRUE),"external_gene_name"])


