library(biomaRt)
library(GO.db)
library(limma)

#citation("mogene10sttranscriptcluster.db")
rm_a <- function(x){
  return(x[1])
}
#################
dat <- read.table(as.matrix("~/Documents/SalkProjects/ME/ShortLongSingature/raw/GSE55591-GPL1261_series_matrix2.txt"),row.names = 1)
probes <- rownames(dat)
#### Set up probes
x <- mogene10sttranscriptclusterACCNUM
mapped_probes <- mappedkeys(x)
xx <- as.list(x[mapped_probes])
#### Set up BiomaRt
mart <- useEnsembl(dataset = "mmusculus_gene_ensembl", biomart = "ensembl", version = 83)
#filter <- c("affy_mouse430a_2") #TRY: listFilters(ensembl)
#attrib <- c("affy_mouse430a_2","external_gene_name")
filter <- c("affy_moe430a") #TRY: listFilters(ensembl)
attrib <- c("affy_mouse430_2","external_gene_name")
res = getBM(attributes=attrib,filters=filter,values=probes,mart=mart)
#####
meta <- data.frame(region = c(rep("MGE",12),rep("CTX",12)),
                   mouse = c(rep(c("WT","KO"),each = 3)),
                   time = c(0,1,6)
                   )

######
conditional <- region == "CTX" & meta$mouse == "WT" & meta$time == 0 | region == "CTX" & meta$mouse == "WT" & meta$time == 1
dat2 <- dat[,conditional]
met2 <- meta[conditional,]
met2$time <- as.factor(met2$time)
design <- model.matrix(~0+met2$time)
colnames(design) <- levels(met2$time)
fit <- lmFit(dat2, design)
fit <- eBayes(fit)
a <-topTable(fit)
a$gene <- res[match(rownames(a), res$affy_mouse430_2),"external_gene_name"]

p <- as.data.frame(fit$p.value)
p <- na.exclude(p[order(p$`1`),])
p$gene <- res[match(rownames(p), res$affy_mouse430_2),"external_gene_name"]
p$f <- p.adjust(p = p$`1`)

fit$p.value["1423100_at",]

res[res$affy_mouse430_2 == "1438625_s_at",]
