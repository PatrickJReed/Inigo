#######
# Retrotransposons in hippocampus
#######
# load raw data from the normal load data for short term r script
#save(list = c("sine.cor.est"),file = "~/Documents/SalkProjects/ME/SINE/SINER/proxC_sine_cors.rda",compress = TRUE)
######
library(ggplot2)
library(reshape)
library(edgeR)
exact.sml <- function(dat, variable, Pair, Lib){
  ########################
  ### Pre-process Data
  ########################
  cds <- DGEList(dat,group=group,lib.size = Lib)
  # Filter out really lowly detected reads (*from EdgeR tutorial)
  cds <- cds[rowSums(1e+06 * cds$counts/expandAsMatrix(cds$samples$lib.size, dim(cds)) > 1) >= 3, ]
  cds <- calcNormFactors( cds )
  ########################
  ### Estimate parameters from data
  ########################
  cds <- estimateCommonDisp( cds )
  cds <- estimateTagwiseDisp( cds )
  cds <- estimateTrendedDisp(cds)
  ## Run the test
  de.cmn <- exactTest( cds , pair = Pair)
  ########################
  ### Format Results
  ########################
  res <- as.data.frame(de.cmn$table)
  res <- res[order(res$PValue),]
  #f <- fdrtool(x=res$PValue,statistic="pvalue",plot=FALSE)
  res$f <- p.adjust(res$PValue, method = "fdr")
  return(res)
}

##########
## Run EdgeR
##########

#####
condition <- which(sine_col_meta$Brain_Region == "DG" & sine_col_meta$value == "tso_elements" & sine_col_meta$Mouse_condition == "EE" & sine_col_meta$FOS != "L" )
dat <- sine_count[,condition] 
met1 <- sine_col_meta[condition,]
dat <- dat[rowSums((dat)) > 0,]
Lib <- sine_col_meta[condition,"alignable"]
group <- met1$FOS
Pair <- levels(as.factor(as.character(group)))

res <- exact.sml(dat, group, Pair, Lib)
##########
## Plot Cummulative results
##########

res$Class <- FALSE
res[grep("B2_",rownames(res)),"Class"] <- "B2"
res[grep("B1_",rownames(res)),"Class"] <- "B1"
res[grep("ID_",rownames(res)),"Class"] <- "ID"
res[grep("L1",rownames(res)),"Class"] <- "L1"
res[res$Class == "FALSE","Class"] <- "misc"

ggplot(res, aes(-logFC, -log(PValue), colour = Class))+
  geom_point(size = 2)+
  xlim(c(-3,3))+
  theme_bw(base_size = 12)+
  geom_hline(yintercept = -log(0.03), linetype  = "dashed")+
  labs(title = "NE Negs\nRepeat Expression")

##########
## Plot SINEs one at a time
##########
dat2 <- sine_tpm[,sine_col_meta$value == "tso_elements"] 
met2 <- sine_col_meta[sine_col_meta$value == "tso_elements",]
tmp <- melt(t(dat2))
tmp$Mouse_condition <- met2$Mouse_condition
tmp$Brain_region <- met2$Brain_Region
tmp$FOS <- met2$FOS
tmp$PROX1 <- met2$PROX1
tmp$CTIP2 <- met2$CTIP2
tmp$element <- as.character(rep(sine_row_meta$V1,each = ncol(dat2)))
tmp$tso <- met2$value 

TE <- "L1MdA_IV"
ggplot(na.exclude(tmp[tmp$element == TE & tmp$value !=0 & tmp$value < 1000, ]), aes(FOS,value))+
  geom_violin()+
  geom_point()+
  facet_grid(Mouse_condition ~ Brain_region)+
  labs(title = TE)


##########
## Correlations with SINEs
##########

samples <- unique(sine_col_meta[sine_col_meta$value == "tso_elements","Sample_ID"] )
a <- table(c(colnames(tpmProxC), samples,samples))
samples <- names(a[a==3])
dat.sine <- sine_tpm[,match(samples, sine_col_meta$Sample_ID)] 
dat.gene <- na.exclude(tpmProxC[,samples])

sine.cor.est <- vector()
nm <- rownames(res)
for(i in nm){
  sine <- as.numeric(dat.sine[i,])
  tmp.cor <- 
  sine.cor.est <- cbind(sine.cor.est, tmp.cor[,1])
}

colnames(sine.cor.est) <- nm[1:ncol(sine.cor.est)]
#colnames(sine.cor.p) <- nm[1:ncol(sine.cor.p)]
rownames(sine.cor.est) <- rownames(dat.gene)
#rownames(sine.cor.p) <- rownames(dat.gene)

####### Plot some of the results
a <- rownames(na.exclude(sine.cor.est[abs(sine.cor.est$B2_Mm1a) > 0.3,]))

gene <- "Snf8"
te <- "B2_Mm1t"
tmp <- data.frame(gene = as.numeric(dat.gene[gene,]),
                  te = as.numeric(dat.sine[te,]),
                  sine_col_meta[samples,])
ggplot(tmp[tmp$te != 0,], aes(gene, te, colour = Brain_Region))+
  geom_point()+
  geom_smooth(method = "lm")+
  xlab(gene)+
  ylab(te)

