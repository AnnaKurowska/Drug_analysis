setwd('/Users/Ania/Desktop/Szkoła/ENS /2nd_year/Stage/Practical/EMBL')
getwd()

EMBL <- read.csv(file = 'ViraminD_GSE94138_series_matrix.txt',sep = "\t")

### Retrieving genes for probes
BiocManager::install("illuminaHumanv4.db")
library(illuminaHumanv4.db)
xx <- as.list(illuminaHumanv4ALIAS2PROBE)
wanted_genes <- c("B2M", "HBD", "HBD-2","HBD-3","HBD-4","LL37", "REG3A", "PLA2G2A", "IL1B", "IL8", "TNF")

identifiers_genes <- xx[wanted_genes]
identifiers_genes_unlisted <-unlist(identifiers_genes)
identifiers_genes_unlisted <- as.data.frame(identifiers_genes_unlisted)
identifiers_to_retrieve <- identifiers_genes_unlisted$identifiers_genes_unlisted
EMBL_selected <- EMBL[EMBL$ID_REF %in% identifiers_to_retrieve,] 
## don't know how to make it loop through both 

# to conduct differential analysis:
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = colData,
                              design = ~condition)


## DE analysis
BiocManager::install("DESeq2")
library(DESeq2)
library(dplyr)
library(reshape2)
library(ggplot2)
setwd("/Users/Ania/Desktop/Szkoła/ENS /2nd_year/Stage/Practical/try_out_datasets")
macrophages <-read.csv(file = 'GSE56583_series_matrix.txt',sep = "\t")
macrophages_ID_REF <- macrophages$ID_REF
macrophages <- dplyr::select(macrophages, 2:ncol(macrophages))
rownames(macrophages) <- macrophages_ID_REF

#replace NA values with 0s
which(is.na(macrophages), arr.ind=TRUE)
macrophages <- macrophages[1:nrow(macrophages)-1,] # last row wasn't part of the matrix
macrophages <- as.integer(macrophages)
macrophages_integers <-mutate_all(macrophages, function(x) as.integer(x))
#forget about this section
# pre <- c()
# post <- c()
# rownamepre <- c()
# rownamepost <- c()
# for (i in 1:22){
#   pre[i] <- "pre"
#   post[i] <- "post"
#   rownamepre[i] <- paste0("Macrophage-VitaminD-Subject", i, "Pre")
#   rownamepost[i] <- paste0("Macrophage-VitaminD-Subject", i, "Post")
# }

macrophages_design <- data.frame(Conditions = c(pre, post))
rownames(macrophages_design) <- colnames(macrophages)

dds <- DESeqDataSetFromMatrix(countData = macrophages_integers,
                              colData = macrophages_design,
                              design = ~Conditions)

colData(dds)$Conditions <- factor(colData(dds)$Conditions,
                                  levels=c("pre","post"))
#Normalization

dds <- estimateSizeFactors(dds)
sizeFactors(dds)

dds <- estimateDispersions(dds,fitType= 'mean')
dds <- nbinomWaldTest(dds)
summary(results(dds))
res <- results(dds)
head(res)
mcols(res, use.names=TRUE)
res <- as.data.frame(res)
# p value distribution
ggplot(res, aes(x = pvalue)) +
  geom_histogram(bins = 100) + 
  ggtitle(label = "P-value distibution")
# so it this analysis was done right, the differences are very tiny between these groups

post_values <-apply(macrophages[,22:ncol(macrophages)],1, mean)
pre_values <-apply(macrophages[,1:22],1, mean) # values in pre and post conditions are simply too similar!
