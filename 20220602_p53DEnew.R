library(tximport)
library(readr)
library(tximportData) #this is dummy data just in case your code broke and you really need help
library(DESeq2)
library(biomaRt)
library(ggplot2)
library(tximeta)
library(devtools)
library(tidyverse)

#provide file path to quant files
dir <- "/Volumes/LaCie"

#read in sample table to select quant files from dir
samples <- read.csv(file.path("/Volumes/LaCie/sampleinfo_p53signiture.csv"), header= TRUE)
rownames(samples) <- samples$Sample
head(samples)

files <- file.path(dir, samples$Name,"salmon.quant", "quant.sf")
names(files) <- samples$Name
head(files)

#get data for tx2gene
mart <- useMart(biomart = 'ensembl', dataset = 'mmusculus_gene_ensembl')
mart_res <- getBM(attributes = c('ensembl_transcript_id', 'external_gene_name', "chromosome_name"), mart = mart)
head(mart_res)

txi <- tximport(files, type="salmon", tx2gene = mart_res)

ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = samples,
                                   design = ~ P53)

keep <- rowSums(counts(ddsTxi)) >= 10 #prefilter counts<10
ddsTxi <- ddsTxi[keep,]
ddsTxi <- DESeq(ddsTxi)
res <-results(ddsTxi)
res

sum(res$padj < 0.1, na.rm=TRUE)
vsd <- varianceStabilizingTransformation(ddsTxi)
scaled_vsd <- t(scale(t(assay(vsd)), center = T, scale = T) )
heatmap(scaled_vsd, cexCol = .5, xlab = "Samples")
scaled_vsd <- na.omit(scaled_vsd)
Heatmap(scaled_vsd, row_km = 4)

resultsNames(ddsTxi) 
shrunk <- lfcShrink(res,type = 'ashr') 
DESeq2::plotMA(shrunk, coef = 2, ylim = c(-5,5))

