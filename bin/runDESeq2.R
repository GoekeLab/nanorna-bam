#!/usr/bin/env Rscript

if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager", repos='http://cran.us.r-project.org')
 }
if (!require("DESeq2")){
    BiocManager::install("DESeq2",update = FALSE, ask= FALSE)
    library(DESeq2)
}


args = commandArgs(trailingOnly=TRUE)
if (length(args) < 3) {
  stop("Please input the directory with the featureCounts results and the sample information file", call.=FALSE)
} else if (length(args)==3) {
  # default output file
  args[4] = "DESeq2out.txt"
}
#DeSeq2
transcriptquant <- args[1]
path<-args[2]
#create a dataframe for all samples
if (transcriptquant == "stringtie"){
  #count_files<- grep(list.files(path), pattern='tx_', inv=TRUE, value=TRUE)
  countTab <- data.frame(read.table(dir(path, pattern = "counts_gene.txt$", full.names = TRUE),sep="\t",header=TRUE, skip = 1))
  countTab$Chr <- countTab$Start <- countTab$End <- countTab$Length <- countTab$Strand <- NULL
}
if (transcriptquant == "bambu"){
 
  countTab <- data.frame(read.table(dir(path, pattern = "counts_gene", full.names = TRUE),sep="\t",header=TRUE))
}

#sampInfo <- read.csv("~/Downloads/nanorna-bam-master/two_conditions.csv",row.names = 1)
sampInfo<-read.csv(args[3],row.names=1)
#all(rownames(sampInfo) %in% colnames(countTab))
#all(rownames(sampInfo) == colnames(countTab))
dds <- DESeqDataSetFromMatrix(countData = countTab,colData = sampInfo,design = ~ condition)
dds <- DESeq(dds)
res <- results(dds)
#register(MulticoreParam(6))
resOrdered <- res[order(res$pvalue),]
write.csv(as.data.frame(resOrdered), file=args[4])
