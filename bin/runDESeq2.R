#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args) < 2) {
  stop("Please input the directory with the featureCounts results and the sample information file", call.=FALSE)
} else if (length(args)==2) {
  # default output file
  args[3] = "out.txt"
}
#DeSeq2
library("DESeq2")
library("BiocParallel")
path<-args[1]
count_files<- grep(list.files(path), pattern='tx_', inv=T, value=T)
#create a dataframe for all samples 
fullpath<-paste(path,count_files[1],sep='/') 
df <- data.frame(read.table(fullpath,sep="\t",header=T)[,c(1,7)])
for(i in 2:length(count_files)){
  fullpath<-paste(path,count_files[i],sep='/') 
  samp_df <- read.table(fullpath,sep="\t",header=T)[,c(1,7)]
  df<- merge(df,samp_df,by="Geneid",all=TRUE)
}
countTab <- df[,-1]
rownames(countTab) <- df[,1]
sampInfo<-read.csv(args[2],row.names=1)
all(rownames(sampInfo) %in% colnames(countTab))
all(rownames(sampInfo) == colnames(countTab))
dds <- DESeqDataSetFromMatrix(countData = countTab,colData = sampInfo,design = ~ condition)
dds <- DESeq(dds)
res <- results(dds)
register(MulticoreParam(4))
resOrdered <- res[order(res$pvalue),]
write.csv(as.data.frame(resOrdered), file=args[3])
