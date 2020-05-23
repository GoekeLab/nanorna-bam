#!/usr/bin/env Rscript

if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager", repos='http://cran.us.r-project.org')
 }
if (!require("DESeq2")){
    BiocManager::install("DESeq2",update = FALSE, ask= FALSE)
    library(DESeq2)
}


args = commandArgs(trailingOnly=TRUE)
if (length(args) < 2) {
  stop("Please input the directory with the featureCounts results and the sample information file", call.=FALSE)
} else if (length(args)==2) {
  # default output file
  args[3] = "out.txt"
}
#DeSeq2
path<-args[1]
count_files<- grep(list.files(path), pattern='tx_', inv=T, value=T)
#create a dataframe for all samples 
fullpath<-paste(path,count_files[1],sep='/')
count.matrix <- data.frame(read.table(fullpath,sep="\t",header=T)[,c(1,7,8)])
for(i in 2:length(count_files)){
  fullpath<-paste(path,count_files[i],sep='/')
  samp_df <- read.table(fullpath,sep="\t",header=T)[,c(1,8)]
  count.matrix<- merge(count.matrix,samp_df,by="Geneid",all=TRUE)
}
df <- count.matrix[,-1]
rownames(df) <- count.matrix[,1]
df <- aggregate(. ~ gene_id, data=df, FUN=sum)
countTab <- df[,-1]
rownames(countTab) <- df[,1]
sampInfo<-read.csv(args[2],row.names=1)
all(rownames(sampInfo) %in% colnames(countTab))
all(rownames(sampInfo) == colnames(countTab))
dds <- DESeqDataSetFromMatrix(countData = countTab,colData = sampInfo,design = ~ condition)
dds <- DESeq(dds)
res <- results(dds)
register(MulticoreParam(6))
resOrdered <- res[order(res$pvalue),]
write.csv(as.data.frame(resOrdered), file=args[3])
