#!/usr/bin/env Rscript

if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager", repos='http://cran.us.r-project.org')
 }
if (!require("DRIMSeq")){
    BiocManager::install("DRIMSeq",update = FALSE, ask= FALSE)
    library(DRIMSeq)
}
if (!require("DEXSeq")){
    BiocManager::install("DEXSeq",update = FALSE, ask= FALSE)
    library(DEXSeq)
}
if (!require("stageR")){
    BiocManager::install("stageR",update = FALSE, ask= FALSE)
    library(stageR)
}


args = commandArgs(trailingOnly=TRUE)
if (length(args) < 3) {
  stop("Please input the directory with the transcript level featureCounts results and the sample information file", call.=FALSE)
} 
######### import data ###########
##Import featureCounts data
#path <- "~/Downloads/nanorna-bam-master/mod/featureCounts_transcript/"
transcriptquant <- args[1]
path<-args[2]
if (transcriptquant == "stringtie"){
  count.matrix <- data.frame(read.table(dir(path, pattern = "counts_transcript.txt$", full.names = TRUE),sep="\t",header=TRUE, skip = 1))
  count.matrix$Chr <- count.matrix$Start <- count.matrix$End <- count.matrix$Length <- count.matrix$Strand <- NULL
}

if (transcriptquant == "bambu"){
  count.matrix <- data.frame(read.table(dir(path, pattern = "counts_transcript", full.names = TRUE),sep="\t",header=T))
}
colnames(count.matrix)[1] <- "feature_id"
colnames(count.matrix)[2] <- "gene_id"
#sample information
#sampleinfo<-read.table("~/Downloads/nanorna-bam-master/samples_conditions.csv",sep=",",header=T)
sampleinfo<-read.table(args[3],sep=",",header=T)
colnames(sampleinfo)[1] <- "sample_id"
condition_names <- sampleinfo[!duplicated(sampleinfo$condition),]$condition
#condition_names <- c(levels(sampleinfo$condition))
lgcolName <- "log2fold"
for (i in length(condition_names):1){
  lgcolName <- paste(lgcolName,condition_names[i],sep='_')
}
######### Filtering #############

cat(paste(sampleinfo$sample_id, collapse = " "))
cat(paste(colnames(count.matrix), collapse = " "))
d <- dmDSdata(counts=count.matrix, samples=sampleinfo)
# include genes expressed in minimal min_samps_gene_expr samples with min_gene_expr
# include transcripts expressed in min_samps_feature_expr samples with min_feature_expr;  
# include transcripts expressed in min_samps_feature_prop samples with min_feature_prop;
n_samp_gene <- length(sampleinfo$sample_id)/2
n_samp_feature <- length(sampleinfo$sample_id)/2
min_count_gene <- 5 
min_count_feature <- 5
dFilter <- dmFilter(d,
                    min_samps_feature_expr = n_samp_feature, 
                    min_samps_feature_prop = n_samp_feature,
                    min_samps_gene_expr = n_samp_gene, 
                    min_feature_expr = min_count_feature,
                    min_gene_expr = min_count_gene,
                    min_feature_prop=0.1)

########## DEXSeq #########
formulaFullModel <- as.formula("~sample + exon + condition:exon")
#formulaReducedModel <- as.formula("~sample + exon + covariate:exon")
dxd <- DEXSeqDataSet(countData=round(as.matrix(counts(dFilter)[,-c(1:2)])),
                     sampleData=DRIMSeq::samples(dFilter),
                     design=formulaFullModel,
                     featureID = counts(dFilter)$feature_id,
                     groupID=counts(dFilter)$gene_id)

dxd <- estimateSizeFactors(dxd)
print('Size factor estimated')
dxd <- estimateDispersions(dxd, formula = formulaFullModel)
print('Dispersion estimated')
#reducedModel = formulaReducedModel, 
dxd <- testForDEU(dxd, fullModel = formulaFullModel)
print('DEU tested')
dxd <- estimateExonFoldChanges(dxd)
print('Exon fold changes estimated')
##Got an error with the MulticoreParam on Docker, but the MulticoreParam was fine locally.
#system.time({
 # BPPARAM <- MulticoreParam(6)
#  dxd <- estimateSizeFactors(dxd)
#  print('Size factor estimated')
#  dxd <- estimateDispersions(dxd, formula = formulaFullModel, BPPARAM=BPPARAM)
#  print('Dispersion estimated')
#  dxd <- testForDEU(dxd, fullModel = formulaFullModel, BPPARAM=BPPARAM)
#  print('DEU tested')
#  dxd <- estimateExonFoldChanges(dxd)
 # dxd <- estimateExonFoldChanges(dxd, BPPARAM=BPPARAM)
#  print('Exon fold changes estimated')
#}) 
# Extract DEXSeq results 
dxr <- DEXSeqResults(dxd, independentFiltering=FALSE)

### Pinpoint genes with transcripts differentially used ### 
##### and transcripts that are differentially used ########
library(stageR)
strp <- function(x) substr(x,1,15)
qval <- perGeneQValue(dxr)
dxr.g <- data.frame(gene=names(qval),qval)

columns <- c("featureID","groupID","pvalue",lgcolName)
dxr_pval <- as.data.frame(dxr[,columns])
#head(dxr_pval)
pConfirmation <- matrix(dxr_pval$pvalue,ncol=1)
dimnames(pConfirmation) <- list(strp(dxr_pval$featureID),"transcript")
pScreen <- qval
names(pScreen) <- strp(names(pScreen))
tx2gene <- as.data.frame(dxr_pval[,c("featureID", "groupID")])
for (i in 1:2) tx2gene[,i] <- strp(tx2gene[,i])

stageRObj <- stageRTx(pScreen=pScreen, pConfirmation=pConfirmation,
                      pScreenAdjusted=TRUE, tx2gene=tx2gene)
stageRObj <- stageWiseAdjustment(stageRObj, method="dtu", alpha=0.05)
suppressWarnings({
  dex.padj <- getAdjustedPValues(stageRObj, order=FALSE,
                                 onlySignificantGenes=TRUE)
})
write.csv(dxr_pval, file="DEXseq_out.txt")
