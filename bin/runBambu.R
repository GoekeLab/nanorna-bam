#!/usr/bin/env Rscript

##install bambu if bambu is not installed
##load bambu
if (!requireNamespace("devtools", quietly = TRUE)){
    install.packages("devtools", repos='http://cran.us.r-project.org')
}
if (!require("bambu")){
    devtools::install_github("GoekeLab/bambu")
    library(bambu)
}

args = commandArgs(trailingOnly=TRUE)
if (length(args) < 3) {
  stop("Please input the fullpath for the present directory and the sample sheet.", call.=FALSE)
} else if (length(args)==3) {
  if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager", repos='http://cran.us.r-project.org')
  }
  if (!require("BSgenome.Hsapiens.NCBI.GRCh38")){
    BiocManager::install("BSgenome.Hsapiens.NCBI.GRCh38",update = FALSE, ask= FALSE)
    library(BSgenome.Hsapiens.NCBI.GRCh38)
  }
  genomeseq <- "BSgenome.Hsapiens.NCBI.GRCh38"  #use BSgenome if fasta file is not provided
} else {
  genomeseq <- args[4]
}

pwd <- args[1]
sampInfo <- read.table(args[2], header = TRUE, sep=",")
annotation_list <- unique(sampInfo$transcriptome)
output_tag <- args[3]
for(i in 1:length(annotation_list)){
    readlist <- sampInfo[sampInfo$transcriptome %in% annotation_list[i]]$bam
    annot_gtf <- paste(pwd,annotation_list[i],sep="/")
    readlist <- paste(pwd,readlist,sep="/")
    grlist <- prepareAnnotationsFromGTF(annot_gtf)
    se <- bambu(reads = readlist, annotations=grlist,genomeSequence = genomeseq)
    writeBambuOutput(se,output_tag)
}
