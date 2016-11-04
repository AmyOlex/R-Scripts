## Sample Code for DESeq2 Analysis
## Amy Olex
## 6/27/16
## 
## Please Cite/Acknowledge the CTSA Award No. UL1TR000058 if this code aided in generating data analysis that was included in a publication.
## Award acknowledgement is the only way we can track our performance in facilitating research.
## More information can be found at http://www.cctr.vcu.edu/informatics/biomedical/cite.html
## 
## DISCLAIMER!!  This code is generic and only an example.  Once comfortable you will need to modify it to suit your individual analysis needs!
## Please DO NOT use this code as-is.


##Set working directory
setwd("~/Desktop/CCTR/Data/GloriaMuday_Arabidopsis_0814/100914_RNASeq_gene_exp_analysis")

##Load needed libraries
library("DESeq2")

##import sample table
sample_table <- read.table("/Users/alolex/Desktop/CCTR/Data/GloriaMuday_Arabidopsis_0814/original_htseq-count_files/to_merge3.txt", sep="\t",header=TRUE)

##check to see which condition is listed first
##the condition listed first will be used as the reference/control sample in the DE analysis.
levels(sample_table$condition)

##import data from HTSeq-count files
dds <- DESeqDataSetFromHTSeqCount(sample_table, directory="/Users/alolex/Desktop/CCTR/Data/GloriaMuday_Arabidopsis_0814/original_htseq-count_files", design = ~condition)

##can type in just the name to view a summary of the object
dds

##Run a DESeq2 analysis, which includes normalization.  This step can take A WHILE depending on how much data you are running.
dds <- DESeq(dds)

##Get the results from the analysis (type in "res" at the command prompt to see a summary)
res <- results(dds)

##Get the significantly expressed genes based on adjusted p-value.
resSig <- as.data.frame(res[which(res$padj < 0.01),])

##Write gene list to file
write.table(resSig, file="my_sig_gene_list.txt", quote=FALSE, sep="\t")

##Create a histogram of the significant expression fold changes.
##After filtering out all the insignificant genes you will most likely see a bi-modal distribution centered around zero.
plot(hist(resSig$log2FoldChange, breaks=100), main="My Histogram", xlab="Fold Change", ylab="Number of Genes")
