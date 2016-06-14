## Amy Olex
## 6/13/16
## Script to calculate PCA of a VCF file.
## This script is meant to run using Rscript command on the terminal
##
## USAGE:  >> Rscript create_VCF_PCA_plots_051316.R <~/working/dir> <infile.vcf> <batches.file> <out.prefix> <pval.cutoff>


#args <- commandArgs(TRUE)
args <- c("/Users/alolex/Desktop/CCTR/Data/MikhailDozmorov_VCUSamples_103015/VCF Batch Removal 061316", "Before-After VCFs/VCU22-chr22-merged.recode.vcf", "VCU22-SampleBatches.txt", "mytest", "0.05")
setwd(args[1])
vcf.fn <- args[2]
batch.file <- args[3]
out.prefix <- args[4]  ## like "TCGA_test24_chr22_LowBatchCorrSNPs"
pval.cutoff <- args[5]

#setwd("~/Desktop/CCTR/Data/MikhailDozmorov_VCUSamples_103015/TCGA test 24 cohort/SNP-subtracted VCFs")
library("SNPRelate")
library("calibrate")
library("plyr")


#vcf.fn<-"TCGA_test24Cohort_merged.chr22.LowBatchCorrSNPs.recode.vcf"
snpgdsVCF2GDS(vcf.fn, "ccm.gds",  method="biallelic.only")
genofile <- openfn.gds("ccm.gds")
ccm_pca<-snpgdsPCA(genofile)
ids <- ccm_pca$sample.id

tmp = read.delim(file=batch.file, row.names=1, header=TRUE)
batches <- tmp[ids,]
colors <- as.numeric(batches$Batch)
legend.labels <- unique(batches[,1:2])$Label
legend.cols <- unique(batches[,1:2])$Batch


outfile=paste(out.prefix, "-pval", pval.cutoff,"-PCA.png", sep="")
png(filename=outfile, width=70, height=150, res=300, units="mm", bg="white", pointsize=9)

par(mfrow=c(3,1),mar=c(3,3,2,2),mgp=c(2,1,0))
plot(ccm_pca$eigenvect[,1],ccm_pca$eigenvect[,2],pch=20, col=colors,xlim=c(-1,1),ylim=c(-1,1),xlab="PC1",ylab="PC2")
#textxy(ccm_pca$eigenvect[,1],ccm_pca$eigenvect[,2],ids, offset=0)
legend(-1,1,legend.labels, pch=20, col=legend.cols,cex=.7)
title(paste("PCA of ", out.prefix, "\n Pval cutoff = ", pval.cutoff, sep=""))

plot(ccm_pca$eigenvect[,1],ccm_pca$eigenvect[,3] , pch=20, col=colors, xlim=c(-1,1),ylim=c(-1,1),xlab="PC1",ylab="PC3")
#textxy(ccm_pca$eigenvect[,1],ccm_pca$eigenvect[,3],ids, offset=2)

plot(ccm_pca$eigenvect[,2],ccm_pca$eigenvect[,3] , pch=20, col=colors, xlim=c(-1,1),ylim=c(-1,1),xlab="PC2",ylab="PC3")
#textxy(ccm_pca$eigenvect[,2],ccm_pca$eigenvect[,3],ids, offset=2)

dev.off()


