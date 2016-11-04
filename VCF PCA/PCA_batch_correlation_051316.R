## Amy Olex
## 5/13/16
## Script to calculate PCA of a VCF file.
## This script is meant to run using Rscript command on the terminal
##
## USAGE:  >> Rscript create_VCF_PCA_plots_051316.R <~/working/dir> <infile.vcf> <batches.file> <out.prefix> <pval.cutoff>


args <- commandArgs(TRUE)
setwd(args[1])
vcf.fn <- args[2]
batch.file <- args[3]
out.prefix <- args[4]  ## like "TCGA_test24_chr22"
pval.cutoff <- args[5]

#setwd("~/Desktop/CCTR/Data/MikhailDozmorov_VCUSamples_103015/TCGA test 24 cohort/SNP-subtracted VCFs")
library("SNPRelate")
library("calibrate")
library("plyr")


#vcf.fn<-"TCGA_test24Cohort_merged.chr22.LowBatchCorrSNPs.recode.vcf"
snpgdsVCF2GDS(vcf.fn, "ccm.gds",  method="biallelic.only")
genofile <- openfn.gds("ccm.gds")
ccm_pca<-snpgdsPCA(genofile)

#ids <- sapply(strsplit(ccm_pca$sample.id, split="_"),"[",1)
#colors <- revalue(factor(sapply(strsplit(ccm_pca$sample.id, split="_"),"[",2)),c("MERGED"=1,"val.pe.sam.sort"=2) )
#colors <- revalue(factor(c(2012,2013,2014,2012,2013,2013,2014,2012,2013,2012,2012,2014,2014,2014,2012,2013,2014,2013,2011,2011,2011,2011,2011,2011)), c("2011"=1,"2012"=2,"2013"=3,"2014"=4))
#NOTE, somehow in the imputed file the columns got moved around, so I created the new color list.
ids <- ccm_pca$sample.id
batches = read.delim(file=batch.file, row.names=1, header=TRUE)
batches2 <- factor(t(as.data.frame(t(batches))[,ids])[,1])
batch = revalue(factor(batches2), c("2011"=1,"2012"=2,"2013"=3,"2014"=4))
colors <- as.numeric(batch)

outfile=paste(out.prefix, "_LowBatchCorrSNPs-pval", pval.cutoff,"-PCA.png", sep="")
png(filename=outfile, width=70, height=150, res=300, units="mm", bg="white", pointsize=9)

par(mfrow=c(3,1),mar=c(3,3,2,2),mgp=c(2,1,0))
plot(ccm_pca$eigenvect[,1],ccm_pca$eigenvect[,2],pch=20, col=colors,xlim=c(-1,1),ylim=c(-1,1),xlab="PC1",ylab="PC2")
#textxy(ccm_pca$eigenvect[,1],ccm_pca$eigenvect[,2],ids, offset=0)
legend(-1,1,c("2011","2012","2013","2014"), pch=20, col=c(1,2,3,4),cex=.7)
title(paste("PCA of LowBatchCorr SNPs of ", out.prefix, "\n Pval cutoff = ", pval.cutoff, sep=""))

plot(ccm_pca$eigenvect[,1],ccm_pca$eigenvect[,3] , pch=20, col=colors, xlim=c(-1,1),ylim=c(-1,1),xlab="PC1",ylab="PC3")
#textxy(ccm_pca$eigenvect[,1],ccm_pca$eigenvect[,3],ids, offset=2)

plot(ccm_pca$eigenvect[,2],ccm_pca$eigenvect[,3] , pch=20, col=colors, xlim=c(-1,1),ylim=c(-1,1),xlab="PC2",ylab="PC3")
#textxy(ccm_pca$eigenvect[,2],ccm_pca$eigenvect[,3],ids, offset=2)

dev.off()


