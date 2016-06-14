## Amy Olex
## 6/13/16
## Script to calculate PCA of 2 VCF files - one before batch correction and one after.
## It will then calculate the correlation of the batches to the Principle Components to see if there was improvements.
## This script is meant to run using Rscript command on the terminal
##
## USAGE:  >> Rscript PCA_batch_correlation_061316.R <~/working/dir> <before.vcf> <after.vcf> <batches.file> <out.prefix>


args <- commandArgs(TRUE)
#args <- c("~/Desktop/CCTR/Data/MikhailDozmorov_VCUSamples_103015/VCF Batch Removal 061316", "Before-After VCFs/VCU22-chr22-merged.recode.vcf", "Before-After VCFs/VCU22-chr22-merged-BatchEffFiltered.recode.vcf", "VCU22-SampleBatches.txt", "VCU22_PCAquantification_chr22")

setwd(args[1])
vcf1.fn <- args[2]
vcf2.fn <- args[3]
batch.file <- args[4]
out.prefix <- args[5]  ## like "TCGA_test24_chr22"

#setwd("~/Desktop/CCTR/Data/MikhailDozmorov_VCUSamples_103015/TCGA test 24 cohort/SNP-subtracted VCFs")
library("SNPRelate")
library("calibrate")
library("plyr")


#####
#### Import Batches
#####
batches = read.delim(file=batch.file, row.names=1, header=TRUE)
colors <- as.numeric(batches$Batch)

########
#### Get PCA for first VCF file
########
snpgdsVCF2GDS(vcf1.fn, "ccm1.gds",  method="biallelic.only")
genofile1 <- openfn.gds("ccm1.gds")
ccm_pca1<-snpgdsPCA(genofile1)
ids1 <- ccm_pca1$sample.id
ordered_batches <- batches[ids1,]
PCs1 <- ccm_pca1$eigenvect

########
#### Get PCA for first VCF file
########
snpgdsVCF2GDS(vcf2.fn, "ccm2.gds",  method="biallelic.only")
genofile2 <- openfn.gds("ccm2.gds")
ccm_pca2<-snpgdsPCA(genofile2)
ids2 <- ccm_pca2$sample.id
ordered_batches <- batches[ids2,]
PCs2 <- ccm_pca2$eigenvect

## Now calculate the correlation of each PC and the ordered_batches
PCs1_b <- cbind(ordered_batches, PCs1)
colnames(PCs1_b) <- c("batches",paste("PC",1:dim(PCs1)[2], sep=""))
pval = matrix(ncol=4, nrow=dim(PCs1_b)[1])
colnames(pval) <- c("Before_pval","Before_cor", "After_pval", "After_cor")
for(i in 1:(length(colnames(PCs1_b))-1)){
  
  form <- as.formula(paste("PC",i," ~ batches", sep=""))
  pval[i,1] <- summary(lm(form, as.data.frame(PCs1_b)))$coefficients[2,4]
  pval[i,2] <- cor(PCs1_b[,"batches"], PCs1_b[,paste("PC",i,sep="")])
}

PCs2_b <- cbind(ordered_batches, PCs2)
colnames(PCs2_b) <- c("batches",paste("PC",1:dim(PCs2)[2], sep=""))
for(i in 1:(length(colnames(PCs2_b))-1)){
  
  form <- as.formula(paste("PC",i," ~ batches", sep=""))
  pval[i,3] <- summary(lm(form, as.data.frame(PCs2_b)))$coefficients[2,4]
  pval[i,4] <- cor(PCs2_b[,"batches"], PCs2_b[,paste("PC",i,sep="")])
}


outfile=paste(out.prefix, "_PCA_Batch_Correlations_Matrix.txt", sep="")
write.table(pval, file=outfile, quote=FALSE)

