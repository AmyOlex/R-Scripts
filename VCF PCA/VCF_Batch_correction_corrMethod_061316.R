## Amy Olex
## 6/13/16
## This script is the first of 3 associate with removing batch effects in VCF files
## by calculating the correlation of SNPs called to a given set of batches.
## Any SNP found to NOT significantly correlate with the given batches is written to a text file.
## The user then needs to use this output file to extract the SNPs from the original VCF file.
## Once you have your new VCF files you can move on to performing the PCA analysis and creating the plots.
## This script is meant to run using Rscript command on the terminal
##
## USAGE: >> Rscript VCF_remove_batch_corrMethod.R <infile.vcf> <sampleBatchFile.txt>

library(VariantAnnotation)
library("plyr")
args <- commandArgs(TRUE)
print(args)
#setwd("/var/bioinformatics-tcga-projects/alolex/clients/MDozmorov/BRCA_WXS_VCU_Wenhu/PREPROCESSING/6_vcfBatchCorrection")
fl <- args[1]
f2 <- args[2]

vcf <- readVcf(fl, "hg19")
GT <- geno(vcf)$GT
rownames(GT) <- make.names(rownames(GT), allow_=F)
colnames(GT) <- make.names(colnames(GT), allow_=F)
# class(GT[,1]) is a character

batches = read.delim(file=f2, row.names=1, header=TRUE)
colors <- as.numeric(batches$Batch)

#sort matrix based on batches order, now all columns are factors.
GT_sorted <- GT[,row.names(batches)]
# class(GT_sorted[,1]) is a character

#revalue everything to 0 = missing genotype, 1 = genotype present
set_to_ones = GT_sorted != "."
set_to_zeros = GT_sorted == "."

GT_sorted_reval <- GT_sorted
GT_sorted_reval[set_to_ones] = 1
GT_sorted_reval[set_to_zeros] = 0
# class(GT_sorted_reval[,1]) is a character
# Convert to numeric
class(GT_sorted_reval) <- "numeric"

#transpose
GT_tr <- t(GT_sorted_reval)
# class(GT_tr[,1]) is now numeric

#Add batch to end
GT_tr <- cbind(GT_tr, batch)
GT_tr_preservation <- GT_tr

#GT_tr <- GT_tr[,1:100]
#GT_tr <- cbind(GT_tr, batch)

# calculating the correlation/linear relationship between my batches and the presence of a SNP.
pval = list()
for(i in 1:(length(colnames(GT_tr))-1)){

  form <- as.formula(paste(colnames(GT_tr)[i]," ~ batch", sep=""))
  form
  pval[colnames(GT_tr)[i]] <- summary(lm(form, as.data.frame(GT_tr)))$coefficients[2,4]
}
pval <- t(as.matrix(unlist(pval)))
rownames(pval) <- "corr_pval"

# add the pval data to the GT_tr matrix as the first row
GT_tr <- GT_tr[,1:(ncol(GT_tr)-1)]
GT_tr_p <- rbind(pval,GT_tr)

saveRDS(GT_tr_p, file=paste(fl,".CorrPvals.RDS",sep=""))

# find all the SNPs with a significant p-value <=0.01
# Significant p-value means there is a high correlation, which we don't want.
to_keep <- !(GT_tr_p["corr_pval",] <= 0.05)

# remove all those with a significant correlation
GT_tr_keep <- GT_tr_p[,to_keep]

# get a list of SNP positions from the column names
positions_to_keep <- colnames(GT_tr_keep)
pos_keep <- lapply(positions_to_keep, strsplit, "\\.")
pos_keep <- lapply(pos_keep,"[[",1)
chr_keep <- unlist(lapply(strsplit(unlist(lapply(pos_keep,"[",1)),"X"),"[",2))
coord_keep <- unlist(lapply(pos_keep,"[",2))
SNPs_to_keep <- data.frame(chr = chr_keep, pos = coord_keep)

write.table(file =paste(fl,".LowBatchCorrSNPs.positions.txt",sep=""), SNPs_to_keep, quote=F, row.names=F,sep="\t")

#system(mutt -s “VCF Batch Removal Completed” alolex@vcu.edu < Routput.txt)

