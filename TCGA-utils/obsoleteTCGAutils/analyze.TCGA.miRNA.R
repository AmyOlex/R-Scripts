#####
# Amy Olex @ CCTR
# 6/10/14 Working file
# Analyzing TCGA miRNA HNSCC data.

#libraries
library(DESeq)
library(plyr)
library(gplots)
library(RColorBrewer)

#Commands run to get input data:
hpv.tumor<-import.TCGA.hpv_status("/Users/alolex/Desktop/CCTR/Data/06-2014_IainMorgan_TCGA/061014_miRexp_HNSCC_tumor_allBatches_wClinical/Clinical/Biotab/nationwidechildrens.org_auxiliary_hnsc.txt")
hpv.normal<-import.TCGA.hpv_status("/Users/alolex/Desktop/CCTR/Data/06-2014_IainMorgan_TCGA/061014_miRexp_HNSCC_norm_allBatches_wClinical/Clinical/Biotab/nationwidechildrens.org_auxiliary_hnsc.txt")
tumor.rdcounts <- import.TCGAmiRNA("/Users/alolex/Desktop/CCTR/Data/06-2014_IainMorgan_TCGA/061014_miRexp_HNSCC_tumor_allBatches_wClinical/miRNASeq/BCGSC__IlluminaHiSeq_miRNASeq/Level_3", "TCGA_miRNA_tumor_allBatch_wClinical_readcounts.txt", "read_count")
normal.rdcounts <- import.TCGAmiRNA("/Users/alolex/Desktop/CCTR/Data/06-2014_IainMorgan_TCGA/061014_miRexp_HNSCC_norm_allBatches_wClinical/miRNASeq/BCGSC__IlluminaHiSeq_miRNASeq/Level_3", "TCGA_miRNA_normal_allBatch_wClinical_readcounts.txt", "read_count")

# Crop the column headers in the read counts data so it is just the bar code
names(normal.rdcounts)<-apply(X=as.array(names(normal.rdcounts)), MARGIN=1, FUN=substr, start=1, stop=12)
names(tumor.rdcounts)<-apply(X=as.array(names(tumor.rdcounts)), MARGIN=1, FUN=substr, start=1, stop=12)

# Reorder the hpv status dataframe to be in the same order as the columns of the recounts
new.hpv.normal <- reorderRows.data.frame(hpv.normal, names(normal.rdcounts))
# Cannot run the below command because of duplicate column headers after cropping.  
# Need to fix this and expand the factor vector as well to compensate.
#new.hpv.tumor <- reorderRows.data.frame(hpv.tumor, names(tumor.rdcounts))

# Remove the normal columns of read counts where the HPV status is Positive.
new.normal.rdcounts <- normal.rdcounts[,new.hpv.normal$hpv_status=="Negative"]

######
# Now we have 2 data sets of read counts (new.norm.rdcounts and tumor.rdcounts). We will leave the tumor read count data
# alone for this first analysis as we are doing tumor vs normal/HPV-.
#####
# Combine the tumor and normal data, and create a factor vector for DESeq analysis.
combined.ds <- merge(new.normal.rdcounts, tumor.rdcounts, by=0, all=TRUE)
row.names(combined.ds) <- t(combined.ds[1])
combined.ds <- combined.ds[,!(names(combined.ds) %in% names(combined.ds)[1])]

conditions.ds <- factor(c(rep("normal", ncol(new.normal.rdcounts)), rep("tumor", ncol(tumor.rdcounts))))

# Now run through the DESeq Pipeline
# Create a countDataSet
design <- data.frame(row.names = names(combined.ds), condition = conditions.ds, libType = rep("single-end", ncol(combined.ds)))
cds <- newCountDataSet(combined.ds, design)
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

plotDispEsts( cds )
# Filter data on overall read counts and remove the lowest 40% (from section 5 of the DESeq vingnette)
rs = rowSums ( counts ( cds ))
theta = 0.4
use = (rs > quantile(rs, probs=theta))
cdsFilt = cds[ use, ]

#return(list(t.hpvStatus=hpv.tumor, n.hpvStatus=hpv.normal, t.rdcounts=tumor.rdcounts, n.rdcounts=new.norm.rdcounts, combined.rdcount=combined.ds, cond.fact=conditions.ds, cds.deseq=cds, cdsFilt.deseq=cdsFilt))


# There are 38 normal samples and 476 tumor samples
count(conditions(cdsFilt))
    
# Trying to run the expression calculations on all 514 columns crashed my computer.  Even half that crashed it.
# I can safely do about 30 columns total.  
# I select a random 10 normal samples and 20 tumor samples and run the expression analysis 5 different times.
# I then plot each and will figure out a way to summarize the results.

#res = nbinomTest( cdsFilt, "normal", "tumor" )
#plotMA(res)
#hist(res$pval, breaks=100, col="skyblue", border="slateblue", main="")
#resSig = res[ res$padj < 0.1, ]

resL <- list(0)
par(mfrow=c(5,2), mar=c(3,3,1,1))
for(i in 1:5) {
  resL[[i]] <- nbinomTest( cdsFilt[,c(sample(1:38, 10, replace=F), sample(39:514, 20, replace=F))], "normal", "tumor" )
  plotMA(resL[[i]], main=paste("cdsFilt Run ", i), cex = ifelse(resL[[i]]$padj>=0.1, .5, 1), col = ifelse(resL[[i]]$padj>=0.1, "gray32", "red3"))
  hist(resL[[i]]$pval, breaks=100, col="skyblue", border="slateblue", main=paste("cdsFilt Run ", i))
  print(paste("Completed ", i))
}

# Now I want to filter on the significant miRNAs and then look at the differences in the lists.
# Using the plyr function llply I am able to do this easily!
resSig <- llply(resL, subset, padj<0.1)

# Use venn to visualize overlap and to get the counts of intersection genes.
# I cannot find a pre-made function to intersect all the sets and return a list of data frames with just the significant genes.
# I will have to write this when it becomes necessary using a loop over Reduce(intersect, list(...list of character vectors here...))
venn(llply(resSig, subset, select = "id"))

# I can extract the fold-change and pval columns and concatenate them into a single data frame for clustering however of just the significant miRNAs.
# This process is similar to that of the import.TCGAmiRNA.R function I wrote.
sigFoldChanges.df <- resSig[[1]][which(names(resSig[[1]]) %in% c("id", "log2FoldChange", "padj"))]
sig.names <- c("id", apply(X=as.array(names(sigFoldChanges.df)[2:3]), MARGIN=1, FUN=paste, "_run1", sep=""))
for(df.idx in 2:length(resSig)){
  
  # Now extract the columns to parse and merge that with the sigFoldChanges 
  tmp.sigData <- resSig[[df.idx]][which(names(resSig[[df.idx]]) %in% c("id", "log2FoldChange", "padj"))]
  # Now merge this tmp data with the out.df data
  sigFoldChanges.df <- merge(sigFoldChanges.df, tmp.sigData, by = "id", all=TRUE)
  # Now update the column labels
  sig.names <- c(sig.names, apply(X=as.array(names(tmp.sigData))[2:3], MARGIN=1, FUN=paste, paste("_run",df.idx,sep=""), sep=""))
}

names(sigFoldChanges.df) <- sig.names

# OR, to get just the fold changes
sigFoldChanges.df <- resSig[[1]][which(names(resSig[[1]]) %in% c("id", "log2FoldChange"))]
for(df.idx in 2:length(resSig)){   
   # Now extract the columns to parse and merge that with the sigFoldChanges 
   tmp.sigData <- resSig[[df.idx]][which(names(resSig[[df.idx]]) %in% c("id", "log2FoldChange"))]
   # Now merge this tmp data with the out.df data
   sigFoldChanges.df <- merge(sigFoldChanges.df, tmp.sigData, by = "id", all=TRUE)
}
names(sigFoldChanges.df) <- c("id","run1","run2","run3","run4","run5")
#Make the ID column the row names
row.names(sigFoldChanges.df) <- t(sigFoldChanges.df[1])
sigFoldChanges.df <- sigFoldChanges.df[,!(names(sigFoldChanges.df) %in% names(sigFoldChanges.df)[1])]

#Now use heatmap.2 to cluster and get a heatmap.
hmcol=colorRampPalette(brewer.pal(9,"RdBu"))(100)
hmcol <- rev(hmcol) #reverses the color vector so red in large values and blue is small values.
# Had to alter the data frame to remove NA and Inf and change to zero for clustering.
# This is just for THIS example, I should NOT do this all the time!  Need to figure out how to get the intersection
# of the data frames.
sigFoldChanges.df[is.na(sigFoldChanges.df)] <- 0
sigFoldChanges.df[sigFoldChanges.df == "-Inf"] <- 0
heatmap.2(as.matrix(sigFoldChanges.df[,]), col=hmcol, trace="none", margin=c(10,6), na.color="black")

#Attempt 2 at extracting only the intersection set fo miRNAs.
# Strategy is to merge everything, replace any -Inf or +Inf with NA, then remove all these rows from the dataframe.
sigFoldChanges.df[sigFoldChanges.df == "-Inf"] <- NA  #I will need to write a better more general function than this for larger data sets.
overlap <- sigFoldChanges.df[apply(X=(!is.na(sigFoldChanges.df)), MARGIN=1, FUN=all),] #This extracts out all rows that have no NA values.



#Now do the heatmap
heatmap.2(as.matrix(sigFoldChanges.df[,2:ncol(overlap)]), col=hmcol, trace="none", margin=c(10,6))


