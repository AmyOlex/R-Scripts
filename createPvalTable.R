## Amy Olex
## 12/10/14
## This script creates the pvalue tables that acompany the RMarkdown documents with Takabe's survival analyses.

# Load needed libraries
library(plyr)
library(survival)
library(gridExtra)

createPvalFile <- function(INFILE, CODE, EXP_FILE, PVAL_ADJUST, alias, STEP_SIZE=0.1, UNIQUE_ID="pvals"){

  #set global varaibles
  #STEP_SIZE=.1
  #INFILE="tcga_DFS_091914.txt"
  #CODE=c("free"=0, "recurrence"=1)
  
  #INFILE="tcga_OS_092214.txt"
  #CODE=c("alive"=0, "dead"=1)
  
  
  # Load data from file and add a type column
  tcga <- read.delim(INFILE, header=TRUE, sep="\t", row.names=1)
  
  # Format the tcga object and create the survival object
  tcga$type = revalue(tcga$event, CODE)

  
  ## Load the gene expression data for tumor samples
  load(EXP_FILE)
  
  ndata <- qnorm_data
  
  # transpose the vsd_data matrix
  ndata_t = as.data.frame(t(ndata))
  
  #rename the genes to remove the pipe | in the name for easier parsing.
  newnames = gsub("\\|", ".", colnames(ndata_t))
  colnames(ndata_t) = newnames
  
  #Extract out the genes of interest from vsd_data_t and transpose it again at the same time so the genes are in the rows again
  ndata_gene = as.data.frame(t(ndata_t[,row.names(alias)]))
  
  #Now merge this with the Alias object and get the alias names as the row names
  ndata_gene2 = merge(alias, ndata_gene, by=0, all=TRUE)
  row.names(ndata_gene2) = ndata_gene2$alias
  ndata_gene2 = ndata_gene2[,-(1:2)]
  
  ## Transpose the matrix so genes are now the columns
  gene_exp <- as.data.frame(t(ndata_gene2))
  
  ## Remove any rows that are not solid "-01" tumors because they will introduce duplicate row names.
  ## First identify any entries that are not solid tumors with a -01 at end.
  tmp = substr(row.names(gene_exp), start=14, stop=15)
  cropped_gene_exp = gene_exp[which(tmp=="01"),]
  
  
  ## Re calculate the means
  means <- lapply(cropped_gene_exp, mean)
  
  ## Crop the row names so I can merge data frames
  row.names(cropped_gene_exp) = substr(row.names(cropped_gene_exp), start=1, stop=12)
  
  ## Create the high/low discretization using the mean expression for a gene across all tumor samples (N=273)
  high_low <- as.data.frame(ifelse(cropped_gene_exp <= means, "low", "high"))
  
  ## Merge the high/low values with the survival data, and only keep those that are in both sets (N=206)
  tcga_surv_genes = merge(tcga, high_low, by=0, all=FALSE)
  row.names(tcga_surv_genes) <- tcga_surv_genes$Row.names
  tcga_surv_genes <- tcga_surv_genes[,-1]
  
  print(dim(tcga_surv_genes))
  
  # Create the Surv() object using time and type, which is really alive or dead status.
  tcga_surv <- Surv(as.numeric(tcga$time), as.numeric(tcga$type))
  
  ### Now we are ready to plot stuff
  
  ## Get the transpose of the vsd data for the density plots
  n_t <- as.data.frame(t(ndata_gene2))
  pval_table_all <- data.frame()
  ## make a KM Plot for each Gene of Interest and calculate statistics
  for(gene in names(cropped_gene_exp)){
    # Find the min and max values for this gene
    # create a vector of values to loop over
    # Loop over values
    #     create high_low matrix
    #     merge high/low with tcga_surv_genes
    #     create survival object
    #     calculate pvalues and test for significance
    #     if significant, plot survival and distribution with data line
    #     if not significant move on to next value
    
    # create vector of 10 values between min and max
    vals <- seq(range(cropped_gene_exp[gene])[1],range(cropped_gene_exp[gene])[2], by=STEP_SIZE)
    pval_table <- data.frame(row.names=seq(length(vals)))
    i=0
    # Loop over vals
    for(v in vals){
      i=i+1
      
      HL <- as.data.frame(ifelse(cropped_gene_exp[gene] <= v, "low", "high"))
      
      ## Merge the high/low values with the survival data, and only keep those that are in both sets (N=206)
      tcga_surv_genes2 = merge(tcga, HL, by=0, all=FALSE)
      row.names(tcga_surv_genes2) <- tcga_surv_genes2$Row.names
      tcga_surv_genes2 <- tcga_surv_genes2[,-1]
      
      # Test to see if we have multiple groups
      if((length(count(tcga_surv_genes2[gene][[1]])$freq)<=1)){
        pval_table$gene[i] <- gene
        pval_table$cutoff[i] <- round(v,3)
        pval_table$pval[i] <- "NA"
        pval_table$nHigh[i] <- "NA"
        pval_table$nLow[i] <- "NA"
        pval_table$pval_adj[i] <- "NA"
        next
      }
      ## Create the Surv() object using time and type, which is really alive or dead status.
      tcga_surv <- Surv(as.numeric(tcga_surv_genes2$time), as.numeric(tcga_surv_genes2$type))
      
      ## Calculate the pvalue and test for significance at p=0.01
      sdf <- survdiff(tcga_surv~tcga_surv_genes2[gene][[1]])
      p.val_adjust <- p.adjust( (1 - pchisq(sdf$chisq, length(sdf$n) - 1)), method=PVAL_ADJUST, n=dim(alias)[1] )
      p.val <- (1 - pchisq(sdf$chisq, length(sdf$n) - 1))
      highN <- count(tcga_surv_genes2[,gene])$freq[1]
      lowN <- count(tcga_surv_genes2[,gene])$freq[2]
      
      pval_table$gene[i] <- gene
      pval_table$cutoff[i] <- round(v,3)
      pval_table$pval[i] <- round(p.val,5)
      pval_table$nHigh[i] <- highN
      pval_table$nLow[i] <- lowN
      pval_table$pval_adj[i] <- round(p.val_adjust,5)
    } #end for v in vals
    
    pval_table_all <- do.call(rbind, list(pval_table_all,pval_table))
  
  } #end for each gene
  
  filename=paste(Sys.Date(), "_survival_pvalTables_stepSize", as.character(STEP_SIZE), "_", UNIQUE_ID, "_", INFILE, sep="")
  
  # Save pval table to a file
  write.table(pval_table_all, file=filename, quote=FALSE, sep="\t", row.names=FALSE)
}
