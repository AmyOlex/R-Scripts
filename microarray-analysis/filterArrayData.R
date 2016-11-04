# Amy Olex @ VCU CCTR
# 6/2/14
# Purpose: Filters microarray data by SLR and p-value using Amy Olex's method.
# Inputs: 
#   pvalCutoff: detection pval cutoff for filtering, usually set of Affy Marginal value of .06
#   SLRCutoff: Signal Log Ratio cutoff for gene expression levels, usually set to .5, which is equivalent to a 1.5-fold change.
#   numConditions: The number of conditions or time point per replicate.
#   title: unique name for this analysis
#   savePath: location to save result files in.
#

# Import needed libraries
library("ggplot2")

# Source other functions being used.
source('~/Desktop/CCTR/SourceCode/R/Graphics/createDistFigure.R')
source('~/Desktop/CCTR/SourceCode/R/Graphics/myImagePlot.R')

filterArrayData <- function(pvalCutoff, SLRCutoff, numConditions, title, dataFile, savePath)
{
  #Check to ensure input is all in the write format.
  if(!is.numeric(pvalCutoff) || !is.numeric(SLRCutoff) || !is.numeric(numConditions)){
    stop("Error in input: pvalCutoff, SLRCutoff and numConditions must be numeric.", call.=TRUE)
  }
  if(!is.character(title)){
    stop("Error in input: title must be a character string.", call.=TRUE)
  }
  if(!file.exists(dataFile)){
    stop(paste("Error: File does not exist: ", dataFile), call.=TRUE)
  }
  if(!file.exists(savePath)){
    stop(paste("Error: Path does not exist: ", savePath), call.=TRUE)
  }
  
  #Import data file into a data.frame
  arrayData <- read.table(dataFile, header = TRUE, sep="\t")
  
  #Set up variables and such
  pvalIdx = numConditions+1
  
  ####
  # Future improvement suggested by Aaron is to filter the columns by the names., e.g. all columns containing pval are pval columns.
  # This would make it a bit more flexable with the formatting as well as they would not have to be in the same order, AND one would not have to specify the number of experimental conditions.
  ####
  ####
  # P-value Filter
  ####
  
  #Filter the data on pvalue and save this information in a new data.frame object.
  sig.vals <- arrayData[pvalIdx:ncol(arrayData)] <= pvalCutoff ##returns a matrix of TRUE/FALSE values for the specified range in the data matrix indicating whether each individual pvalue has passed the threshold.
  sig.cols <- apply(sig.vals, MARGIN = 1, all) ##returns a vector of TRUE/FALSE values. apply() applying the all() function to each row of sig.vals (MARGIN=1) to ensure all values pass the threshold.
  
  PvalFilteredData <- arrayData[sig.cols, ] ##returns a data.frame object containing all columns of the filtered data.
  
  #Write the filtered data to a file.
  filename = paste(savePath, title, "_PvalFiltered_", nrow(PvalFilteredData), ".txt")
  write.table(PvalFilteredData, file = filename, append = FALSE, quote = FALSE, sep = "\t", eol = "\n")
  
  ####
  # SLR Filter
  ####
  # Take the PvalFilteredData frame and filter it for SLR significance.
  # First identify those values meeting either the positive or negative filter and create a boolean matrix.
  sig.pos.slrs <- PvalFilteredData[1:numConditions] >= SLRCutoff
  sig.neg.slrs <- PvalFilteredData[1:numConditions] <= -SLRCutoff
  
  # Second apply the any function to each row of the positive and negative matrix to identify genes that
  # met either criteria in at least one condition. These will be vectors of boolean values.
  sig.pos <- apply(sig.pos.slrs, MARGIN=1, any)
  sig.neg <- apply(sig.neg.slrs, MARGIN=1, any)
  
  # Third, concatenate the above two boolean vectors and use the apply function again to get one vector 
  # indicating if the gene in a given row passed either positive or negative filter for any condition.
  sig.gene.rows <- apply(cbind(sig.pos, sig.neg), MARGIN=1, any)
  
  # Finally, use this vector to identify the rows in the data frame that we want to keep.
  # This is another data.frame object that we now write to a file.
  SLRFilteredData <- PvalFilteredData[sig.gene.rows, 1:numConditions]
  
  #Write the filtered data to a file.
  filename <- paste(savePath, title, "_SLRFiltered_", nrow(SLRFilteredData), ".txt")
  write.table(SLRFilteredData, file = filename, append = FALSE, quote = FALSE, sep = "\t", eol = "\n")
  
  ##Generate a histogram of SLR filtered data.
  #bins <- seq(-5, 7, by=.1)
  #hist.data <- lapply(as.list(SLRFilteredData), hist, breaks = bins, plot=FALSE)
  #counts.data <- getCounts(hist.data)
  ##open figure file, create plot, close stream to figure file
  #pdf(paste(title, "_SLRFilteredData_Histogram.pdf"))
  #matplot(hist.data[[1]]$mids, counts.data, type="l",lwd=2, lty=1:8, col=1:8, main=paste(title, " SLR Filtered Data Distribution"), xlab="Bins", ylab="SLR Counts")
  #legend("topright", legend = names(myData), lty=1:8, lwd=2, col=1:8, bty="n", ncol=2, title="Legend")
  #dev.off()
  
  createDistFigure(data.df = SLRFilteredData, title = paste(title, "SLR Filtered Data Distribution"), file.path = paste(savePath, title, "_SLRFilteredData_Histogram") )
  
  ##Generate the heatmap of the array correlations
  pdf(paste(savePath, title, "_SLRFilteredData_Correlations.pdf"))
  myImagePlot(as.matrix(cor(SLRFilteredData)), title=paste(title, "SLR Filtered Data Correlation Heatmap"))
  dev.off()
  return(SLRFilteredData)
}

