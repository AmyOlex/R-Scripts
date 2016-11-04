# Amy Olex @ VCU CCTR
# 6/3/14
# Purpose: creates a line version of a histogram to display the 
# distribution of values in each column (list element) in the data.frame object.
# Inputs:
#     data.df: a data frame containing the data I want to create a distribution plot for.  
#              If not already a data frame, it tries to convert the data to a data frame.
#     title: The title to put on the plot
#     file.path: The full path and file name to write minus the extension.

createDistFigure <- function(data.df, title, file.path){
  
  if(!is.data.frame(data.df)){
    data.df <- as.data.frame(data.df)
  }
  
  #create Bins vector from max and min in data.df
  bins <- seq(floor(min(data.df)), ceiling(max(data.df)), by=.1)
  #coerce data.df to a list (each column in the data will be a list element, and then apply the hist function 
  #to each element to get the distribution/counts data per bin.
  hist.data <- lapply(as.list(data.df), hist, breaks = bins, plot=FALSE)
  #now extract the counts vectors for each bin.
  counts.data <- getCounts(hist.data)
  #open figure file, create plot, close stream to figure file
  #note for plotting we use the $mids vector for the X axis.
  #pdf(paste(file.path, ".pdf"))
  matplot(hist.data[[1]]$mids, counts.data, type="l",lwd=2, lty=1:8, col=1:8, main=title, xlab="Bins", ylab="Counts")
  #legend("topright", legend = names(data.df), lty=1:8, lwd=2, col=1:8, bty="n", ncol=2, title="Legend")
  #dev.off()
}

#Helper function for the above function.
getCounts <- function(histList){
  tmp = numeric()
  for(i in histList){
    tmp = cbind(tmp, i$counts) 
  }
  return(tmp)
}