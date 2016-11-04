# Amy Olex @ CCTR
# 6/18/14
# This function was written to help with TCGA file importing.
# When importing TCGA data, the column headers are set to be the full file names.
# The FILE_SAMPLE_MAP.txt file distributed with the TCGA meta data contains a map from 
# the file names to the unique bar codes for each sample.
# This function takes the file names from the column headers (input as a character vector x)
# and loops through each one individually, finds the index of the same entry in the first column of the sampleMap
# and inserts the associated bar code in the newX vector.  This way the column names are returned in the same order
# as they were input.  Finally, the newX vector elements are cropped to the first 15 characters.
# This includes the unique patient barcode plus 2 digits that specify whether the sample is from 
# tumor (-01) or normal (-11) tissue.

mapBarcodes <- function(sampleMap, x){
  map <- read.table(sampleMap, header=T)
  newX <- dim(x)
  for(i in 1:length(x)){
    newX[i] = as.character(map$barcode.s.[map$filename %in% x[i]])
  }
  newX <-apply(X=as.array(newX), MARGIN=1, FUN=substr, start=1, stop=15)
  return(newX)
}