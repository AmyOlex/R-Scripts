## Amy Olex
## 1/17/17
##
## Function to calculate the TPM values given the raw read counts and transcript lengths.

## Input is 2 dataframes.  
## data : DF that contains feature names in first column, raw read counts, and samples names as column headers.
## feature_lengths : DF that contains feature names (in same order as data) in first column, and the header "length".

calcTPM <- function(data, feature_length){
  # initilize RPK matrix
  RPK <- matrix(0, dim(data)[1], dim(data)[2])
  # Step 1: calculate first step in TPM calculation: count/feature_length
  for(i in 1:dim(data)[1]){
    for(j in 1:dim(data)[2]){
      RPK[i,j] <- data[i,j]/feature_length[i,1]
    }
  }
  
  # Step 2: get the sum of all columns and divide by 1,000,000 to get scaling factor
  scaling_factor <- colSums(RPK)/1000000
  
  # Step 3: Divide all RPK values by the scaling factor
  TPM <- as.data.frame(t(t(RPK)/scaling_factor))
  
  # Convert to Dataframe and add in row names and column names
  names(TPM) <- names(data)
  row.names(TPM) <- row.names(data)
  
  return TPM
}
