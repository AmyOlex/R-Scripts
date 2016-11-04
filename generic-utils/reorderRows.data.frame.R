#####
# Amy Olex @ CCTR
# 6/10/2014
#####
# Purpose: To reorder the rows of a data.frame using row.labels and a pre-defined ordereding specificed by an input vector.
#             
# Usage: >> new.df <- reorderRows.data.frame(dataframe, newOrderVector)
#
######

reorderRows.data.frame <- function(df, new.order){
  # extract row names from DF
  r.names <- row.names(df)
  
  # Compare to see if already in the same order
  if(all(r.names == new.order)){
    return(df)
  }
  else{
    # Get the reordered indicies
    reorder.idx <- match(new.order, r.names)
    
    # Test to see if the reorder worked
    if(!all(row.names(df)[reorder.idx] == new.order)){
      stop("Error: reorder failed")
    }
    
    # Reconstruct the DF with the new ordering
    reordered.df <- as.data.frame(df[reorder.idx,])
    row.names(reordered.df) <- row.names(df)[reorder.idx]
    names(reordered.df) <- names(df)
    
    # Return the new df
    return(reordered.df)
    
  }
  
  
}