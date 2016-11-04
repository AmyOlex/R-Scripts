# Amy Olex @ CCTR
# 6/18/14
# This function was written to help with TCGA file importing.
# It takes in the column headers, which should be bar codes with a -01 or -11 
# extension to indicate tumor or normal tissue, respectively, 
# and turns it into a factor.

createTissueFactor <- function(sample.names){
  newFactor <- t(apply(X=as.array(sample.names), MARGIN=1, FUN=gsub, pattern="TCGA-..-....-(..)", replacement="\\1", x=sample.names)[,1])
  newFactor[newFactor=="01"] <- "tumor"
  newFactor[newFactor=="11"] <- "normal"
  return(as.factor(newFactor))
}