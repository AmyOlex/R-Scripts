loadTCGA <- function(computer="mac"){
  
  if(computer=="hershey"){
    source("~/src/R/TCGA/tcgaRC.R")
    source("~/src/R/TCGA/TCGA_condenseLvl3.R")
    source("~/src/R/TCGA/TCGA_importBiotab.R")
  }
  else{
    source("/Users/alolex/Desktop/CCTR/SourceCode/R/github_TCGA/tcgaRC.R")
    source("/Users/alolex/Desktop/CCTR/SourceCode/R/github_TCGA/TCGA_condenseLvl3.R")
    source("/Users/alolex/Desktop/CCTR/SourceCode/R/github_TCGA/TCGA_importBiotab.R")
  }
  
  
}