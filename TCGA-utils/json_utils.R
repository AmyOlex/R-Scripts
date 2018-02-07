# 3/7/17
# Amy Olex
# Methods to manipulate JSON Files from GDC

library(jsonlite)
library(rjson)

convertJSON <- function(json_file, outfile){
  # read in file
  #my_json <- readLines(json_file)
  # convert to data frame
  json.df <- jsonlite::fromJSON(json_file)
  json.df$caseID <- unlist(json.df$cases)
  json.df$cases <- NULL
  # write to outfile
  write.table(json.df, file=outfile, quote=FALSE, sep="\t", row.names=FALSE)
}


JSON.to.df <- function(json_file){
  # read in file
  #my_json <- readLines(json_file)
  # convert to data frame
  json.df <- jsonlite::fromJSON(json_file)
  json.df$caseID <- unlist(json.df$cases)
  json.df$cases <- NULL
  # write to outfile
  return(json.df)
}

