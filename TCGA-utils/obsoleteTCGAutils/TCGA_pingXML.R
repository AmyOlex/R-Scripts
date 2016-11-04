TCGA_pingXML <- function(file_dir = stop("Error: Must enter directory where files are stored."), file_type = "clinical"){
  
  # Get file names in directory and throw an error if the list comes back empty
  files <- list.files(path=file_dir, pattern=file_type)
  if(length(files)==0) stop(paste("Error: No files were found with \"", file_type, "\" as part of the file name in the directory under ", file_dir, sep=""))
  
  # import the first file to an XML object
  data <- xmlTreeParse(paste(file_dir,files[1],sep=""), useInternal=TRUE)
  
  # Get the root node
  top <- xmlRoot(data)
  # Print out the XML node names under the patient category.
  row.names(as.matrix(names(top[[2]])))
    
}