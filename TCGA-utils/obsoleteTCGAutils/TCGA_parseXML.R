# Amy Olex @ CCTR
# 7/23/2014
#
# This script is designed to parse out specific information from TCGA clinical files in XML format.
# Once familiar with the syntax I'll need to make this more generic so we can get any data we want.
#
# For now, to make it work with any file, we will have to look up the namespaces and include those in the input array.
#
# Sample code to run on Takabe data.
# >> TCGA_parseXML(file_dir="/Users/alolex/Desktop/CCTR/Data/KazuakiTakabe_S1P_0714_TCGA/Clinical_COAD_all/Clinical/XML/", ...
#                   out_file="/Users/alolex/Desktop/CCTR/Projects/KazuakiTakabe_S1P_0714/TCGA_Parsed_Clinical/<file_name>")

prognosis=c("shared:days_to_last_known_alive", "shared:days_to_death", "shared:days_to_last_followup", "shared:days_to_initial_pathologic_diagnosis",
            "shared:residual_tumor", "shared:primary_lymph_node_presentation_assessment", "shared:lymph_node_examined_count", "shared:number_of_lymphnodes_positive_by_he",
            "shared:number_of_lymphnodes_positive_by_ihc", "coad_read_shared:non_nodal_tumor_deposits", "shared:venous_invasion","shared:lymphatic_invasion",
            "coad_read_shared:perineural_invasion_present", "coad_read_shared:synchronous_colon_cancer_present", "shared:days_to_last_followup"
            "shared:days_to_death", "")



library("XML")

TCGA_parseXML <- function(file_dir = stop("Error: Must enter directory where files are stored."), 
                          file_type = "clinical", 
                          to_parse=c("shared:tumor_tissue_site", "shared:histological_type", "shared:gender", "shared:days_to_birth", "shared:days_to_last_followup", 
                                     "shared:age_at_initial_pathologic_diagnosis", "shared:year_of_initial_pathologic_diagnosis", "shared:person_neoplasm_cancer_status",
                                     "shared:weight", "shared:height", "/shared_stage:pathologic_M", "shared:anatomic_neoplasm_subdivision",
                                     "/cqcf:tumor_type", "/cqcf:diagnosis_subtype", "/nte:new_tumor_event_after_initial_treatment"), 
                          out_file="parsedXML"){
    
  # create list of files to parse (Each patient data is in a seperate file)
  #files = c("sample_input/Alive.xml", "sample_input/Dead.xml", "sample_input/Died.xml")
  
  # Get file names in directory and throw an error if the list comes back empty
  files <- list.files(path=file_dir, pattern=file_type)
  if(length(files)==0) stop(paste("Error: No files were found with \"", file_type, "\" as part of the file name in the directory under ", file_dir, sep=""))
  
  print(length(files))
  # create an empty data.frame
  
  df <- data.frame()
  
  # now loop through files
  for(x in files){
    # import to an XML object
    data <- xmlTreeParse(paste(file_dir,x,sep=""), useInternal=TRUE)
    
    # now get the disease code for the namespace
    code <- tolower(xpathSApply(data, "//admin:disease_code", xmlValue))
    
    # now parse out the information needed one by one
    # First get the bar code, because that will be the row name.
    xpath <- paste("//", code, ":patient/shared:bcr_patient_barcode", sep="")
    barcode <- xpathSApply(data, xpath, xmlValue)
    
    # Now get all the other information and add it to the data frame
    for(p in to_parse){
      xpath <- paste("//", code, ":patient/", p, sep="")
      value <- xpathSApply(data, xpath, xmlValue)
      if(length(value)==0){
        df[barcode, sub(".*:", "", p)] = "Field not found"
      }
        
      else if(length(value)==1){
        if(value==""){
          status <- xpathSApply(data, xpath, fun=xmlGetAttr, name="procurement_status")
          df[barcode, sub(".*:", "", p)] = status
        }
        else
          df[barcode, sub(".*:", "", p)] = value
      }
      else{
        tmp_val = ""
        for(v in value){
          tmp_val = paste(tmp_val, ",",v)
        }
        df[barcode, sub(".*:", "", p)] = tmp_val
      }
        
      
      
    }
    xpath <- paste("//", code, ":patient//shared_stage:pathologic_stage", sep="")
    df[barcode, "pathologic_stage"] = xpathSApply(data, xpath, xmlValue)
    xpath <- paste("//", code, ":patient/shared:vital_status", sep="")
    df[barcode, "initial_vital_status"] = xpathSApply(data, xpath, xmlValue)
    
    if(df[barcode, "initial_vital_status"]=="Dead")
      xpath <- paste("//", code, ":patient/shared:days_to_death[@procurement_status='Completed']", sep="")
    else
      xpath <- paste("//", code, ":follow_ups//shared:days_to_death[@procurement_status='Completed']", sep="")
    death = xpathSApply(data, xpath, xmlValue)
    if(length(death)==0){
      df[barcode, "days_to_death"] = "Field Not Found"
      df[barcode, "final_vital_status"] = "Field Not Found"
    }
    else{
      df[barcode, "days_to_death"] = death[1]
      df[barcode, "final_vital_status"] = "Dead"
    }
  }
  
  
  ########################################### Biospecimen #########################################
  # Now do something similar for the biospecimen files
  
  to_parse_bio=c("/shared:histological_type_other", "/cqcf:maximum_tumor_dimension")
  # Get file names in directory and throw an error if the list comes back empty
  files2 <- list.files(path=file_dir, pattern="biospecimen")
  if(length(files2)==0) stop(paste("Error: No files were found with \"", "biospecimen", "\" as part of the file name in the directory under ", file_dir, sep=""))
  
  print(length(files2))
  # create an empty data.frame
  
  #df <- data.frame()
  
  # now loop through files
  for(x in files2){
    # import to an XML object
    data <- xmlTreeParse(paste(file_dir,x,sep=""), useInternal=TRUE)
    
    # now get the disease code for the namespace
    code <- "bio"
    
    # now parse out the information needed one by one
    # First get the bar code, because that will be the row name.
    xpath <- paste("//", code, ":patient/shared:bcr_patient_barcode", sep="")
    barcode <- xpathSApply(data, xpath, xmlValue)
    
    # Now get all the other information and add it to the data frame
    for(p in to_parse_bio){
      xpath <- paste("//", code, ":patient/", p, sep="")
      value <- xpathSApply(data, xpath, xmlValue)
      if(length(value)==0){
        df[barcode, sub(".*:", "", p)] = "Field not found"
      }
      else if(length(value)==1){
        if(value==""){
          status <- xpathSApply(data, xpath, fun=xmlGetAttr, name="procurement_status")
          df[barcode, sub(".*:", "", p)] = status
        }
        else
          df[barcode, sub(".*:", "", p)] = value
      }
      else{
        tmp_val = ""
        for(v in value){
          print(paste(barcode, "--", xpath, "--", v))
          tmp_val = paste(tmp_val, ",",v)
        }
        df[barcode, sub(".*:", "", p)] = tmp_val
      }
      
      
      
    }
    
  }
  
  
  
  # Now write the dataframe to a file in the out directory provided, default is working directory
  write.table(df, file=paste(out_file, ".txt", sep=""), quote=FALSE, sep="\t")
  
  
} # end function