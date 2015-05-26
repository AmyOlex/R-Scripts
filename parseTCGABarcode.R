parseTCGABarcode <- function(sample_names, detailed=FALSE){
  new_factor <- t(apply(X=as.array(sample_names), MARGIN=1, FUN=gsub, pattern="TCGA-..-....-(..)", replacement="\\1", x=sample_names)[,1])
  
  if(detailed){
    new_factor[new_factor=="01"] <- "TP"
    new_factor[new_factor=="02"] <- "TR"
    new_factor[new_factor=="03"] <- "TB"
    new_factor[new_factor=="04"] <- "TRBM"
    new_factor[new_factor=="05"] <- "TAP"
    new_factor[new_factor=="06"] <- "TM"
    new_factor[new_factor=="07"] <- "TAM"
    new_factor[new_factor=="08"] <- "THOC"
    new_factor[new_factor=="09"] <- "TBM"
    new_factor[new_factor=="10"] <- "NB"
    new_factor[new_factor=="11"] <- "NT"
    new_factor[new_factor=="12"] <- "NBC"
    new_factor[new_factor=="13"] <- "NEBV"
    new_factor[new_factor=="14"] <- "NBM"
    new_factor[new_factor=="20"] <- "CELLC"
    new_factor[new_factor=="40"] <- "TRB"
    new_factor[new_factor=="50"] <- "CELL"
    new_factor[new_factor=="60"] <- "XP"
    new_factor[new_factor=="61"] <- "XCL"
    return(as.factor(new_factor))
  }
  else{
    new_factor[new_factor=="01"] <- "tumor"
    new_factor[new_factor=="02"] <- "tumor"
    new_factor[new_factor=="03"] <- "tumor"
    new_factor[new_factor=="04"] <- "tumor"
    new_factor[new_factor=="05"] <- "tumor"
    new_factor[new_factor=="06"] <- "tumor"
    new_factor[new_factor=="07"] <- "tumor"
    new_factor[new_factor=="08"] <- "tumor"
    new_factor[new_factor=="09"] <- "tumor"
    new_factor[new_factor=="10"] <- "normal"
    new_factor[new_factor=="11"] <- "normal"
    new_factor[new_factor=="12"] <- "normal"
    new_factor[new_factor=="13"] <- "normal"
    new_factor[new_factor=="14"] <- "normal"
    new_factor[new_factor=="20"] <- "control"
    new_factor[new_factor=="40"] <- "tumor"
    new_factor[new_factor=="50"] <- "cell line"
    new_factor[new_factor=="60"] <- "xenograft"
    new_factor[new_factor=="61"] <- "xenograft"
    return(as.factor(new_factor))
  }
}
