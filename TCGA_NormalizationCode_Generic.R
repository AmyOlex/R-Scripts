########## Generic Normalization Code for TCGA Data Analysis ###########
## Amy Olex ##
## Normalization of various subsets of TCGA RNASeq data.
## So I only have to do it once!
## Each block of normalization has an associated .RData file that needs to be saved.
## These data files can be imported by markdown docs etc for faster generation.
##
## Directions for use:
##    - Import all TCGA Data for a selected cancer into a tcgaRC object and save this to an RData file.
##    - Save this file into a "normalization" subdirectory of the main working directory for this cancer type.
##    - For example: ~/Desktop/CCTR/Data/KazuakiTakabe_S1P_0714_TCGA/TCGA Analyses/RNAseq_COAD_all_082714/normalization/coad_082714_geneExp.tcgaRC.RData
##    - This will be the file needed for this script to run.  
##    - some of the variables are generalized, but some are not.  Scan through file to change these variables for each dataset.

####################################################################
########### Common steps for all normalizations #############
####################################################################
######### Define Global Variables #########
working_dir <- "/Users/alolex/Desktop/CCTR/Data/KazuakiTakabe_S1P_0714_TCGA/TCGA Analyses/RNAseq_PAAD_all_111914/normalization"
main_file <- "paad_111914_tcgaRC.RData"
fact_file <- "paad_111914_clinical_facts.txt"
out_file_prefix <- "paad_111914"

######### Set up my environment ###########
library("DESeq2")
library("limma")
library("TeachingDemos")
library("scales")
library("plyr")
source('/Users/alolex/Desktop/CCTR/SourceCode/R-Scripts/loadTCGA.R')
loadTCGA()
setwd(working_dir)

############ Load Data #############
# Get the RNASeq Gene Expression TCGA object
load(main_file)

####CHANGE!
############ Load Additional Factors (Must have at least one other that is not Tumor dependent) ###########
added_facts <- createFactorFromFile(tab_file=fact_file, sample_order=paad_111914$getDataSamples())
for(i in names(added_facts)){
  paad_111914$setCustomFactor(new_fact=as.factor(added_facts[,i]), new_fact_name=i)
}

# Round to nearest Integer for DESeq2 normalization
data = paad_111914$getData()  ###CHANGE
rounded_data1 = round(data, digits=0)

###CHANGE the below!  Generally I'll just have tumor/normal
#coldata1 <- data.frame(tumor_status=coad_082714$getCustomFactor("tissue2")[[1]], 
#                      DFS=coad_082714$getCustomFactor("DFS")[[1]],
#                      OS=coad_082714$getCustomFactor("OS")[[1]],
#                      anatomic_site=coad_082714$getCustomFactor("anatomic_site")[[1]],
#                      histological_type=coad_082714$getCustomFactor("histological_type")[[1]],
#                      row.names=coad_082714$getDataSamples())
coldata1 <- data.frame(tumor_status=paad_111914$getTissueType(), 
                         tissue_source_site=paad_111914$getCustomFactor("tissue_source_site"),
                         histological_type=paad_111914$getCustomFactor("histological_type"),
                         anatomic_neoplasm_subdivision=paad_111914$getCustomFactor("anatomic_neoplasm_subdivision"),
                         gender=paad_111914$getCustomFactor("gender"),
                         race=paad_111914$getCustomFactor("race"),
                         pathologic_stage=paad_111914$getCustomFactor("pathologic_stage"),
                         person_neoplasm_cancer_status=paad_111914$getCustomFactor("person_neoplasm_cancer_status"),
                         row.names=paad_111914$getDataSamples())

####################################################################
########### Normalization over all patient samples #############
####################################################################
#### DESeq2 Normalizations ####
# Now I can create the DESeq2 object and run the normalization
rounded_data <- rounded_data1
coldata <- coldata1
dds <- DESeqDataSetFromMatrix(countData = rounded_data, colData = coldata, design = ~ tumor_status)
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
vsd_data <- as.data.frame(assay(vsd))

#### Other normalizations ####
logboxcox_data <- bct(as.matrix(rounded_data+1),0)
logshift_data <- as.data.frame(log2(rounded_data + 1))
qnorm_data <- normalizeQuantiles(log2(rounded_data+1))
voom_cyc_data <- voom(rounded_data, normalize.method="cyclicloess")
voom_scale_data <- voom(rounded_data, normalize.method="scale")
voom_quant_data <- voom(rounded_data, normalize.method="quantile")
voom_none_data <- voom(rounded_data, normalize.method="none")

######### Saving objects to a file ############
## Create Out file name
out_file <- paste(out_file_prefix, "_allSamplesNormalized_n", dim(rounded_data)[2], ".RData", sep="")
save(list=c("data", "rounded_data", "coldata", "dds", "vsd", "vsd_data", "logboxcox_data", "logshift_data", "qnorm_data", "voom_cyc_data", "voom_scale_data", "voom_none_data", "voom_quant_data"), file = out_file)




#####################################################################
########### Normalization over Tumor Samples only ###############
#####################################################################
rounded_data = as.data.frame(rounded_data1[coldata1$tumor_status=="tumor"])
coldata=coldata1[coldata1$tumor_status=="tumor",]
dds <- DESeqDataSetFromMatrix(countData = rounded_data, colData = coldata, design = ~ gender)
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
vsd_data <- as.data.frame(assay(vsd))

#### Other normalizations ####
logboxcox_data <- bct(as.matrix(rounded_data+1),0)
logshift_data <- as.data.frame(log2(rounded_data + 1))
qnorm_data <- normalizeQuantiles(log2(rounded_data+1))
voom_cyc_data <- voom(rounded_data, normalize.method="cyclicloess")
voom_scale_data <- voom(rounded_data, normalize.method="scale")
voom_quant_data <- voom(rounded_data, normalize.method="quantile")
voom_none_data <- voom(rounded_data, normalize.method="none")

######### Saving objects to a file ############
out_file <- paste(out_file_prefix, "_tumorOnlySamplesNormalized_n", dim(rounded_data)[2], ".RData", sep="")
save(list=c("data", "rounded_data", "coldata", "dds", "vsd", "vsd_data", "logboxcox_data", "logshift_data", "qnorm_data", "voom_cyc_data", "voom_scale_data", "voom_none_data", "voom_quant_data"), file = out_file)





