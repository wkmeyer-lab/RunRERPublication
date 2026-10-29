.libPaths("/share/ceph/wym219group/shared/libraries/R4") #add path to custom libraries to searched locations
#Sys.setenv('R_MAX_VSIZE'=20000000000)
library(RERconverge) #load RERconverge package
library(RERconverge)
library("tools")
library(stringr)
source("Src/Reu/cmdArgImport.R")
source("Src/Reu/convertLogiToNumeric.R")
source("Src/Reu/permPValCorReport.R")

#---- USAGE -----
#used to calculate p values on combined permulations files made by CombinePermualtions.R. 

# ARUGMENTS: 
#If an argument contains a '(' it is evaluated as code.
# 'r="filePrefix"'              This is the prefix attached to all files; a required argument. 
# 't = <s OR f OR p>'           This sets which permulation filetype to look for. s is for slow, f is for fast, and p is for pruned-fast
               

#-------
#Debug setup defaults
#permulationNumberValue = 3
#Testing args:
args = c('r=demoinsectivory', 'n=3', 'e=F', 't=p')
args = c('r=allInsectivory', 'e=F', 's=1', 't=s', 'i=1')
args = c('r=carnvHerbs', "t=p")
args = c('r=CVHRemake', "t=p")
#------


# --- Import prefix ----
args = commandArgs(trailingOnly = TRUE)
paste(args)
message(args)

#File Prefix
if(!is.na(cmdArgImport('r'))){
  filePrefix = cmdArgImport('r')
}else{
  stop("THIS IS AN ISSUE MESSAGE; SPECIFY FILE PREFIX")
}

#------------

#------ Make Output directory -----

#Make output directory if it does not exist
if(!dir.exists("Output")){
  dir.create("Output")
}
#Make a specific subdirectory if it does not exist 
outputFolderNameNoSlash = paste("Output/",filePrefix, sep = "")
#create that directory if it does not exist
if(!dir.exists(outputFolderNameNoSlash)){
  dir.create(outputFolderNameNoSlash)
}
outputFolderName = paste("Output/",filePrefix,"/", sep = "")

#----------

#----- Default values -------

runInstanceValue = NULL

#-------


# -- Import the instance number of the script --- 
if(!is.na(cmdArgImport('i'))){
  runInstanceValue = cmdArgImport('i')
}else{
  message("This script does not have a run instance value")
}

# -- Import which filetype to use  --- 
if(!is.na(cmdArgImport('t'))){
  fileType = cmdArgImport('t')
  fileType = as.character(fileType)
}else{
  message("No fileType specified, defaulting to slow (PermulationsData)")
}

# ------- append the files ------- 

# ---- Get filetype string ---
#Do the first combination:
if(fileType == 's'){
  fileTypeString = ''
}else if(fileType == 'f'){
  fileTypeString = "UnprunedFast"
}else if(fileType == 'p'){
  fileTypeString = "PrunedFast"
}else{
  stop( "THIS IS AN ISSUE MESSAGE, IMPROPER FILETYPE ARGUMENT")
}

#Figure out the file names
filesInFolder = list.files(outputFolderName)
paste(filesInFolder)
fileSetStep1 = grep(paste("Combined", fileTypeString, sep =""), filesInFolder, value = TRUE)
paste(fileSetStep1)
fileSetStep2 = grep("PermulationsPValue", fileSetStep1, value = TRUE)
paste(fileSetStep2)
OtherAppendOutputFiles = grep("Appended", fileSetStep2)
OtherAllOutputFiles= grep("All", fileSetStep2)
OtherAppendOutputFiles = append(OtherAppendOutputFiles, OtherAllOutputFiles)
paste(OtherAppendOutputFiles)
if(is.integer(OtherAppendOutputFiles) & length(OtherAppendOutputFiles) != 0){ 
  fileSetStep3 = fileSetStep2[-OtherAppendOutputFiles]
}else{fileSetStep3 = fileSetStep2}
message(fileSetStep3)



#sort the files correctly 
fileSetNamesSort = fileSetStep3
fileSetNamesSort = sapply(fileSetNamesSort, str_sub, end = -23)
fileSetNamesSort = sapply(fileSetNamesSort, gsub, pattern = ".*-", replacement = "")
fileSetNamesSort = fileSetNamesSort[order(as.integer(fileSetNamesSort))]

fileSetStep4 = names(fileSetNamesSort)
message(fileSetStep4)

#append the files
message("Appending files:")
appenedPermPValues = NULL
#import all of the files
for(i in 1:length(fileSetStep4)){
  currentFile = readRDS(paste(outputFolderName, fileSetStep4[i], sep= ""))
  appenedPermPValues = append(appenedPermPValues, currentFile)
  message(paste(outputFolderName, fileSetStep4[i], sep= ""))
}

while((tail(names(appenedPermPValues), n=1))==""){
  appenedPermPValues = appenedPermPValues[-length(appenedPermPValues)]
}



appendedPerPValuesFilename = paste(outputFolderName, filePrefix, "Combined", fileTypeString, "Appended", "PermulationsPValue", runInstanceValue, ".rds", sep= "")
saveRDS(appenedPermPValues, file = appendedPerPValuesFilename)

