#Library setup 
.libPaths("/share/ceph/wym219group/shared/libraries/R4") #add path to custom libraries to searched locations
library(RERconverge) #load RERconverge package
library(RERconverge)
library("tools")
source("Src/Reu/cmdArgImport.R")
source("Src/Reu/convertLogiToNumeric.R")
source("Src/Loc/permPValCorReport.R")


#---- USAGE -----
#used to combine permulations files made by RunPermulationsManual.R. 

#--IMPORTANT NOTE: This script takes exponentially longer the more files it is used to combine at once. 
#    It is strongly recommended that you ONLY COMBINE UP TO 200 permulation files in a single instance of the script, then switch to meta-combining those files. 


# ARUGMENTS: 
#If an argument contains a '(' it is evaluated as code.
# 'r="filePrefix"'              This is the prefix attached to all files; a required argument. 
# 'n=numberOfPermulations'      This is the number of permulation files the script will try to combine
# 'e=F' OR 'e=T'                This is if the permulations being combined are enriched or not. Accepts 'T', 'F', 'TRUE', 'FALSE', '0', and '1'. 
# 's=<number>'                  This is the permulation number to start at. Used for parrallelization. 
# 'i=<number>'                  This is used to generate unique filenames for each instance of the script. Typically fed in by for loop used to run script in parallel.
# 'c=F' OR 'c=T'                This is used to set if the script is being run to combine previous combinations. Called "metacombination". Used for parrallelization. 
# 't = <s OR f OR p>            This sets which permulation filetype to look for. s is for slow, f is for fast, and p is for pruned-fast
# 'p = <T or F>'                This value determines whether or not to calculate the pValues. [DEPRECATED]
#-------
#Debug setup defaults
#permulationNumberValue = 3
#Testing args:
args = c('r=demoinsectivory', 'n=3', 'e=F', 't=p')
args = c('r=allInsectivory','n=3', 'e=F', 's=1', 't=p')
args = c('r=carnvHerbs','n=3', 'e=F', 's=1', 't=p')
#------
timeStart = Sys.time()
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

permulationNumberValue = 100
enrichValue = F
startValue = 1
runInstanceValue = NULL
metacombineValue = FALSE
fileType = 's'
calulateValue = TRUE


#-------

# ---- Command args import ----
{
#Import number of permulations to combine
if(!is.na(cmdArgImport('n'))){
  permulationNumberValue = cmdArgImport('n')
  permulationNumberValue = as.numeric(permulationNumberValue)
}else{
  paste("Number of permulations not specified, using 100")
}

#Import if enriched or not 
if(!is.na(cmdArgImport('e'))){
  enrichValue = cmdArgImport('e')
  enrichValue = as.logical(enrichValue)
  if(is.na(enrichValue)){
    enrichValue = FALSE
    paste("Enrichment value not interpretable as logical. Did you remember to capitalize? Using FALSE.")
  }
}else{
  paste("enrichment not specified, using FLASE")
}

#Import the permulation number to start at
if(!is.na(cmdArgImport('s'))){
  startValue = cmdArgImport('s')
  startValue = as.numeric(startValue)
}else{
  paste("Start value not specified, using 1")
}



#Import if this is being run to combine combinations 
if(!is.na(cmdArgImport('c'))){
  metacombineValue = cmdArgImport('c')
  metacombineValue = as.logical(metacombineValue)
  if(is.na(metacombineValue)){
    metacombineValue = FALSE
    paste("Metacombination value not interpretable as logical. Did you remember to capitalize? Using FALSE.")
  }
}else{
  paste("Metacombination value not specified, using FALSE. If you aren't parrallelizing, don't worry about this.")
}

#Import which filetype to use  
if(!is.na(cmdArgImport('t'))){
  fileType = cmdArgImport('t')
  fileType = as.character(fileType)
}else{
  paste("No fileType specified, defaulting to slow (PermulationsData)")
}

#Import if this is being run to combine combinations
if(!is.na(cmdArgImport('p'))){
  calulateValue = cmdArgImport('p')
  calulateValue = as.logical(calulateValue)
  if(is.na(calulateValue)){
    calulateValue = TRUE
    paste("p-value calulation value not interpretable as logical. Did you remember to capitalize? Using FALSE.")
  }
}else{
  paste("p-value calulation value not specified, using TRUE. If you aren't parrallelizing, don't worry about this.")
}
}
# -- Combine permulations data files ----

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

if(metacombineValue == F){
    basePermulationsFilename = paste(outputFolderName, filePrefix, fileTypeString, "PermulationsData",  sep="")
}else{
  basePermulationsFilename = paste(outputFolderName, filePrefix, "Combined", fileTypeString, "PermulationsData",  sep="")
}

firstPermulationsFilename = paste(basePermulationsFilename, startValue, ".rds", sep="")
firstPermulationsData = readRDS(firstPermulationsFilename)
firstPermulationsData = convertLogiToNumericList(firstPermulationsData)

fpTime = Sys.time()
fpLoad = fpTime - timeStart
message("First permulation load time: ", fpLoad, attr(fpLoad, "units"))

secondPermulationsFilename = paste(basePermulationsFilename, (startValue+1), ".rds", sep="")
secondPermulationsData = readRDS(secondPermulationsFilename)
secondPermulationsData = convertLogiToNumericList(secondPermulationsData)

spTime = Sys.time()
spLoad = spTime - fpTime
message("Second permulation load time: ", spLoad, attr(spLoad, "units"))

combinedPermulationsData = combinePermData(firstPermulationsData, secondPermulationsData, enrich = enrichValue)

fcTime = Sys.time()
fcCombine = fcTime - spTime
message("Initial permulation combination time: ", fcCombine, attr(fcCombine, "units"))

rm(firstPermulationsData)
rm(secondPermulationsData)


#Do all subsequent combinations
if((startValue+2) < (startValue+permulationNumberValue-1)){
for(i in (startValue+2):(startValue+permulationNumberValue-1)){
  message(i)
  iteratingPermulationsFilename = paste(basePermulationsFilename, i, ".rds", sep="")
  
  if(file.exists(iteratingPermulationsFilename)){
    ipStartTime = Sys.time()
    
    iteratingPermulationsData = readRDS(iteratingPermulationsFilename)
    iteratingPermulationsData = convertLogiToNumericList(iteratingPermulationsData)
    ipLoadTime= Sys.time()
    ipLoad = ipLoadTime - ipStartTime
    message("Iterating permulation load time: ", ipLoad, attr(ipLoad, "units"))
    
    combinedPermulationsData = combinePermData(combinedPermulationsData, iteratingPermulationsData, enrich = enrichValue)
    ipCombineTime = Sys.time()
    ipCombine = ipCombineTime - ipLoadTime
    message("Iterating permulation combination time: ", ipCombine, attr(ipCombine, "units"))
    
    rm(iteratingPermulationsData)
    ipRemoveTime = Sys.time()
    ipRemove = ipRemoveTime - ipCombineTime
    message("Iterating permulation removal time: ", ipRemove, attr(ipRemove, "units"))
    
    message("Added file ", i, " to combination.")
  }else{
    message("Permulation file number ", i, " does not exist. Combining other files.")
  }
}
}

permEndTime = Sys.time()
totalPermTime = permEndTime - timeStart
message(" Total Permulation Combination time: ", totalPermTime, attr(totalPermTime, "units"))

# ---- Save Combined Permulations output as a file ---
if(metacombineValue == F){
  combinedDataFileName = paste(outputFolderName, filePrefix, "Combined", fileTypeString,"PermulationsData", runInstanceValue, ".rds", sep="")
}else{
  combinedDataFileName = paste(outputFolderName, filePrefix, "MetaCombined", fileTypeString, "PermulationsData", runInstanceValue, ".rds", sep="")
}
saveStartTime = Sys.time()
saveRDS(combinedPermulationsData, file = combinedDataFileName)
saveEndTime = Sys.time()
permSavingDuration = saveEndTime - saveStartTime
message("Time to save combine permulations: ", permSavingDuration, attr(permSavingDuration, "units"))



# ---- Calculate the pValues (Deprecated)-----

if(calulateValue){
  message("P Value Calculation Functionality deprecated. Use calculateCombinePermulationsPValues.R")
}
{  #This bracket is to allow for the collapsing of the deprecated function
# -- Get clades correlation --
#This relies on a correlationFile produced by the runPermulations initial scripts. If that file doesn't exist, it skips the step and goes directly to saving the combined permulations.

#cladesCorellationFileName = paste(outputFolderName, filePrefix, "CladesCorrelationFile", sep= "")
#if(file.exists(paste(cladesCorellationFileName, ".rds", sep=""))& calulateValue){
#  cladesCorrelation = readRDS(paste(cladesCorellationFileName, ".rds", sep=""))
#
#  # -- run pValue calculation --
#  permulationPValues = permPValCorReport(cladesCorrelation, combinedPermulationsData)
#  pvalTime = Sys.time()
#  pvalCalcDuration = pvalTime - saveEndTime
#  message("Time to calculate pValues: ", pvalCalcDuration, attr(pvalCalcDuration, "units"))
#  
#  
#  #save the permulations p values
#  permulationPValueFileName = paste(outputFolderName, filePrefix, "CombinedPermulationsPValue", runInstanceValue, ".rds", sep= "")
#  saveRDS(permulationPValues, file = permulationPValueFileName)
#  pValSaveTime = Sys.time()
#  pvalSaveDuration = pValSaveTime - pvalTime
#  message("Time to save pValues: ", pvalSaveDuration, attr(pvalSaveDuration, "units"))
#  
#}else{
#  message("Clades Correlation file does not exist, p-values not calculated. (sisterListGeneration.R)")
#}
}


finalEndTime = Sys.time()
grandTotalTime = finalEndTime - timeStart
message(" Total time: ", grandTotalTime, attr(grandTotalTime, "units"))










# ---- Debug Code ----
#testCombinedPermsData = combinePermData(firstPermulationsData, iteratingPermulationsData, enrich = enrichValue)
#testCombinedPermsDataTwo = combinePermData(testCombinedPermsData, secondPermulationsData, enrich = enrichValue)
#i=3

# --- Permulations cleaning code; for use elsewhere ----

#filtering version
#cleanedfirstPermulationsData = firstPermulationsData
#cleanedfirstPermulationsData$corP = Filter(function(x)!all(is.na(x)), firstPermulationsData$corP)
#cleanedfirstPermulationsData$corRho = Filter(function(x)!all(is.na(x)), firstPermulationsData$corRho)
#cleanedfirstPermulationsData$corStat = Filter(function(x)!all(is.na(x)), firstPermulationsData$corStat)


#converting version
#convertLogiToNumericDataframe = function(inputDataframe){
#  collumnsLogical = sapply(inputDataframe, is.logical)
#  inputDataframe[,collumnsLogical] = lapply(inputDataframe[,collumnsLogical], as.numeric)
#  return(inputDataframe)
#}

#convertLogiToNumericList = function(inputList){
#  framesInList = sapply(inputList, is.data.frame)
#  inputList[framesInList] = lapply(inputList[framesInList], convertLogiToNumericDataframe)
#  return(inputList)
#}
