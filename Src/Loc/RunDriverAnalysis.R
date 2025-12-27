clusterRun = F
clusterRun = T
if(clusterRun){.libPaths("/share/ceph/wym219group/shared/libraries/R4")} #add path to custom libraries to searched locations
library(RERconverge)
source("Src/Reu/cmdArgImport.R")
source("Src/Reu/AssessRERDriver.R")
source("Src/Reu/AssessGoCategoryDriver.R")

# - Usage:
# This script runs the driver analysis code on a pairwise comparision of a categorical phenotype, once the appropriate binary analyses have been compelted. 
# This script assumes that the binary analysis are named with the format "CategoricalBinary[category]Tree" as their file prefix. 
# 
# Arguments: 
# r = filePrefix                                                               This is a prefix used to organize and separate files by analysis run. Always required. 
# c = c('Category1', 'Category2')     This determines which two categories are in the pairwise comparision, must be entered in alphabetical order 
# g = c("geneset", "optionaladditionalgeneset") OR NULL                     This determines which Go Category set has directionality assessed. only one can be run at a time. 


args = c("r=CategoricalInsVertivoreTree", "c=c('Herbivore', 'Vertivore')", "g=KeggReactome")
args = c("r=CategoricalInsVertivoreTree", "c=c('Herbivore', 'Insectivore')", "g=KeggReactome")
args = c("r=CategoricalInsVertivoreTree", "c=c('Carnivore', 'Herbivore')", "g=KeggReactome")



# -- Standard Startup code -- 
if(clusterRun)args = commandArgs(trailingOnly = TRUE)
{  # Bracket used for collapsing purposes
  #File Prefix
  if(!is.na(cmdArgImport('r'))){
    filePrefix = cmdArgImport('r')
  }else{
    stop("THIS IS AN ISSUE MESSAGE; SPECIFY FILE PREFIX")
  }
  
  #  Output Directory 
  if(!dir.exists("Output")){                                      #Make output directory if it does not exist
    dir.create("Output")
  }
  outputFolderNameNoSlash = paste("Output/",filePrefix, sep = "") #Set the prefix sub directory
  if(!dir.exists(outputFolderNameNoSlash)){                       #create that directory if it does not exist
    dir.create(outputFolderNameNoSlash)
  }
  outputFolderName = paste("Output/",filePrefix,"/", sep = "")
  
  #  Force update argument
  forceUpdate = FALSE
  if(!is.na(cmdArgImport('v'))){                                 #Import if update being forced with argument 
    forceUpdate = cmdArgImport('v')
    forceUpdate = as.logical(forceUpdate)
  }else{
    message("Force update not specified, not forcing update")
  }
}




# -- Argument imports -- 
categories = NULL
geneSets = "KeggReactome"

#Categories
if(length(cmdArgImport('c'))==2){
  categories = cmdArgImport('c')
}else{
  stop("Provide exactly two categories; this script must be run seperately for each pairwise comparision.")
}

#genesets 
if(!all(is.na(cmdArgImport('g')))){
  geneSets = cmdArgImport('g')
}else{
  message("No geneSet provided, uses default")
}


# ---- Generate the function inputs --- 

pairwiseFolderName = paste(categories, collapse = '-')
binaryNameOne = paste0("CategoricalBinary", categories[1], "Tree")
binaryNameTwo = paste0("CategoricalBinary", categories[2], "Tree")




driverTable = AssessRERDriver(filePrefix, pairwiseFolderName, binaryNameOne, categories[1], binaryNameTwo, categories[2])

driverTableFilename = paste0(outputFolderName, pairwiseFolderName, "/", filePrefix, pairwiseFolderName, "DriverTable")
write.csv(driverTable, paste0(driverTableFilename, ".csv"))
saveRDS(driverTable, paste0(driverTableFilename, ".rds"))



for(i in 1:length(geneSets)){
  GoDriver= AssessGoCategoryDriver(filePrefix, pairwiseFolderName, geneSets[i], 0.1, 1, F)
  GoDriverTableFilename = paste0(outputFolderName, pairwiseFolderName, "/", filePrefix, pairwiseFolderName, "GoDriverTable-", geneSets[i])
  
  
  write.csv(GoDriver, paste0(GoDriverTableFilename, ".csv"))
  saveRDS(GoDriver, paste0(GoDriverTableFilename, ".rds"))
  
}


