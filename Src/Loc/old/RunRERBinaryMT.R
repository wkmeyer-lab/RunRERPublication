#Library setup 
.libPaths("/share/ceph/wym219group/shared/libraries/R4") #add path to custom libraries to searched locations
library(RERconverge) #load RERconverge package
library(RERconverge)
library("tools")
source("Src/Reu/cmdArgImport.R")

# ---- USAGE,README ----
# Command line arguments: 
#All strings must be in quotes. 
#If an argument contains a '(' it is evaluated as code.
# 'm=mainTreeFilename.txt or .rds'
# 'p=phenotypeTreeFilename.txt or .rds'
# 'r="filePrefix"'
# 'f=speciesFilterText'
#  v = <T or F>                         This value forces the regeneration of output files which would otherwise not be (RERFile.rds, PathsFile.rds). 


#Takes an inputed main trees multiphylo, binary phenotype tree, file prefix, and optional species filter list
#Outputs a path, RER residuals, and correlation files. 

#Test args: 
args = c('m=data/RemadeTreesAllZoonomiaSpecies.rds', 'r=demoInsectivory')
args = c('m=Data/FishTree.rds', 'r=carnvHerbsOldRebuild', 'v=F')
args = c('m=Data/RemadeTreesAllZoonomiaSpecies.rds', 'r=Domestication', 'v=T')

# ---- Default values if no arguments

#default maintree and phylo location:
mainTreesLocation = "/share/ceph/wym219group/shared/projects/MammalDiet/Zoonomia/RemadeTreesAllZoonomiaSpecies.rds"
binaryPhenotypeTreeLocation = ""

#prefix for output files. Typically the phenotype of interest. 
filePrefix = "test" 

#Put a filter for species here. Default is NULL, which means useSpecies is not applied, and all species are used.
speciesFilter = NULL

# --- Import prefix --- 

args = commandArgs(trailingOnly = TRUE)

#File Prefix
if(!is.na(cmdArgImport('r'))){
  filePrefix = cmdArgImport('r')
}else{
  stop("THIS IS AN ISSUE MESSAGE; SPECIFY FILE PREFIX")
}


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

# ----- load force update argument
forceUpdate = FALSE

#Import if update being forced with argument 
if(!is.na(cmdArgImport('v'))){
  forceUpdate = cmdArgImport('v')
  forceUpdate = as.logical(forceUpdate)
}else{
  paste("Force update not specified, not forcing update")
}

# ---- Command args import ----
{
paste(args)
message(args)

#write.csv(args, file = "Output/args.csv")

#Main Tree Location
mTreesCommandline = grep("^m=",args, value = TRUE) #get a string based on the identifier

#Process the input
if(!is.na(cmdArgImport('m'))){
  mainTreesLocation = cmdArgImport('m')
}else{
  paste("No maintrees arg, using default")                          #Report using default
  message("No maintrees arg, using default")
}



#phenotype tree location
binaryPhenotypeTreeFilename = paste(outputFolderName, filePrefix, "BinaryForegroundTree.rds", sep="") #Make the name of the location a pre-made phenotype tree would have to test for it

if(!is.na(cmdArgImport('p'))){
  binaryPhenotypeTreeLocation = cmdArgImport('p')
}else if(file.exists(paste(binaryPhenotypeTreeFilename))){                  #See if a pre-made binaryTree for this prefix exists 
  binaryPhenotypeTreeLocation = binaryPhenotypeTreeFilename                      #if so, use it 
  paste("Pre-made Phenotype tree found, using pre-made tree.")
  
}else{
  #paste("THIS IS AN ISSUE MESSAGE; SPECIFY PHENOTYPE TREE")
  stop("THIS IS AN ISSUE MESSAGE; SPECIFY PHENOTYPE TREE")
}

#speciesFilter
speciesFilterFileName = paste(outputFolderName, filePrefix, "SpeciesFilter.rds",sep="") #Make the name of the location a pre-made filter would have to test for it

if(!is.na(cmdArgImport('f'))){
  speciesFilter = cmdArgImport('f')
}else if (file.exists(paste(speciesFilterFileName))){                  #See if a pre-made filter for this prefix exists 
  speciesFilter = readRDS(speciesFilterFileName)                       #if so, use it 
  paste("Pre-made filter found, using pre-made filter.")
}else{                                                    
  paste("No speciesFilter arg, using NULL")                           #if not, use no filter
}
}

# ---- MAIN ----

#Read in main tree file 
if(file_ext(mainTreesLocation) == "rds"){
  mainTrees = readRDS(mainTreesLocation)
}else{
mainTrees = readTrees(mainTreesLocation) 
}
#Read in the phenotype tree -- branch length is phenotype value, binary 
if(file_ext(binaryPhenotypeTreeLocation) == "rds"){
  binaryPhenotypeTree = readRDS(binaryPhenotypeTreeLocation)
}else{
  binaryPhenotypeTree = readTrees(binaryPhenotypeTreeLocation) 
}






# ---- RERs ----

RERFileName = paste(outputFolderName, filePrefix, "RERFile.rds", sep= "")

if(!file.exists(paste(RERFileName)) | forceUpdate){
  RERObject = getAllResiduals(mainTrees, useSpecies = speciesFilter, plot = F)
  saveRDS(RERObject, file = RERFileName)
}else{
  RERObject = readRDS(RERFileName)
}

# ---- PATHS ----

pathsFileName = paste(outputFolderName, filePrefix, "PathsFile.rds", sep= "")

if(!file.exists(paste(pathsFileName)) | forceUpdate){
  pathsObject = tree2Paths(binaryPhenotypeTree, mainTrees, binarize=T, useSpecies = speciesFilter)
  saveRDS(pathsObject, file = pathsFileName)
}else{
  pathsObject = readRDS(pathsFileName)
}


# ---- Correlate ----

outputFileName = paste(outputFolderName, filePrefix, "CorrelationFile", sep= "")

correl = correlateWithBinaryPhenotype(RERObject, pathsObject, min.sp =35)

#head(correl[order(correl$p.adj),])

write.csv(correl, file= paste(outputFileName, ".csv", sep=""), row.names = T, quote = F)
saveRDS(correl, paste(outputFileName, ".rds", sep=""))
   
paste("Script successful.")
message("Script successful.")
# ----- Arguments integration testing -----




