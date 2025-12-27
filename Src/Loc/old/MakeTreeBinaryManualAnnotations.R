.libPaths("/share/ceph/wym219group/shared/libraries/R4") #add path to custom libraries to searched locations
library(RERconverge) #load RERconverge package
library(RERconverge)
library(tools)
library(devtools)
source(file = "Src/Reu/cmdArgImport.R")

# ---- USAGE,README ----
# Command line arguments: 
#All strings must be in quotes. 
#If an argument contains a '(' it is evaluated as code.
# 'm=mainTreeFilename.txt or .rds'       This is the main tree's file location
# 'a="annotCollumn"'                     This is the column in the manual annotations spreadsheet to use
# 'r="filePrefix"'                       This is the file prefix for this execution
# 't=<"uni"/"bi>"                        Sets if transitions are unidirectional or bidirectional         
# 'c=<"ancestral"/"all"/"terminal">'     Sets the clade type 
# 'w=<T/F>"                              Sets if the foreground2tree is weighted or not 
#  v = <T or F>                         #This value forces the regeneration of output files Which otherwise would not be (SpeciesFilter.rds). 
#  s = "screenCollumn"                   This is a collumn which must have a value of 1 for the species to be included. 

# ---- Default Arguments ----
#Default main tree file location

mainTreesLocation = "/share/ceph/wym219group/shared/projects/MammalDiet/Zoonomia/RemadeTreesAllZoonomiaSpecies.rds"

#local computer debug version:
# mainTreesLocation = "data/RemadeTreesAllZoonomiaSpecies.rds"
args = c("m=data/RemadeTreesAllZoonomiaSpecies.rds", "r=Domestication", "a=DomesticationBinary", "t=bi", "c=terminal", "w=F", "v=T")
args = c("m=data/RemadeTreesAllZoonomiaSpecies.rds", "r=CVHRemake", "a=CarnFish_Herbs", "s=Laurasiatheria", "t=bi", "c=terminal", "w=F", "v=T")

#Other defaults if not specified
annotCollumn = "ERRORDEFAULT"
filePrefix = "ERRORDEFAULT"
transisitionValue = "Default"
cladeValue = "Default"
weightValue = FALSE
useScreen = F
screenCollumn = NA




# ---- Command args import ----

args = commandArgs(trailingOnly = TRUE)
paste(args)
message(args)

#write.csv(args, file = "Output/MakePathsArgs.csv")

#Main Tree Location
if(!is.na(cmdArgImport('m'))){
  mainTreesLocation = cmdArgImport('m')
}else{
  paste("No maintrees arg, using default")                          #Report using default
  message("No maintrees arg, using default")
}

#File Prefix
if(!is.na(cmdArgImport('r'))){
  filePrefix = cmdArgImport('r')
}else{
  stop("THIS IS AN ISSUE MESSAGE; SPECIFY FILE PREFIX")
}

#Annots Collumn
if(!is.na(cmdArgImport('a'))){
  annotCollumn = cmdArgImport('a')
}else{
  stop("THIS IS AN ISSUE MESSAGE; SPECIFY ANNOTATION COLLUMN")
}

#Transition value
if(!is.na(cmdArgImport('t'))){
  transitionValue = cmdArgImport('t')
}else{
  message("Using default bidirectional transistion")
}

#clade value
if(!is.na(cmdArgImport('c'))){
  cladeValue = cmdArgImport('c')
}else{
  message("Using default all clade")
}

#Weight Value
if(!is.na(cmdArgImport('w'))){
  weightValue = cmdArgImport('w')
}else{
  message("Weight = false")
}

#Screen Collumn
if(!is.na(cmdArgImport('s'))){
  useScreen = T
  screenCollumn = cmdArgImport('s')
}else{
  message("No screen collumn used.")
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


# -------- Make Paths Main Code ---------------


# -- Import a tree --
if(file_ext(mainTreesLocation) == "rds"){
  mainTrees = readRDS(mainTreesLocation)
}else{
  mainTrees = readTrees(mainTreesLocation) 
}

#make a Results directory if one doesn't exist
if(!dir.exists("Results")){
  dir.create("Results")
}

# -- Import manual annotations CSv -- 

manualAnnots = read.csv("Data/manualAnnotationsSheet.csv")

# ---- if a pre-defined species filter exists, use that. If not, create one. 
speciesFilterFilename = paste(outputFolderName, filePrefix, "SpeciesFilter.rds",sep="")

if(!file.exists(speciesFilterFilename) | forceUpdate){
  # --- subset the manual annots to only those with data in the collumn
  relevantSpecies = manualAnnots[ manualAnnots[[annotCollumn]] %in% c(0,1),]
  if(useScreen){
    relevantSpecies = relevantSpecies[ relevantSpecies[screenCollumn] %in% 1, ]
  }
  relevantSpeciesNames = relevantSpecies$FaName
  
  # -- output a species filter for this fileprefix 
  
  saveRDS(relevantSpeciesNames, file = speciesFilterFilename)
}else{
  relevantSpecieslist = readRDS(speciesFilterFilename)
  relevantSpeciesNames = relevantSpecieslist
  relevantSpecies = manualAnnots[ manualAnnots[["FaName"]] %in% relevantSpecieslist,]
}
# -- Setup foreground Species -- 

foregroundSpeciesAnnot = relevantSpecies[ relevantSpecies[[annotCollumn]] %in% 1,]

foregroundNames = foregroundSpeciesAnnot$FaName

foregroundFilename = paste(outputFolderName, filePrefix, "BinaryTreeForegroundSpecies.rds", sep="")
saveRDS(foregroundNames, file = foregroundFilename)


# -- set arguments for foreground2Trees --
f2tInputList = list(foregroundNames, mainTrees, useSpecies = relevantSpeciesNames)
#Transition
if(transitionValue == "uni"){
  f2tInputList[["transition"]]= "unidirectional"
  message("Unidirectional transition")
}else{
  f2tInputList[["transition"]]= "bidirectional"
  message("Bidirectional transition")
}

#clade
if(cladeValue == "ancestral"){
  f2tInputList[["clade"]] = "ancestral"
  message("ancestral clade")
}else if(cladeValue == "terminal"){
  f2tInputList[["clade"]] = "terminal"
  message("terminal clade")
}else{
  f2tInputList[["clade"]] = "all"
  message("all clade")
}

#Weight
if(weightValue ==TRUE || weightValue == 't'){
  f2tInputList[["weighted"]] = TRUE
  message("Weighted = TRUE")
}else{
  f2tInputList[["weighted"]] = FALSE
  message("Weighted = False")
}


#Make binary tree output 
#paste(f2tInputList)
#f2tInputList
binaryForegroundTreeOutput = do.call(foreground2Tree, f2tInputList)



#Save that output 
binaryTreeFilename = paste(outputFolderName, filePrefix, "BinaryForegroundTree.rds", sep="")
saveRDS(binaryForegroundTreeOutput, file = binaryTreeFilename)
binaryTreeFilename = paste(outputFolderName, filePrefix, "BinaryTree.rds", sep="")
saveRDS(binaryForegroundTreeOutput, file = binaryTreeFilename)

#Read back in the outputted tree as a test 
readTest = readRDS(binaryTreeFilename)
testTreeDisplayable = readTest
testTreeDisplayable$edge.length = replace(testTreeDisplayable$edge.length, testTreeDisplayable$edge.length==0, 0.5)
testTreeDisplayable$edge.length = replace(testTreeDisplayable$edge.length, testTreeDisplayable$edge.length==1, 4)

binaryTreePdfname = paste(outputFolderName, filePrefix, "BinaryForegroundTree.pdf", sep="")
pdf(binaryTreePdfname, width=8, height = 14)
plotTreeHighlightBranches(testTreeDisplayable, hlspecies=which(readTest$edge.length==1), hlcols="blue",)
dev.off()