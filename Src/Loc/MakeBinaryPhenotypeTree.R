# -- Libraries 
.libPaths("/share/ceph/wym219group/shared/libraries/R4") #DO NOT RUN LOCALLY. add Cluster path to custom libraries to searched locations

library(RERconverge)
library(tools)
source("Src/Reu/cmdArgImport.R")

# -- Usage:
# This Script creates a binary tree from a csv file containing (at least) the following columns: 
# FaName: The name of the species as it appears in the Mutliphylo 
# Common.Name.Or.Group: The common name of the species
# Species.Name: The scientific name of the species
# [phenotypeCollumn]: A column (name set by argument) which contains the pehnotype data 
# (optional)[ScreenColumn]: A column which can be used to select a subset of the data. If a screen column is specified by argument, only species with a "1" in that column will be used to make the tree. 
# If this file is not named "manualAnnotationsSheet.csv", some functions may break. 

# -- Command arguments list
# r = filePrefix                           This is a prefix used to organize and separate files by analysis run. Always required. 
# v = <T or F>                             This prefix is used to force the regeneration of the script's output, even if the files already exist. Not required, not always used.
# m = mainTreeFilename.txt or rds          This is the location of the main trees file 
# a = annotationFileLocation.csv           This is the location of the annotation csv. If not the default "manualAnnotationsSheet.csv", some functions in other scripts may break. 
# p = phenotypeColumn                      This is the name of the column contianing the phenotype data 
# t = <"uni" OR "bi">                      This set unidirectional or bidirectional transition 
# c = <"ancestral" OR "all" OR "terminal"> This sets the clade type to be used in tree creation
# w = <T or F>                             This sets if the tree should be weighted 
# s = "screenColumn"                       This sets the column with the screening data
# n = "nameColumn"                         This sets the column with the tip names as they appear in the maintrees file. 
# z = <minimum branch length>              This sets the minimum branch length for terminal branches in the master tree. Branches shorter than this will be removed. 
# x = "pruningPrefrenceColumn"             This sets a column, where if the value is 1, the tip will be preferentially kept. If the value is TRUE, the tip will never be pruned.
# y = "c('unprunedtip1', 'unprunedtip2')"  This allows you to add a list of specific tips to not be dropped during pruning. Must use the tip name, not common name. 


#----------------
args = c('m=data/RemadeTreesAllZoonomiaSpecies.rds', "r=EcholocationUpdate2", "t=bi", "p=Echolocation", "c=all", "v=T", "s=Laurasiatheria")
args = c("m=data/FirstExpressionTrees.rds", "r=LiverExpression", "p=Carnivory", "t=bi", "c=all", "w=F", "v=T", "a=Data/ExpressionAnnots.csv")
args = c("m=data/RemadeTreesAllZoonomiaSpecies.rds", "r=CVHApplesToApples", "p=Carnivory", "t=bi", "c=all", "w=F", "v=T", "a=Data/ExpressionAnnots.csv")

args = c("m=Data/UNICORNsDemo.txt", "r=ActueLoafs", "a=Results/ToothData.csv", "p=FCT_AL")


args = c("m=data/RemadeTreesAllZoonomiaSpecies.rds", "r=CvHNew", "p=CarnFish_Herbs", "t=bi", "c=all", "w=F", "v=T", "s=Laurasiatheria", "a=Data/MergedData.csv", "n=Zoonomia")

args = c("m=data/newHillerMainTrees.rds", "r=CvHNewHiller", "p=CarnFish_Herbs", "t=bi", "c=all", "w=F", "v=T", "a=Data/MergedData.csv", "n=HillerName")
args = c("m=data/newHillerMainTrees.rds", "r=CvHNewHiller", "p=CarnFish_Herbs", "t=bi", "c=all", "w=F", "v=T", "s=Laurasiatheria", "a=Data/MergedData.csv", "n=HillerName")

args = c("m=data/newHillerMainTrees.rds", "r=CvHNewHiller", "p=CarnFish_Herbs", "t=bi", "c=all", "w=F", "v=T", "a=Data/MergedData.csv", "n=HillerName")
args = c("m=data/newHillerMainTrees.rds", "r=CvHNewHiller", "p=CarnFish_Herbs", "t=bi", "c=all", "w=F", "v=T", "s=Laurasiatheria", "a=Data/MergedData.csv", "n=HillerName")

args = c("m=data/RemadeTreesAllZoonomiaSpecies.rds", "r=CvAllZoonomia", "p=CarnBinary", "t=bi", "c=all", "w=F", "v=T", "s=LaurasiatheriaFilter", "a=Data/MergedData.csv", "n=Zoonomia")
args = c("m=data/RemadeTreesAllZoonomiaSpecies.rds", "r=HvAllZoonomia", "p=HerbBinary", "t=bi", "c=all", "w=F", "v=T", "s=LaurasiatheriaFilter", "a=Data/MergedData.csv", "n=Zoonomia")
args = c("m=data/RemadeTreesAllZoonomiaSpecies.rds", "r=CIvAllZoonomia", "p=CarnInsectBinary", "t=bi", "c=all", "w=F", "v=T", "s=LaurasiatheriaFilter", "a=Data/MergedData.csv", "n=Zoonomia")

args = c("m=data/RemadeTreesAllZoonomiaSpecies.rds", "r=CvHNew", "p=CarnFish_Herbs", "t=bi", "c=all", "w=F", "v=T", "s=Laurasiatheria", "a=Data/MergedData.csv", "n=Zoonomia", "z=0.01", "x=HillerZoonomiaOverlap")

args = c("r=CVHDemo", "m=data/RemadeTreesAllZoonomiaSpecies.rds", "a=Data/MergedData.csv", "p=CarnFish_Herbs", "t=bi", "n=Zoonomia")

args = c("r=harshalBats", "m=Data/batDemoMaintrees.rds", "a=Data/harshalpheno.csv", "p=Phenotype", "t=bi", "n=Species", "v=T")


#pruned reruns
args = c("m=data/RemadeTreesAllZoonomiaSpecies.rds", "r=CvHZoonomia", "p=CarnFish_Herbs", "t=bi", "c=all", "w=F", "v=T", "s=Laurasiatheria", "a=Data/MergedData.csv", "n=Zoonomia", "z=0.01", "x=HillerZoonomiaOverlap")
args = c("m=data/RemadeTreesAllZoonomiaSpecies.rds", "r=CIvAllZoonomia", "p=CarnInsectBinary", "t=bi", "c=all", "w=F", "v=T", "s=LaurasiatheriaFilter", "a=Data/MergedData.csv", "n=Zoonomia", "z=0.01", "x=HillerZoonomiaOverlap")
args = c("m=data/RemadeTreesAllZoonomiaSpecies.rds", "r=CvAllZoonomia", "p=CarnBinary", "t=bi", "c=all", "w=F", "v=T", "s=LaurasiatheriaFilter", "a=Data/MergedData.csv", "n=Zoonomia", "z=0.01", "x=HillerZoonomiaOverlap")
args = c("m=data/RemadeTreesAllZoonomiaSpecies.rds", "r=HvAllZoonomia", "p=HerbBinary", "t=bi", "c=all", "w=F", "v=T", "s=LaurasiatheriaFilter", "a=Data/MergedData.csv", "n=Zoonomia", "z=0.01", "x=HillerZoonomiaOverlap")

args = c("m=data/newHillerMainTrees.rds", "r=CvHHiller", "p=CarnFish_Herbs", "t=bi", "c=all", "w=F", "v=T", "s=Laurasiatheria", "a=Data/MergedData.csv", "n=HillerName", "z=0.01", "x=HillerZoonomiaOverlap")
args = c("m=data/newHillerMainTrees.rds", "r=CIvAllHiller", "p=CarnInsectBinary", "t=bi", "c=all", "w=F", "v=T", "s=LaurasiatheriaFilter", "a=Data/MergedData.csv", "n=HillerName", "z=0.01", "x=HillerZoonomiaOverlap")
args = c("m=data/newHillerMainTrees.rds", "r=CvAllHiller", "p=CarnBinary", "t=bi", "c=all", "w=F", "v=T", "s=LaurasiatheriaFilter", "a=Data/MergedData.csv", "n=HillerName", "z=0.01", "x=HillerZoonomiaOverlap")
args = c("m=data/newHillerMainTrees.rds", "r=HvAllHiller", "p=HerbBinary", "t=bi", "c=all", "w=F", "v=T", "s=LaurasiatheriaFilter", "a=Data/MergedData.csv", "n=HillerName", "z=0.01", "x=HillerZoonomiaOverlap")


# --- Standard start-up code ---
args = commandArgs(trailingOnly = TRUE)
{  # Bracket used for collapsing purposes
  #File Prefix
  if(!is.na(cmdArgImport('r'))){                                                #This cmdArgImport script is a way to import arguments from the command line. 
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

# --- Argument Imports ---
{ # Bracket used for collapsing purposes
  # Defaults
  mainTreesLocation = "/share/ceph/wym219group/shared/projects/MammalDiet/Zoonomia/RemadeTreesAllZoonomiaSpecies.rds"
  annotationsLocation = "Data/manualAnnotationsSheet.csv"
  phenotypeColumn = "ERRORDEFAULT"
  transitionValue = "Default"
  cladeValue = "Default"
  weightValue = FALSE
  useScreen = F
  screenCollumn = NA
  nameCollumn = "tipName"
  usingPruning = F
  pruningCutoff = NA
  pruningPrefrenceColumn = NA
  
  #Main Tree Location
  if(!is.na(cmdArgImport('m'))){
    mainTreesLocation = cmdArgImport('m')
  }else{
    message("No maintrees arg, using default")
  }
  
  #Annotations Location
  if(!is.na(cmdArgImport('a'))){
    annotationsLocation = cmdArgImport('a')
  }else{
    message("No maintrees arg, using default 'manualAnnotationsSheet.csv' ")
  }
  
  #Phenotype Column
  if(!is.na(cmdArgImport('p'))){
    phenotypeColumn = cmdArgImport('p')
  }else{
    stop("THIS IS AN ISSUE MESSAGE; SPECIFY PHENOTYPE COLLUMN")
  }
  
  #Transition value
  if(!is.na(cmdArgImport('t'))){
    transitionValue = cmdArgImport('t')
  }else{
    message("Using default bidirectional transistion")
  }
  
  #Clade value
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
  
  #Screen Column 
  if(!is.na(cmdArgImport('s'))){
    useScreen = T
    screenCollumn = cmdArgImport('s')
  }else{
    message("Screen Column not specified, not using screen column.")
  }

  #Name column 
  if(!is.na(cmdArgImport('n'))){
    nameCollumn = cmdArgImport('n')
  }else{
    message("Name Column not specified, using 'tipName'.")
  }
  
  #Pruning cutoff
  if(!is.na(cmdArgImport('z'))){
    usingPruning = T
    pruningCutoff = cmdArgImport('z')
  }else{
    message("Pruning Cutoff not specified, not pruning tree.")
  }
  
  #PruningPrefrenceColumn 
  if(!is.na(cmdArgImport('x'))){
    pruningPrefrenceColumn = cmdArgImport('x')
  }else{
    if(usingPruning){message("No pruning prefrence column specified")}
  }
  
  #ManualPruningProtections
  if(!is.na(cmdArgImport('y'))){
    manualPruningProtections = cmdArgImport('y')
  }else{
    if(usingPruning){message("No manually protected species specified")}
  }
}

#                   ------- Code Body -------- 

# - Import Files -
if(file_ext(mainTreesLocation) == "rds"){
  mainTrees = readRDS(mainTreesLocation)
}else{
  mainTrees = readTrees(mainTreesLocation) 
}
manualAnnots = read.csv(annotationsLocation)

# - Species Filter - 
speciesFilterFilename = paste(outputFolderName, filePrefix, "SpeciesFilter.rds",sep="")

if(!file.exists(speciesFilterFilename) | forceUpdate){                          #If not using exisitng species filter
  relevantSpecies = manualAnnots[ manualAnnots[[phenotypeColumn]] %in% c(0,1),] #Set the relevant species to species with a 1 or a 0 in the annotation column 
  if(useScreen){                                                                #if using a screening column, 
    relevantSpecies = relevantSpecies[relevantSpecies[screenCollumn] == 1,]     #only include species with a 1 in the screen column
  }                                                                             
  relevantSpeciesNames = relevantSpecies[[nameCollumn]]                                 #use the species names 
  relevantNamesNA = which(is.na(relevantSpeciesNames))
  if(length(relevantNamesNA) > 0){  relevantSpeciesNames = relevantSpeciesNames[-relevantNamesNA]}

  
  if(usingPruning){
    source("Src/Reu/autoPruner.R")
    pruningProtectionSpecies = NA
    if(!is.na(pruningPrefrenceColumn)){
      pruningProtectionRows = manualAnnots[which(as.logical(manualAnnots[[pruningPrefrenceColumn]])),]
      pruningProtectionSpecies = pruningProtectionRows[[nameCollumn]]
      if(all(is.logical(manualAnnots[[pruningPrefrenceColumn]]))){
        pruningProtection = T
      }else{ 
        pruningProtection = F
      }
    }
    allProtectedSpecies = append(pruningProtectionSpecies, manualPruningProtections)
    
    
    workingTree = mainTrees$masterTree
    workingTree = drop.tip(workingTree, which(!workingTree$tip.label %in% relevantSpeciesNames))
    
    fewGeneSpecies = dropFewGeneSpecies(mainTrees, workingTree, nameConversionColumn = nameCollumn, nameConversionData = annotationsLocation)
    workingTree = drop.tip(workingTree, fewGeneSpecies)
    
    pruningFilename = paste(outputFolderName, filePrefix, "PruningTree.pdf", sep="")
    pdf(pruningFilename, width=16, height = 14)
    prunedTree = autopruner(workingTree, dropValue = pruningCutoff, tipsToKeep = pruningProtectionSpecies, nameConversionColumn = nameCollumn, nameConversionData = annotationsLocation, preDroppedTips = fewGeneSpecies)
    if(!pruningProtection){
      prunedTree = autopruner(prunedTree, dropValue = pruningCutoff, tipsToKeep = manualPruningProtections, nameConversionColumn = nameCollumn, nameConversionData = annotationsLocation, preDroppedTips = droppedTips, originalTree = workingTree)
    }
    dev.off()
    
    prunedSpecies = relevantSpeciesNames[!relevantSpeciesNames %in% prunedTree$tip.label]
    relevantSpeciesNames = relevantSpeciesNames[-which(relevantSpeciesNames %in% prunedSpecies)]
  }
  
  saveRDS(relevantSpeciesNames, file = speciesFilterFilename)
  
}else{
  relevantSpecieslist = readRDS(speciesFilterFilename)
  relevantSpeciesNames = relevantSpecieslist
  relevantSpecies = manualAnnots[ manualAnnots[[nameCollumn]] %in% relevantSpecieslist,]
}

# - Setup foreground Species --

foregroundSpeciesAnnot = relevantSpecies[ relevantSpecies[[phenotypeColumn]] %in% 1,] #set the foreground species to spcies with a 1 in the annotation column

foregroundNames = foregroundSpeciesAnnot[[nameCollumn]]
foregroundNamesNA = which(is.na(foregroundNames))
if(length(foregroundNamesNA) > 0){  foregroundNames = foregroundNames[-foregroundNamesNA]}


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

# - Make the tree - 
binaryForegroundTreeOutput = do.call(foreground2Tree, f2tInputList)


# ---- Saving and visualization ----

# - Save the Tree - 
binaryTreeFilename = paste(outputFolderName, filePrefix, "BinaryForegroundTree.rds", sep="")
saveRDS(binaryForegroundTreeOutput, file = binaryTreeFilename)
binaryTreeFilename = paste(outputFolderName, filePrefix, "BinaryTree.rds", sep="")
saveRDS(binaryForegroundTreeOutput, file = binaryTreeFilename)
binaryNewickFilename = paste(outputFolderName, filePrefix, "NewickTree.nwk", sep="")
write.tree(binaryForegroundTreeOutput, file = binaryNewickFilename)

# - Read back in tree and print to pdf - 
readTest = readRDS(binaryTreeFilename)
testTreeDisplayable = readTest
testTreeDisplayable$edge.length = replace(testTreeDisplayable$edge.length, testTreeDisplayable$edge.length==0, 0.5)
testTreeDisplayable$edge.length = replace(testTreeDisplayable$edge.length, testTreeDisplayable$edge.length==1, 4)

source("Src/Reu/plotBinaryTree.R")

binaryTreePdfname = paste(outputFolderName, filePrefix, "BinaryForegroundTree.pdf", sep="")
pdf(binaryTreePdfname, width=8, height = 14)
plotBinaryTree(mainTrees, readTest, foregroundNames, mainTitle = paste(filePrefix, "Binary", "Foreground", "Tree"), tipColumn = nameCollumn)
plotBinaryTree(mainTrees, readTest, foregroundNames, convertNames = F, mainTitle = paste(filePrefix, "Binary", "Foreground", "Tree"), tipColumn = nameCollumn)
plotTree(testTreeDisplayable)
dev.off()


