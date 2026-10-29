







#  SCRIPT DOES NOT WORK!!!!
# The categorical drop tip is unable to handle sufficiently large clades of a single dropped phenotype













clusterRun = F
clusterRun = T
if(clusterRun){.libPaths("/share/ceph/wym219group/shared/libraries/R4")} #add path to custom libraries to searched locations
library(RERconverge)
library(tools)
source("Src/Reu/cmdArgImport.R")



# r = filePrefix                                         This is a prefix used to organize and separate files by analysis run. Always required. 
# v = <T or F>                                           This prefix is used to force the regeneration of the script's output, even if the files already exist. Not required, not always used.
# d = CategoriesToDrop                                   This sets the location of the maintrees file
# m = mainTreesLoaction
args = c('r=ComplexDietCentralAnalysisSimplifyStrictPred2', 'v=F', 'd=zMixedPredator', 'm=data/zoonomiaAllMammalsTrees.rds')

# --- Standard start-up code ---
if(clusterRun){args = commandArgs(trailingOnly = TRUE)}
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



# --- Load arguments ----
mainTreesLocation = "/share/ceph/wym219group/shared/projects/MammalDiet/Zoonomia/RemadeTreesAllZoonomiaSpecies.rds"


#Read in the categories to drop 
if(!is.na(cmdArgImport('d'))){
  droppedCats = cmdArgImport('d')
}else{
  stop("THIS IS AN ISSUE MESSAGE; SPECIFY DROPPED PHENOTYPES")
}
#load Maintrees 
#MainTrees Location
if(!is.na(cmdArgImport('m'))){
  mainTreesLocation = cmdArgImport('m')
}else{
  message("No maintrees arg, using default")
}

if(file_ext(mainTreesLocation) == "rds"){
  if(!exists("mainTrees")){mainTrees = readRDS(mainTreesLocation)}
}else{
  if(!exists("mainTrees")){mainTrees = readTrees(mainTreesLocation)} 
}


# --- Main code ---
producePrunedFiles = function(droppedCategories = droppedCats, plot=F){
  source("Src/Reu/CategoricalDropTip.R")
  #Read in the main versions
  categoricalTreeFilename = paste(outputFolderName, filePrefix, "CategoricalTree.rds", sep="") #make a filename based on the prefix
  phenotypeTree = readRDS(categoricalTreeFilename)
  
  #categoricalCommonTreeFilename = paste(outputFolderName, filePrefix, "CategoricalCommonTree.rds", sep="") #make a filename based on the prefix
  #scientificCategoricalTreeFilename = paste(outputFolderName, filePrefix, "CategoricalScientificTree.rds", sep="") #make a filename based on the prefix
  #commonTree = readRDS(categoricalCommonTreeFilename) #these are canceled to avoid needing to pass in the arugments for the name convert function
  #scientificTree = readRDS(scientificCategoricalTreeFilename)
  
  
  phenotypeVectorFilename = paste(outputFolderName, filePrefix, "CategoricalPhenotypeVector.rds",sep="") #make a filename based on the prefix
  phenotypeVector = readRDS(phenotypeVectorFilename)                       #save the phenotype vector
  
  speciesFilterFilename = paste(outputFolderName, filePrefix, "SpeciesFilter.rds",sep="") #set a filename for the species filter based on the prefix 
  speciesFilter = readRDS(speciesFilterFilename)
  
  
  speciesToDrop = names(phenotypeVector)[which(phenotypeVector %in% droppedCategories)]
  
  #Drop category from tree and vectors 
  prunedCategoricalTree = categoricalDropTip(phenotypeTree, speciesToDrop)
  
  #prunedCommonTree = categoricalDropTip(phenotypeTree, speciesToDrop)
  #prunedScientificTree = categoricalDropTip(phenotypeTree, speciesToDrop)
  
  
  prunedPhenotypeVector = phenotypeVector[-which(names(phenotypeVector)%in% speciesToDrop)]
  prunedSpeciesFilter = speciesFilter[-which(speciesFilter %in% speciesToDrop)]
  
  
  
  # ---- Save the files ----- 
  
  categoricalPrunedTreeFilename = paste(outputFolderName, filePrefix, "PrunedCategoricalTree.rds", sep="") #make a filename based on the prefix
  #categoricalPrunedCommonTreeFilename = paste(outputFolderName, filePrefix, "PrunedCategoricalCommonTree.rds", sep="") #make a filename based on the prefix
  #scientificPrunedCategoricalTreeFilename = paste(outputFolderName, filePrefix, "PrunedCategoricalScientificTree.rds", sep="") #make a filename based on the prefix
  
  saveRDS(prunedCategoricalTree, categoricalPrunedTreeFilename)
  #saveRDS(prunedCommonTree, categoricalPrunedCommonTreeFilename)
  #saveRDS(prunedScientificTree, scientificPrunedCategoricalTreeFilename)
  
  phenotypeVectorPrunedFilename = paste(outputFolderName, filePrefix, "PrunedCategoricalPhenotypeVector.rds",sep="") #make a filename based on the prefix
  speciesFilterPrunedFilename = paste(outputFolderName, filePrefix, "PrunedSpeciesFilter.rds",sep="") #set a filename for the species filter based on the prefix 
  
  saveRDS(prunedPhenotypeVector, phenotypeVectorPrunedFilename)
  saveRDS(prunedSpeciesFilter, speciesFilterPrunedFilename)
  
  if(plot){
    source("Src/Reu/ZoonomTreeNameToCommon.R")
    #commonMainTrees = mainTrees
    #commonMainTrees$masterTree = ZoonomTreeNameToCommon(commonMainTrees$masterTree, manualAnnotLocation = spreadSheetLocation, tipCol = nameColumn)
    palette(c( "darkgreen", "darkblue","black", "red", "yellow"))
    
    treeImageFilename = paste(outputFolderName, filePrefix, "CategoryPrunedCategoricalTree.pdf", sep="") #make a filename based on the prefix
    
    categoryLabels = unique(prunedPhenotypeVector)[order(unique(prunedPhenotypeVector))]
    
    pdf(treeImageFilename, height = length(prunedPhenotypeVector)/18, width = 10)                     #make a pdf to store the plot, sized based on tree size
    
    
    #plotTreeCategorical(prunedCommonTree, categoryLabels, master = commonMainTrees$masterTree)
    plotTreeCategorical(prunedCategoricalTree, categoryLabels, master = mainTrees$masterTree)
    dev.off()  
    
  }
}

producePrunedFiles(droppedCats, plot = T)

prunedCategoricalTree$ed

commonMainTrees = mainTrees
commonMainTrees$masterTree = ZoonomTreeNameToCommon(commonMainTrees$masterTree, manualAnnotLocation = spreadSheetLocation, tipCol = nameColumn)
