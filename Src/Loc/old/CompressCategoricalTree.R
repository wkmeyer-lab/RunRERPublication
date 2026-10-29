# -- Libraries 
clusterRun = F
clusterRun = T
if(clusterRun){.libPaths("/share/ceph/wym219group/shared/libraries/R4")} #add path to custom libraries to searched locations
library(RERconverge)
library(tools)
source("Src/Reu/cmdArgImport.R")
source("Src/Reu/ZoonomTreeNameToCommon.R")

# -- Usage:
# This script takes a categorical tree with a large number of categories, and combines them into fewer categories based on used specifications. 
# In theory, this script could be used on any spreadsheet, so long as the column containing the tip.labels is specified using the n argument, and the column with common names is named "CommonName". 

# -- Command arguments list
# r = filePrefix                                         This is a prefix used to organize and separate files by analysis run. Always required. 
# v = <T or F>                                           This prefix is used to force the regeneration of the script's output, even if the files already exist. Not required, not always used.
# m = mainTreeFilename.txt or .rds                       This sets the location of the maintrees file
# c = filePrefixOfRefrenceTree                           This sets the location of the Refrence Tree to be compressed 
# u = list(c("replace1", "with1"),c("replace2, with2"))  This sets the phenotypes to be merge together 


#---------------------------------
args = c('r=TestCompress', 'm=data/zoonomiaAllMammalsTrees.rds', 'u=list( c("Frugivore", "Glucivore"), c("Nectarivore", "Glucivore"))')
args = c('r=TestCompress', 'm=data/zoonomiaAllMammalsTrees.rds', 'u=list( c("Frugivore", "Glucivore"), c("Nectarivore", "Glucivore"))')




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

# --- Argument Imports ---
{ # Bracket used for collapsing purposes
  # Defaults
  mainTreesLocation = "/share/ceph/wym219group/shared/projects/MammalDiet/Zoonomia/RemadeTreesAllZoonomiaSpecies.rds"
  refrenceTreePrefix = "TrueCategoricalRefrenceTree"
   
  #MainTrees Location
  if(!is.na(cmdArgImport('m'))){
    mainTreesLocation = cmdArgImport('m')
  }else{
    message("No maintrees arg, using default")
  }
  #read in the tree based on filetype extension
  if(file_ext(mainTreesLocation) == "rds"){
    if(!exists("mainTrees")){mainTrees = readRDS(mainTreesLocation)}
  }else{
    if(!exists("mainTrees")){mainTrees = readTrees(mainTreesLocation)} 
  }
  
  
  #reference Tree
  if(!is.na(cmdArgImport('c'))){
    refrenceTreePrefix = cmdArgImport('c')
  }else{
    message("Using TrueCategoricalRefrenceTree as reference tree")
  }
  refrenceFolderName = paste("Output/",refrenceTreePrefix,"/", sep = "")
  
  #Substitution list 
  if(!is.null(cmdArgImport('u'))){
    substitutions = cmdArgImport('u')
  }else{
    message("No substitutions provided")
  }
  
}



# ---- Main code ----


# -- Load in the complex Tree --
refrenceTreeFilename = paste0(refrenceFolderName, refrenceTreePrefix, "CategoricalTree.rds")
refrenceTree = readRDS(refrenceTreeFilename)



refrencePhenotypeVectorFilename = paste0(refrenceFolderName, refrenceTreePrefix, "CategoricalPhenotypeVector.rds")
refrencePhenotypeVector = readRDS(refrencePhenotypeVectorFilename)

refrenceCategories = map_to_state_space(refrencePhenotypeVector)                              #and use it to connect branch lengths to phenotype name
refrenceCategoryNames = refrenceCategories$name2index                                         #store the length-phenotype connection


# --- Combining Categories --- 
phenotypeVector = refrencePhenotypeVector
phenotypeTree = refrenceTree



if(!is.null(substitutions)){
  edgePhenotypes = phenotypeTree$edge.length
  for(i in 1:length(refrenceCategoryNames)){
    edgePhenotypes[edgePhenotypes == i] = names(refrenceCategoryNames)[i]
  }
  for( i in 1:length(substitutions)){
    substitutePhenotypes = substitutions[[i]]
    
    # -- phenotypeVector -- 
    message(paste("replacing", substitutePhenotypes[1], "with", substitutePhenotypes[2]))
    phenotypeVector = gsub(substitutePhenotypes[1], substitutePhenotypes[2], phenotypeVector)
    edgePhenotypes = gsub(substitutePhenotypes[1], substitutePhenotypes[2], edgePhenotypes)
  }
  phenotypeCategories = map_to_state_space(phenotypeVector)                              #and use it to connect branch lengths to phenotype name
  phenotypeCategoryNames = phenotypeCategories$name2index                                         #store the length-phenotype connection
  phenotypeCategoryNames
  
  for(i in 1:length(phenotypeCategoryNames)){
    edgePhenotypes[edgePhenotypes == names(phenotypeCategoryNames)[i]] = i
  }
  phenotypeTree$edge.length = edgePhenotypes
  phenotypeTree$edge.length = as.integer(phenotypeTree$edge.length)
}


treeImageFilename = paste(outputFolderName, filePrefix, "CategoricalTree.pdf", sep="") #make a filename based on the prefix
pdf(treeImageFilename, height = length(phenotypeVector)/18)                     #make a pdf to store the plot, sized based on tree size
plotTreeCategorical(phenotypeTree, names(phenotypeCategoryNames))
dev.off()   

# -- handle paths -- 

refrencePathsFilename = paste0(refrenceFolderName, refrenceTreePrefix, "CategoricalPathsFile.rds")
refrencePaths = readRDS(refrencePathsFilename)
testPaths = tree2Paths(refrenceTree, mainTrees, categorical = T)

all.equal(testPaths, refrencePaths)
