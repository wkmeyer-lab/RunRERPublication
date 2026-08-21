# -- Libraries 
clusterRun = F
clusterRun = T
if(clusterRun){.libPaths("/share/ceph/wym219group/shared/libraries/R4")} #add path to custom libraries to searched locations

library(RERconverge)
library(tools)
source("Src/Reu/cmdArgImport.R")
library(data.table)

# -- Usage:
# This script can be used to run RER calculations and phenotype correlations with Binary, continuous, or categorical phenotypes.
# It inputs a phenotype tree made by the appropriate script for the style, and outputs an RER file, a Paths file, and a Correlation file. 

# -- Command arguments list
# r = filePrefix                                                               This is a prefix used to organize and separate files by analysis run. Always required. 
# v = <T or F>                                                                 This prefix is used to force the regeneration of the script's output, even if the files already exist. Not required, not always used.
# m = mainTreeFilename.txt or .rds                                             This sets the location of the maintrees file
# s = < ["b" or "binary"] or ["c" or "continuous"] or ["g" or "categorical"]>  This prefix is used to set the type of phenotype being supplied
# l = <min.sp value>                                                           This sets the min.sp value to be used in the correlation. 
# i = run instance                                                             This sets which alternate(s) to run on for parallelization. if left blank, will run on all alternates. 
#----------------

args = c('r=ComplexDietCentralAnalysis', 's=g', 'm=Data/zoonomiaAllMammalsTrees.rds', 'l=170', 'i=1')



# --- Standard start-up code ---
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

# --- Argument Imports ---
# Defaults
mainTreesLocation = "/share/ceph/wym219group/shared/projects/MammalDiet/Zoonomia/RemadeTreesAllZoonomiaSpecies.rds"  #this is the standard location on Sol 
binaryPhenotypeTreeLocation = NULL
speciesFilter = NULL
phenotypeStyle = "continuous"
continousMetric = "diff"
validMetrics = c("diff", "mean", "last")
minSpValue = 10
alternatesSpecified = F
usedAlternates = NULL

{ # Bracket used for collapsing purposes
  
  #MainTrees Location
  if(!is.na(cmdArgImport('m'))){
    mainTreesLocation = cmdArgImport('m')
  }else{
    message("No maintrees arg, using default")
  }
  
  #Phenotype Style
  if(!is.na(cmdArgImport('s'))){
    phenotypeStyle = cmdArgImport('s')
    #Convert the various input options; so that the style can be specified with several names 
    if(phenotypeStyle == "b" | phenotypeStyle == "B" | phenotypeStyle == "binary"){phenotypeStyle = "Binary"}
    if(phenotypeStyle == "c" | phenotypeStyle == "C" | phenotypeStyle == "continuous"){phenotypeStyle = "Continuous"}
    if(phenotypeStyle == "g" | phenotypeStyle == "G" | phenotypeStyle == "categorical"){phenotypeStyle = "Categorical"}
  }else{
    message("PhenotypeStyle not specified, using continuous")
  }
  
  #min.sp Value
  if(!is.na(cmdArgImport('l'))){
    minSpValue = cmdArgImport('l')
  }else{
    message("min.sp not specified, using 10")
  }
  
  #alternate set
  if(!is.na(cmdArgImport('i'))){
    usedAlternates = cmdArgImport('i')
    alternatesSpecified = T
  }else{
    message("instance not specified, running on all alternates.")
    
  }
}


# -- Read in the trees --
#MainTrees
if(file_ext(mainTreesLocation) == "rds"){
  if(!exists("mainTrees")){mainTrees = readRDS(mainTreesLocation)}
}else{
  if(!exists("mainTrees")){mainTrees = readTrees(mainTreesLocation)} 
}

# ---- Make required adjustments to run on alternates properly ---- 

alternateSets = readRDS(paste0(outputFolderName, filePrefix, "AlternatePruningSpecies.rds"))

outputFolderName = paste0(outputFolderName, "Alternates/")
primaryFilePrefix = filePrefix


#                   ------- Code Body -------- 

if(!alternatesSpecified){
  usedAlternates = 1:length(alternateSets)
}

for(i in usedAlternates){
  currentSet = alternateSets[[i]]
  alternateFilePrefix = paste0("/Alternates/Alternate", i)
  filePrefix = paste0("Alternate", i, primaryFilePrefix)


  #Phenotype tree
  phenotypeTreeFilename = paste(outputFolderName, filePrefix, phenotypeStyle, "Tree.rds", sep="")
  
  if(file_ext(phenotypeTreeFilename) == "rds"){                                   #if the tree is an RDS file
    phenotypeTree = readRDS(phenotypeTreeFilename)                                #Read as RDS
  }else{                                                                          #Otherwise
    phenotypeTree = readTrees(phenotypeTreeFilename)                              #read as text
  }
  
  speciesFilterFileName = paste(outputFolderName, filePrefix, "SpeciesFilter.rds",sep="") #Make the name of the location a pre-made filter would have to test for it
  speciesFilter = readRDS(speciesFilterFileName) 
  #Species filter
  if(all(is.null(speciesFilter))){ # If the species filter is meant to be empty, that is, all of the species should be used
    speciesFilter = mainTrees$masterTree$tip.label #include all of the species on the tree in the filter 
  }
  
  
  
  # --- RERs ---
  
  RERFileName = paste(outputFolderName, filePrefix, "RERFile.rds", sep= "")       #Set a filename for the RERs based on the prefix
  
  if(!file.exists(paste(RERFileName)) | forceUpdate){                             #if it does not exist, or update is forced 
    gc()
    RERObject = getAllResiduals(mainTrees, useSpecies = speciesFilter, plot = F)  #Calculate the RERs
    saveRDS(RERObject, file = RERFileName)                                        #Save them
  }else{                                                                          #Otherwise
    RERObject = readRDS(RERFileName)                                              #Use the existing ones
  }
  
  # --- PATHS ---
  
  pathsFileName = paste(outputFolderName, filePrefix, phenotypeStyle, "PathsFile.rds", sep= "") #Set a filename for the pathss based on the prefix and style
  if(!file.exists(paste(pathsFileName)) | forceUpdate){                           #if it does not exist, or update is forced 
    if(phenotypeStyle == "Binary"){                                               #If binary 
      pathsObject = tree2Paths(phenotypeTree, mainTrees, binarize=T, useSpecies = speciesFilter)  #make path with binary function
    }else if(phenotypeStyle == "Continuous"){                                      #If continous
      continousTraitVector = readRDS(phenotypeTreeLocation)                         #repurposing the phenotype tree argument for this, as it is the equivalent
      pathsObject = char2Paths(phenotypeTree, mainTrees, metric = continousMetric)
    } else if(phenotypeStyle == "Categorical"){                                   #if categorical
      pathsObject = tree2Paths(phenotypeTree, mainTrees, useSpecies = speciesFilter, categorical = TRUE) #do not binarize; the categorical data is already contained in the phenotype tree.
    }
    saveRDS(pathsObject, file = pathsFileName)                                    #Save the paths
  }else{
    pathsObject = readRDS(pathsFileName)                                          #If the file already exists, use the existing one.
  }
  
  
  # --- CORRELATION ---
  correlationFileName = paste(outputFolderName, filePrefix, "CorrelationFile", sep= "") #Make a correlation filename based on the prefix
  
  if(phenotypeStyle == "Categorical"){                                     #if categorical
    categoricalCorrelation = correlateWithCategoricalPhenotype(RERObject, pathsObject, min.sp = minSpValue, min.pos = 2) #Calculate with categorical, min 2 species per category 
    overalCategorical = categoricalCorrelation[[1]]                               #select the results relating to overall difference between all categories
    correlation = overalCategorical                                               # and classify it as the main correlation file
    
    #process the pairwise outputs
    pairwiseCategorical = categoricalCorrelation[[2]]                             #select the group of pairwise comparisons
    
    phenotypeVectorFilename = paste(outputFolderName, filePrefix, "CategoricalPhenotypeVector.rds",sep="") #select the phenotype vector based on prefix
    phenotypeVector = readRDS(phenotypeVectorFilename)                            #load in the phenotype vector 
    categories = map_to_state_space(phenotypeVector)                              #and use it to connect branch lengths to phenotype name
    categoryNames = categories$name2index                                         #store the length-phenotype connection
    
    pairwiseTableNames = names(pairwiseCategorical)                               #Prepare to repalce the number-number titles with phenotype-phenotype titles
    for(i in 1:length(categoryNames)){                                            #for each phenotype
      pairwiseTableNames= gsub(i, names(categoryNames)[i], pairwiseTableNames)                        #replace the number with the phenotype name  
    }
    names(pairwiseCategorical) = pairwiseTableNames                               #update the dataframe titles
    
    pairwiseCorrelationFileName = paste(outputFolderName, filePrefix, "PairwiseCorrelationFile", sep= "") #make a name for the pairwise comparisons based on prefix
    write.csv(pairwiseCategorical, file= paste(pairwiseCorrelationFileName, ".csv", sep=""), row.names = T, quote = F) #save the correlations as a csv
    saveRDS(pairwiseCategorical, paste(pairwiseCorrelationFileName, ".rds", sep="")) #and as an rds 
    
    combinedCategoricalCorrelationFilename = pairwiseCorrelationFileName = paste(outputFolderName, filePrefix, "CombinedCategoricalCorrelationFile", sep= "") # make this file for later functions that want it in combo
    saveRDS(categoricalCorrelation, paste(combinedCategoricalCorrelationFilename, ".rds", sep="")) #and as an rds 
    
    #save the outputs to subdirectories 
    outputSubdirectoryNoslash = paste(outputFolderName, "Overall", sep = "")
    if(!dir.exists(outputSubdirectoryNoslash)){                       #create that directory if it does not exist
      dir.create(outputSubdirectoryNoslash)
    }
    outputSubdirectory = paste(outputSubdirectoryNoslash, "/", sep="")
    
    correlationsOverallFilename = paste(outputSubdirectory, filePrefix, "OverallCorrelationFile.rds", sep= "")
    saveRDS(categoricalCorrelation[[1]], correlationsOverallFilename)
    
    for(i in 1:length(pairwiseTableNames)){
      pairwiseTableNames= gsub(" ", "", pairwiseTableNames)
      
      outputSubdirectoryNoslash = paste(outputFolderName, pairwiseTableNames[i], sep = "")
      if(!dir.exists(outputSubdirectoryNoslash)){                       #create that directory if it does not exist
        dir.create(outputSubdirectoryNoslash)
      }
      outputSubdirectory = paste(outputSubdirectoryNoslash, "/", sep="")
      
      correlationsPairFilename = paste(outputSubdirectory, filePrefix, pairwiseTableNames[i], "CorrelationFile",".rds", sep= "")
      saveRDS(categoricalCorrelation[[2]][[i]], correlationsPairFilename)
    }
    
  }
  
  write.csv(correlation, file= paste(correlationFileName, ".csv", sep=""), row.names = T, quote = F) #Save correlations as csv
  saveRDS(correlation, paste(correlationFileName, ".rds", sep=""))                          #and as an rds 

  
  
  
  
  
  # ------
  # -- Run the merged Nalysis for the alternate --
  # -----
  
  
  #-------
  # Make the merged phenotype
  #-------
  
  toMerge = c("Vertivore", "Insectivore")
  mergeValue = "Carnivore"
  
  
  #Handle the phenotyeVector
  mainPhenotypeVectorFilename = paste(outputFolderName, filePrefix, "CategoricalPhenotypeVector.rds",sep="") #make a filename based on the prefix
  
  mainPhenotypeVector = readRDS(mainPhenotypeVectorFilename)
  mainPhenotypeLabels = names(table(mainPhenotypeVector))
  
  replaceIndexes = which(mainPhenotypeLabels %in% toMerge)
  
  
  mergedPhenotypeLabels = mainPhenotypeLabels[-replaceIndexes]
  mergedPhenotypeLabels = append(mergedPhenotypeLabels, mergeValue)
  mergedPhenotypeLabels = mergedPhenotypeLabels[order(mergedPhenotypeLabels)]
  
  mergedPhenotypeVector = mainPhenotypeVector
  for(i in 1:length(toMerge)){
    mergedPhenotypeVector = gsub(toMerge[i], mergeValue, mergedPhenotypeVector)
  }
  
  mergedPhenotypeVectorFilename = paste(outputFolderName, filePrefix, "CategoricalMergedPhenotypeVector.rds",sep="") #make a filename based on the prefix
  saveRDS(mergedPhenotypeVector, mergedPhenotypeVectorFilename)
  
  #Handle the phenotype tree
  mainCategoricalTreeFilename = paste(outputFolderName, filePrefix, "CategoricalTree.rds", sep="") #make a filename based on the prefix
  mainPhenotypeTree = readRDS(mainCategoricalTreeFilename)
  
  bufferValue = 10
  mergedTree = mainPhenotypeTree
  for(i in 1:length(mainPhenotypeLabels)){
    if(mainPhenotypeLabels[i] %in% toMerge){
      mergedPosition = which(mergedPhenotypeLabels %in% mergeValue)
      mergedTree$edge.length[mergedTree$edge.length == i] = mergedPosition+bufferValue
    }else{
      currentLabel = mainPhenotypeLabels[i]
      mergedPosition = which(mergedPhenotypeLabels == currentLabel)
      mergedTree$edge.length[mergedTree$edge.length == i] = mergedPosition+bufferValue
    }
  }
  mergedTree$edge.length = mergedTree$edge.length-bufferValue
  
  mergedCategoricalTreeFilename = paste(outputFolderName, filePrefix, "MergedCategoricalTree.rds", sep="") #make a filename based on the prefix
  saveRDS(mergedTree, mergedCategoricalTreeFilename)
  
  #-------
  # Run RER and correlation for it 
  #-------
  
  #Phenotype tree
  phenotypeTree = mergedTree
  
  
  #Species filter
  speciesFilterFileName = paste(outputFolderName, filePrefix, "SpeciesFilter.rds",sep="") #Make the name of the location a pre-made filter would have to test for it
  speciesFilter = readRDS(speciesFilterFileName)                       #if so, use it 
  
  
  # --- RERs ---
  
  RERFileName = paste(outputFolderName, filePrefix, "RERFile.rds", sep= "")       #Set a filename for the RERs based on the prefix
  
  if(!file.exists(paste(RERFileName)) | forceUpdate){                             #if it does not exist, or update is forced 
    gc()
    RERObject = getAllResiduals(mainTrees, useSpecies = speciesFilter, plot = F)  #Calculate the RERs
    saveRDS(RERObject, file = RERFileName)                                        #Save them
  }else{                                                                          #Otherwise
    RERObject = readRDS(RERFileName)                                              #Use the existing ones
  }
  
  # --- PATHS ---
  
  pathsObject = tree2Paths(phenotypeTree, mainTrees, useSpecies = speciesFilter, categorical = TRUE) #do not binarize; the categorical data is already contained in the phenotype tree.
  
  
  
  # --- CORRELATION ---
  correlationFileName = paste(outputFolderName, filePrefix, "MergedCorrelationFile", sep= "") #Make a correlation filename based on the prefix
  
  {                                     #if categorical
    categoricalCorrelation = correlateWithCategoricalPhenotype(RERObject, pathsObject, min.sp = minSpValue, min.pos = 2) #Calculate with categorical, min 2 species per category 
    overalCategorical = categoricalCorrelation[[1]]                               #select the results relating to overall difference between all categories
    correlation = overalCategorical                                               # and classify it as the main correlation file
    
    #process the pairwise outputs
    pairwiseCategorical = categoricalCorrelation[[2]]                             #select the group of pairwise comparisons
    
    
    phenotypeVector = mergedPhenotypeVector                            #load in the phenotype vector 
    categories = map_to_state_space(phenotypeVector)                              #and use it to connect branch lengths to phenotype name
    categoryNames = categories$name2index                                         #store the length-phenotype connection
    
    pairwiseTableNames = names(pairwiseCategorical)                               #Prepare to repalce the number-number titles with phenotype-phenotype titles
    for(i in 1:length(categoryNames)){                                            #for each phenotype
      pairwiseTableNames= gsub(i, names(categoryNames)[i], pairwiseTableNames)                        #replace the number with the phenotype name  
    }
    names(pairwiseCategorical) = pairwiseTableNames                               #update the dataframe titles
    
    pairwiseCorrelationFileName = paste(outputFolderName, filePrefix, "MergedPairwiseCorrelationFile", sep= "") #make a name for the pairwise comparisons based on prefix
    write.csv(pairwiseCategorical, file= paste(pairwiseCorrelationFileName, ".csv", sep=""), row.names = T, quote = F) #save the correlations as a csv
    saveRDS(pairwiseCategorical, paste(pairwiseCorrelationFileName, ".rds", sep="")) #and as an rds 
    
    combinedCategoricalCorrelationFilename = pairwiseCorrelationFileName = paste(outputFolderName, filePrefix, "MergedCombinedCategoricalCorrelationFile", sep= "") # make this file for later functions that want it in combo
    saveRDS(categoricalCorrelation, paste(combinedCategoricalCorrelationFilename, ".rds", sep="")) #and as an rds 
    
    #save the outputs to subdirectories 
    outputSubdirectoryNoslash = paste(outputFolderName, "Overall", sep = "")
    if(!dir.exists(outputSubdirectoryNoslash)){                       #create that directory if it does not exist
      dir.create(outputSubdirectoryNoslash)
    }
    outputSubdirectory = paste(outputSubdirectoryNoslash, "/", sep="")
    
    correlationsOverallFilename = paste(outputSubdirectory, filePrefix, "MergedOverallCorrelationFile.rds", sep= "")
    saveRDS(categoricalCorrelation[[1]], correlationsOverallFilename)
    
    mergeUniquePairwiseNames = character()
    
    for(i in 1:length(pairwiseTableNames)){
      pairwiseTableNames= gsub(" ", "", pairwiseTableNames)
      
      outputSubdirectoryNoslash = paste(outputFolderName, pairwiseTableNames[i], sep = "")
      if(!dir.exists(outputSubdirectoryNoslash)){                       #create that directory if it does not exist
        dir.create(outputSubdirectoryNoslash)
      }
      if(length(grep(mergeValue, outputSubdirectoryNoslash)) != 0 ){ #check that this is a merge-unique directory
        outputSubdirectory = paste(outputSubdirectoryNoslash, "/", sep="")
        
        correlationsPairFilename = paste(outputSubdirectory, filePrefix, pairwiseTableNames[i], "CorrelationFile",".rds", sep= "")
        saveRDS(categoricalCorrelation[[2]][[i]], correlationsPairFilename)
        mergeUniquePairwiseNames = append(mergeUniquePairwiseNames, pairwiseTableNames[i])
      }
    }
    
  }
  
  write.csv(correlation, file= paste(correlationFileName, ".csv", sep=""), row.names = T, quote = F) #Save correlations as csv
  saveRDS(correlation, paste(correlationFileName, ".rds", sep=""))                          #and as an rds 
  
  
  #-------
  # Combine the existing Correlation file with the merged one 
  #-------
  
  mainPairwiseCorrelationFileName = paste(outputFolderName, filePrefix, "PairwiseCorrelationFile", sep= "") #make a name for the pairwise comparisons based on prefix
  mainPairwiseCategorical = readRDS(paste0(mainPairwiseCorrelationFileName, ".rds"))
  
  names(mainPairwiseCategorical) = gsub(" ", "", names(mainPairwiseCategorical))
  names(pairwiseCategorical) = gsub(" ", "", names(pairwiseCategorical))
  
  for(i in 1:length(mergeUniquePairwiseNames)){
    currentName = mergeUniquePairwiseNames[i]
    
    if(currentName %in% names(mainPairwiseCategorical)){
      currentIndex = which(names(mainPairwiseCategorical) %in% currentName)
      mainPairwiseCategorical[currentIndex] = pairwiseCategorical[which(names(pairwiseCategorical) %in% currentName)]
    }else{
      mainPairwiseCategorical = append(mainPairwiseCategorical, pairwiseCategorical[which(names(pairwiseCategorical) %in% currentName)])
    }
    pairwiseCategorical
  }
  
  
  
  
  
  write.csv(mainPairwiseCategorical, file= paste(mainPairwiseCorrelationFileName, ".csv", sep=""), row.names = T, quote = F) #save the correlations as a csv
  saveRDS(mainPairwiseCategorical, paste(mainPairwiseCorrelationFileName, ".rds", sep="")) #and as an rds 
  
  
}

