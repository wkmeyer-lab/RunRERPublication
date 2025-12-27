clusterRun = F
clusterRun = T
library(RERconverge)
library(ggvenn)
library(stats)
library(combinat)
library(ggpointdensity)
library(viridis)
library(ggpmisc)
library(eulerr)
library(gridExtra)
library(grid)
library(gridGraphics)
source("Src/Reu/cmdArgImport.R")

args = c("r=CategoricalInsvertivoreTree", "p=NULL", "g=gene")
args = c("r=CategoricalInsvertivoreTree", "p=NULL", "g=KeggReactome")

args = c("r=CategoricalInsvertivoreTreeLiamInference", "p=NULL", "g=gene")
args = c("r=CategoricalInsvertivoreTreeLiamInference", "p=NULL", "g=KeggReactome")

args = c("r=CategoricalInsvertivoreTreeCarnivoreLiamInference", "p=NULL", "g=gene", "c=0.05", "s=F")
args = c("r=CategoricalInsvertivoreTreeCarnivoreLiamInference", "p=NULL", "g=KeggReactome", "c=0.1", "s=T")

args = c("r=CategoricalInsvertivoreTreeLiamInference", "p=NULL", "g=MGI_Mammalian_Phenotype_Level_4", "c=0.05", "s=T")
args = c("r=CategoricalInsvertivoreTreeLiamInference", "p=NULL", "g=GO_Biological_Process_2023", "c=0.05", "s=T")
args = c("r=CategoricalInsvertivoreTreeLiamInference", "p=NULL", "g=DisGeNET", "c=0.05", "s=T")
args = c("r=CategoricalInsvertivoreTreeLiamInference", "p=NULL", "g=tissue_specific", "c=0.05", "s=T")
args = c("r=CategoricalInsvertivoreTreeLiamInference", "p=NULL", "g=EnrichmentHsSymbolsFile2", "c=0.05", "s=T")




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



# --- Argument Imports ---
# Defaults
significanceCutoff = 0.05
pairwiseSets = c("Herbivore-Insectivore", "Herbivore-Vertivore", "Carnivore-Herbivore", "Herbivore-Omnivore", "Insectivore-Vertivore", "Omnivore-Vertivore", "Invertivore-Omnivore")
geneSet = NULL
usingGo = F
saveData = T
usingGene = T
saveCombinedData = T

{ # Bracket used for collapsing purposes
  
  #pairwise Sets
  if(!is.na(cmdArgImport('p'))){
    pairwiseSets = cmdArgImport('p')
  }else{
    folderNames = basename(list.dirs(outputFolderName, full.names = TRUE, recursive = FALSE))
    pairwiseSets = folderNames[grep("-", folderNames)]
    message("Folders not specified, using all folders. Note that if driver analysis is active, this may cause issues.")
  }
  if(is.null(pairwiseSets) | pairwiseSets == "NULL"){
    folderNames = basename(list.dirs(outputFolderName, full.names = TRUE, recursive = FALSE))
    pairwiseSets = folderNames[grep("-", folderNames)]
    message("Null Folders specified, using all folders.")
  }
  
  #geneset
  if(!is.na(cmdArgImport('g'))){
    geneSet = cmdArgImport('g')
    if(is.null(geneSet) | geneSet == "NULL" | geneSet == "gene" | geneSet == "Gene"){
      usingGene = T
      message("Null geneset specified or gene specified, performing gene correlation.")
    }else{
      usingGo = T
      usingGene = F
    }
  }else{
    message("No geneset specified or gene specified, performing gene correlation.")
  }
  
  #cutoffsaveData 
  if(!is.na(cmdArgImport('c'))){
    saveData = cmdArgImport('c')
  }else{
    message("pvalue cuttoff not specified, using 0.05")
  }  
  
  
  #saveData 
  if(!is.na(cmdArgImport('s'))){
    saveData = as.logical(cmdArgImport('s'))
  }else{
    message("saveData not specified, using TRUE")
  }
}

# ---- Main Code ----- 



# -- Make central RERData object 

getComparisionDifference = function(dataframe, colOne, colTwo){
  colOneIndex = names(dataframe)[which(names(dataframe) == colOne)]
  colTwoIndex = names(dataframe)[which(names(dataframe) == colTwo)]
  distanceFromEqual = abs(dataframe[colOneIndex] - dataframe[colTwoIndex]) / sqrt(2)
  distanceFromEqual
}

pairwiseCorrelationFileName = paste(outputFolderName, filePrefix, "PairwiseCorrelationFile.rds", sep= "") #make a name for the pairwise comparisons based on prefix
correlationResults = readRDS(pairwiseCorrelationFileName)

if(usingGene){
  combinedResults = NA
  combinedDrivers = NA
  combinedBinaries = NA
  prefixSet = NULL
  prefixList = NULL
  
  { # load in the data 
    for(i in 1:length(pairwiseSets)){
      currentSet = pairwiseSets[i]
      correlationSubsetName = gsub("-", " - ", currentSet)
      correlationPrefix = paste(substr(strsplit(currentSet, split = "-")[[1]],1,1), collapse = '')
      prefixEntry = correlationPrefix; names(prefixEntry) = currentSet; prefixSet = append(prefixSet, prefixEntry)
      for(j in 1:2){
        prefixSingle = strsplit(correlationPrefix, split = "")[[1]][j]; names(prefixSingle) = strsplit(currentSet, split="-")[[1]][j]; prefixList = append(prefixList, prefixSingle)
      }
      currentDataframe = which(names(correlationResults) == correlationSubsetName) 
      
      currentResults = correlationResults[[currentDataframe]]
      currentResults$significant = currentResults$p.adj < significanceCutoff
      
      
      names(currentResults) = paste0(correlationPrefix, "-", names(currentResults))
      combinedResults = cbind(combinedResults, currentResults)
      
      driverFilename = paste0(outputFolderName, currentSet, "/", filePrefix, currentSet, "DriverTable.rds")
      if(file.exists(driverFilename)){
        driverTable = readRDS(driverFilename)
        driverTable = driverTable[,-grep("main", names(driverTable))]
        names(driverTable)[which(names(driverTable) == "Driver")] = paste0(correlationPrefix, "-", "Driver")
        names(driverTable)[which(names(driverTable) == "DriverNumeric")] = paste0(correlationPrefix, "-", "DriverNumeric")    
        
        combinedDrivers = cbind(combinedDrivers, driverTable[,grep(paste0(correlationPrefix,"-"), names(driverTable))])
        combinedBinaries = cbind(combinedBinaries, driverTable[,-grep(paste0(correlationPrefix,"-"), names(driverTable))])
        rm(driverTable)
      }else{
        message(paste0("Driver file not not exist for ", currentSet, ". Not including driver and continuing."))
      }
      rm(currentResults)
      
    }
  }
  
  combinedResults = combinedResults[,-1]
  if(!is.na(combinedDrivers)){
    combinedDrivers = combinedDrivers[,-1]
    combinedBinaries = combinedBinaries[,-1]
    combinedBinaries = combinedBinaries[,-grep(".1", names(combinedBinaries))]
  
    combinedResults = cbind(combinedResults, combinedDrivers)
    combinedResults = cbind(combinedResults, combinedBinaries)
  }
  rm(combinedDrivers); rm(combinedBinaries)
  
  
  
  # -- Add overlap information -- 
  significanceColumns = names(combinedResults)[grep("significant", names(combinedResults))]
  geneSignificanceResults = combinedResults[, names(combinedResults) %in% significanceColumns]
  
  for(i in 2:length(significanceColumns)){
    combinations = combn(significanceColumns, i, simplify = FALSE)
    for(j in 1:length(combinations)){
      currentCombination = combinations[[j]]
      headers = gsub("-.*","",  currentCombination)
      comboName = paste0(paste0(headers, collapse = "-"), "-Overlap")
      
      colsToCompare = geneSignificanceResults[,names(geneSignificanceResults) %in% currentCombination]
      
      comboValue = apply(colsToCompare, 1, function(row) all(row == TRUE) == 1)
      which(comboValue)
      combinedResults$newOverlapColumn = comboValue
      names(combinedResults)[length(names(combinedResults))] = comboName
      rm(colsToCompare)
    }
    rm(combinations)
  }
  
  # -- Add Delta information -- 
  rhoColumns = names(combinedResults)[grep("-Rho", names(combinedResults))]
  geneRhoColumns = combinedResults[, names(combinedResults) %in% rhoColumns]
  for(i in 2:length(rhoColumns)){
    combinations = combn(rhoColumns, i, simplify = FALSE)
    for(j in 1:length(combinations)){
      currentCombination = combinations[[j]]
      if(length(currentCombination) > 2){next}
      headers = gsub("-.*","",  currentCombination)
      comboName = paste0(paste0(headers, collapse = "-"), "-Delta")
      
      colsToCompare = geneRhoColumns[,names(geneRhoColumns) %in% currentCombination]
      
      deltaValue = abs(colsToCompare[1] - colsToCompare[2]) / sqrt(2)
      
      
      
      combinedResults$newDeltaColumn = deltaValue[,1]
      names(combinedResults)[length(names(combinedResults))] = comboName
      rm(colsToCompare)
    }
    rm(combinations)
  }
  
  # -- save combination -- 
  if(saveData){
    combinedDataFilename = paste0(outputFolderName, filePrefix, "combinedGeneResults")
    write.csv(combinedResults, paste0(combinedDataFilename, ".csv"))
    saveRDS(combinedResults, paste0(combinedDataFilename, ".rds"))
  }
}


# --- Import GO Data --- 
if(usingGo){
  GOResults = list()
  for(i in 1:length(pairwiseSets)){
    currentSet = pairwiseSets[i]
    correlationSubsetName = gsub("-", " - ", currentSet)
    correlationPrefix = paste(substr(strsplit(currentSet, split = "-")[[1]],1,1), collapse = '')
    
    goFilename = paste0(outputFolderName, currentSet, "/", filePrefix, currentSet,"Enrichment-", geneSet, ".rds")
    currentGoData = readRDS(goFilename)[[1]]
    currentGoData$significant = currentGoData$p.adj < significanceCutoff
    
    
    names(currentGoData) = paste0(correlationPrefix, "-", names(currentGoData))
    
    
    driverFilename = paste0(outputFolderName, currentSet, "/", filePrefix, currentSet, "GoDriverTable-", geneSet, ".rds")
    if(file.exists(driverFilename)){
      driverTable = readRDS(driverFilename)
      if(all(rownames(driverTable) == rownames(currentGoData))){
        currentGoData$Driver = driverTable[which(names(driverTable) == "Driver")][[1]]
        currentGoData$DriverNumeric = driverTable[which(names(driverTable) == "DriverNumeric")][[1]]
        
        names(currentGoData)[which(names(currentGoData) == "Driver")] = paste0(correlationPrefix, "-", "Driver")
        names(currentGoData)[which(names(currentGoData) == "DriverNumeric")] = paste0(correlationPrefix, "-", "DriverNumeric")    
        
      }
      rm(driverTable)
    }
    GOResults[[i]] = currentGoData
    names(GOResults)[i] = correlationPrefix
    
  }
  rm(currentGoData)
  
  GoCombinedResults = NA
  for(i in 1:length(GOResults)){
    if(all(rownames(GOResults[[1]]) == rownames(GOResults[[i]]))){
      cat("Combining GO Data", i, "\n")
      GoCombinedResults = cbind(GoCombinedResults, GOResults[[i]])
    }else{
      stop("Rownames of GO results are not the same, there is an issue with the GO data.")
    }
  }
  GoCombinedResults = GoCombinedResults[,-1]
  rm(GOResults)
  
  # -- Add overlap information -- 
  GoSignificanceColumns = names(GoCombinedResults)[grep("significant", names(GoCombinedResults))]
  GoSignificanceResults = GoCombinedResults[, names(GoCombinedResults) %in% GoSignificanceColumns]
  
  for(i in 2:length(GoSignificanceColumns)){
    combinations = combn(GoSignificanceColumns, i, simplify = FALSE)
    for(j in 1:length(combinations)){
      currentCombination = combinations[[j]]
      headers = gsub("-.*","",  currentCombination)
      comboName = paste0(paste0(headers, collapse = "-"), "-Overlap")
      
      colsToCompare = GoSignificanceResults[,names(GoSignificanceResults) %in% currentCombination]
      
      comboValue = apply(colsToCompare, 1, function(row) all(row == TRUE) == 1)
      which(comboValue)
      GoCombinedResults$newOverlapColumn = comboValue
      names(GoCombinedResults)[length(names(GoCombinedResults))] = comboName
      rm(colsToCompare)
    }
    rm(combinations)
  }
  
  if(saveCombinedData){
    combinedGODataFilename = paste0(outputFolderName, filePrefix, "combinedGOResults-", geneSet)
    write.csv(GoCombinedResults, paste0(combinedGODataFilename, ".csv"))
    if(saveData){saveRDS(GoCombinedResults, paste0(combinedGODataFilename, ".rds"))}
  }
}

