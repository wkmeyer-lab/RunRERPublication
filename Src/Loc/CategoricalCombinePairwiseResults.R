#This script combines all of the pairwise comparisons, and the overall comparison, of a categorical RERConverge analysis into a single file. 
#For the older version of the script which included overlap, detla, and driver analyses, see CategoricalCombinePairwiseResults-OLD.R


clusterRun = F
clusterRun = T
if(clusterRun){.libPaths("/share/ceph/wym219group/shared/libraries/R4")} #add path to custom libraries to searched locations
library(RERconverge)
source("Src/Reu/cmdArgImport.R")
source("Src/Reu/trimGmtFilenameToSetname.R")

# -- Command arguments list
# r = filePrefix                            This is a prefix used to organize and separate files by analysis run. Always required. 
# v = <T or F>                              This prefix is used to force the regeneration of the script's output, even if the files already exist. Not required, not always used.
# s = pairwise Sets                         This sets the subdirectories that the code should be run using. 
# g = gene set                              This sets the gene set (set name or file name) that the script should be run on, or "gene" for use on genes. 
# c = significance cutoff                   This sets the cutoff for inclusion in the significant columns
# w = writeData                             This sets if the data should be saved, potentially overwriting old combinations
# p = usePermulations                       This sets if permulation data should be used. If set to "A", it will add permulations data, but use the unpermulationed data for significance. 


args = c('r=ComplexDietCentralAnalysisSimplify', 'g=c("gene", "KeggReactome")')

{
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
}


# --- Argument Imports ---
# Defaults
significanceCutoff = 0.05
pairwiseSets = c("Herbivore-Insectivore", "Herbivore-Vertivore", "Carnivore-Herbivore", "Herbivore-Omnivore", "Insectivore-Vertivore", "Omnivore-Vertivore", "Invertivore-Omnivore")
pairwiseSets = c("Herbivore-Invertivore", "Herbivore-Vertivore", "Carnivore-Herbivore", "Herbivore-Omnivore", "Invertivore-Vertivore", "Omnivore-Vertivore", "Invertivore-Omnivore")
geneSet = NULL
usingGo = F
usingGene = T
validGeneArguments = c("NULL", "gene", "Gene")

writeData = T

usePermulations = F
addPermulations = F

{ # Bracket used for collapsing purposes
  
  #pairwise Sets
  if(!all(is.na(cmdArgImport('s')))){
    if(all(cmdArgImport('s') == "NULL")){
      folderNames = basename(list.dirs(outputFolderName, full.names = TRUE, recursive = FALSE))
      pairwiseSets = folderNames[grep("-", folderNames)]
    }else{
      pairwiseSets = cmdArgImport('s') 
    }
  }else{
    folderNames = basename(list.dirs(outputFolderName, full.names = TRUE, recursive = FALSE))
    pairwiseSets = folderNames[grep("-", folderNames)]
    message("Folders not specified, using all folders. Note that if driver analysis is active, this may cause issues.")
  }

  #geneset
  if(!all(is.na(cmdArgImport('g')))){
    geneSets = cmdArgImport('g')
    if(is.null(geneSets) | any(geneSets %in% validGeneArguments)){
      usingGene = T
      message("Null geneset specified or gene specified, performing gene correlation.")
    }else{
      usingGene = F
    }
    
    if(any(!geneSet %in% validGeneArguments)){
      usingGo = T
      geneSet = trimGmtFilenameToSetname(geneSet)
    }else{
      usingGo = F
    }
    
  }else{
    message("No geneset specified or gene specified, performing gene correlation.")
  }
  
  #cutoff 
  if(!is.na(cmdArgImport('c'))){
    significanceCutoff = cmdArgImport('c')
  }else{
    message("pvalue cuttoff not specified, using 0.05")
  }  
  
  
  #writeData 
  if(!is.na(cmdArgImport('w'))){
    writeData = as.logical(cmdArgImport('w'))
  }else{
    message("writeData not specified, using TRUE")
  }
  
  #Use Permulations 
  if(!is.na(cmdArgImport('p'))){
    if(cmdArgImport('p') == "A"){
      usePermulations = T
      addPermulations = T
    }else if(cmdArgImport('p') == "C"){
      usePermulations = T
    }else{
      usePermulations = as.logical(cmdArgImport('p'))
    }

  }else{
    message("use Permulations not specified, using FLASE")
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

if(!usePermulations){
  pairwiseCorrelationFileName = paste(outputFolderName, filePrefix, "PairwiseCorrelationFile.rds", sep= "") #make a name for the pairwise comparisons based on prefix
  overallCorrelationFilename = paste(outputFolderName, filePrefix, "CorrelationFile.rds", sep= "") #make a name for the pairwise comparisons based on prefix
}else{
  pairwiseCorrelationFileName = paste(outputFolderName, filePrefix, "PermulationPairwiseCorrelationFile.rds", sep= "") #make a name for the pairwise comparisons based on prefix
  overallCorrelationFilename = paste(outputFolderName, filePrefix, "PermulationCorrelationFile.rds", sep= "") #make a name for the pairwise comparisons based on prefix
}

correlationResults = readRDS(pairwiseCorrelationFileName)
overallResults = readRDS(overallCorrelationFilename)
overall = list(overallResults)
names(overall)[[1]] = "Overall"
correlationResults = append(correlationResults, overall)
pairwiseSets = append(pairwiseSets, "Overall")

if(usingGene){
  combinedResults = NA
  combinedDrivers = NA
  combinedBinaries = NA
  prefixSet = NULL
  prefixList = NULL
  
  { # load in the data 
    for(i in 1:length(pairwiseSets)){
      currentSet = pairwiseSets[i]
      if(grepl(" ", currentSet)){
        correlationSubsetName = gsub("-", " - ", currentSet)
      }else{
        correlationSubsetName = currentSet
      }
      
      correlationPrefix = paste(substr(strsplit(currentSet, split = "-")[[1]],1,1), collapse = '')
      prefixEntry = correlationPrefix; names(prefixEntry) = currentSet; prefixSet = append(prefixSet, prefixEntry)
      for(j in 1:2){
        prefixSingle = strsplit(correlationPrefix, split = "")[[1]][j]; names(prefixSingle) = strsplit(currentSet, split="-")[[1]][j]; prefixList = append(prefixList, prefixSingle)
      }
      currentDataframe = which(names(correlationResults) == correlationSubsetName) 
      
      currentResults = correlationResults[[currentDataframe]]
      if(usePermulations){
        if(exists(currentResults$permP.adj)){
          currentResults$permSignificant = currentResults$permP.adj < significanceCutoff
          currentResults$unpermSignificant = currentResults$p.adj < significanceCutoff
        }else{
          currentResults$permSignificant = currentResults$permP < significanceCutoff
          currentResults$unpermSignificant = currentResults$p.adj < significanceCutoff
        }
        if(addPermulations){
          currentResults$significant = currentResults$unpermSignificant
        }else{
          currentResults$significant = currentResults$permSignificant
        }
      }else{
        currentResults$significant = currentResults$p.adj < significanceCutoff
      }
      
      
      names(currentResults) = paste0(correlationPrefix, "-", names(currentResults))
      combinedResults = cbind(combinedResults, currentResults)
      
      rm(currentResults)
      
    }
  }
  
  combinedResults = combinedResults[,-1]
  
  # -- save combination -- 
  if(writeData){
    combinedDataFilename = paste0(outputFolderName, filePrefix, "combinedGeneResults")
    write.csv(combinedResults, paste0(combinedDataFilename, ".csv"))
    saveRDS(combinedResults, paste0(combinedDataFilename, ".rds"))
  }
}


# --- Import GO Data --- 
if(usingGo){
  goSets = geneSets[-which(geneSets %in% validGeneArguments)]
  for(j in goSets){
    geneSet = j
    
    GOResults = list()
    for(i in 1:length(pairwiseSets)){
      currentSet = pairwiseSets[i]
      correlationSubsetName = gsub("-", " - ", currentSet)
      correlationPrefix = paste(substr(strsplit(currentSet, split = "-")[[1]],1,1), collapse = '')
      
      goFilename = paste0(outputFolderName, currentSet, "/", filePrefix, currentSet,"Enrichment-", geneSet, ".rds")
      tryCatch({
        currentGoData = readRDS(goFilename)[[1]]

        if(usePermulations){
          if(exists(currentGoData$permP.adj)){
            currentGoData$permSignificant = currentGoData$permP.adj < significanceCutoff
            currentGoData$unpermSignificant = currentGoData$p.adj < significanceCutoff        
          }else{
            currentGoData$permSignificant = currentGoData$permP < significanceCutoff
            currentGoData$unpermSignificant = currentGoData$p.adj < significanceCutoff
          }
          if(addPermulations){
            currentGoData$significant = currentGoData$unpermSignificant
          }else{
            currentGoData$significant = currentGoData$permSignificant
          }
        }else{
          currentGoData$significant = currentGoData$p.adj < significanceCutoff
        }
        
        
        
        names(currentGoData) = paste0(correlationPrefix, "-", names(currentGoData))
        
        
        GOResults[[i]] = currentGoData
        names(GOResults)[i] = correlationPrefix
      }, error = function(e) {
        message("Skipping file: ", goFilename, " | Error: ", conditionMessage(e))
      })
    }
    rm(currentGoData)
    missingGoResults = which(sapply(GOResults, is.null)) # Remove any empty lists generated by missing enrichment files 
    GOResults = GOResults[-missingGoResults]
    
    
    GoCombinedResults = NA
    longestGOResult = which(sapply(GOResults, nrow) == max(sapply(GOResults, nrow)))[1]
    for(i in 1:length(GOResults)){
      if(all(rownames(GOResults[[longestGOResult]]) == rownames(GOResults[[i]]))){
        cat("Combining GO Data", i, "\n")
        GoCombinedResults = cbind(GoCombinedResults, GOResults[[i]])
      }else{
        tryCatch({
          cat("Mismatch in order; reordering current data to first data; Combining GO Data", i, "\n")
          
          which(!rownames(GOResults[[longestGOResult]]) %in% rownames(GOResults[[i]]))
          
          matchedOrder = match(rownames(GOResults[[longestGOResult]]), rownames(GOResults[[i]]))
          #matchedOrder = matchedOrder[-which(is.na(matchedOrder))]
          GOResults[[i]] = GOResults[[i]][matchedOrder,]
          GoCombinedResults = cbind(GoCombinedResults, GOResults[[i]])
        },error = function(e){      
          stop("Rownames of GO results are not the same and cannot be reordered, there is an issue with the GO data.")
        }
        
        )
      }
    }
    GoCombinedResults = GoCombinedResults[,-1]
    rm(GOResults)
    
    
    
    if(writeData){
      combinedGODataFilename = paste0(outputFolderName, filePrefix, "combinedGOResults-", geneSet)
      write.csv(GoCombinedResults, paste0(combinedGODataFilename, ".csv"))
      saveRDS(GoCombinedResults, paste0(combinedGODataFilename, ".rds"))
    }
  }
}

