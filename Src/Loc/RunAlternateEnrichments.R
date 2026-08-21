# -- Libraries 
clusterRun = F
clusterRun = T
if(clusterRun){.libPaths("/share/ceph/wym219group/shared/libraries/R4")} #add path to custom libraries to searched locations
library(RERconverge)
library("tools")
source("Src/Reu/cmdArgImport.R")
if(!clusterRun){library(xlsx)}


# -- Usage:
# This script is used to generate enrichment values from a correlation file. Can include permulation p values if provided. 

# -- Command arguments list
# r = filePrefix                            This is a prefix used to organize and separate files by analysis run. Always required. 
# v = <T or F>                              This prefix is used to force the regeneration of the script's output, even if the files already exist. Not required, not always used.
# m = gmtFileLocation.gmt                   This is the location of the main gmt file
# p = <T or F or C>                         This sets if the code should use permulated or unpermulated values. If set to "C", will use permulations for categorical values, which are stored in the main file.
# s = "subdirectoryName"                    This is used to specify a subdirectory for the analysis to be run in. Primarily used for the components of categorical results. 
# c = "correlationFileOverride.rds"         This can be used to manually set the correlation file. Used primarily to target pairwise comparisons of categorical phenotypes.
# f = "permulationPvalueFileLocation.rds"   This is a manual override to specify the script use a specific Permulation p-value file.
#If using any file other than "CombinedPrunedFastAll" with no run instance number, it must be specified manually.
#----------------


args = c('r=ComplexDietCentralAnalysis', 'm=c("Data/KeggReactome.gmt")', 'p=F', 's=c("Insectivore-Vertivore", "Herbivore-Vertivore", "Insectivore-Omnivore", "Herbivore-Omnivore", "Herbivore-Insectivore", "Overall")' )



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
# Defaults
gmtFileLocation = "Data/enrichmentGmtFile.gmt"
usePermulations = FALSE
usePermulationPValOverride = FALSE
permulationPValOverride = NULL 
useCorrelationOverride = FALSE
correlationOverride = NULL
useCategoricalPermulations = FALSE
useSubdirectory = FALSE
subdirectoryValueList = NULL

{ # Bracket used for collapsing purposes
  #Gmt file location
  if(!any(is.na(cmdArgImport('m')))){
    gmtFileLocation = cmdArgImport('m')
  }else{
    paste("No gmt location arugment, using Data/enrichmentGmtFile.gmt")                          #Report using default
    message("No gmt location arugment, using Data/enrichmentGmtFile.gmt")
  }
  
  #Permulation use
  if(!is.na(cmdArgImport('p'))){
    usePermulations = cmdArgImport('p')
    if(usePermulations == "C"){
      usePermulations = TRUE
      useCategoricalPermulations = TRUE
    }
    usePermulations = as.logical(usePermulations)
    if(is.na(usePermulations)){
      usePermulations = TRUE
      message("Use Permualtions value not interpretable as logical. Did you remember to capitalize? Using TRUE.")
    }
  }else{
    message("Use Permulations value not specified, using TRUE.")
  }
  if(usePermulations){
    #Import permulation p-value Override 
    if(!is.na(cmdArgImport('f'))){
      usePermulationPValOverride = TRUE
      permulationPValOverride = cmdArgImport('f')
      message("Using Manually specified permulation p-value file.")
    }else{
      message("No Permulation pValue override, using standard CombinedPrunedFastAllPermulationsPValue.rds.")
    }
  }
  #Import correlation file override
  if(!is.na(cmdArgImport('c'))){
    useCorrelationOverride = TRUE
    correlationOverride = cmdArgImport('c')
    message("Using Manually specified correlation file.")
  }else{
    message("No corelation override specified.")
  }
  
  #Import subdirectory
  if(!any(is.na(cmdArgImport('s')))){
    useSubdirectory = TRUE
    subdirectoryValueList = cmdArgImport('s')
    #message(paste("Using subdirectory", subdirectoryValue, "."))
    
    #if(length(subdirectoryValueList ==1)){
      #outputFolderName = paste(outputFolderName, subdirectoryValueList[1], "/", sep=""); 
      #subdirectoryValue = subdirectoryValueList[1]
      #}
    
  }else{
    message("No subdirectory specified.")
  }
}


#alternate set
if(!is.na(cmdArgImport('i'))){
  usedAlternates = cmdArgImport('i')
  alternatesSpecified = T
}else{
  message("instance not specified, running on all alternates.")
  
}


# ---- Make required adjustments to run on alternates properly ---- 

alternateSets = readRDS(paste0("Output/",filePrefix,"/", filePrefix, "AlternatePruningSpecies.rds"))

mainOutputFolderName = outputFolderName

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



  #                   ------- Code Body -------- 
  for(i in 1:length(subdirectoryValueList)){
    outputFolderName = paste("Output/",primaryFilePrefix,"/", "Alternates/", sep = "")
    message(paste("Using subdirectory", subdirectoryValueList[i], "."))
    if(useSubdirectory){
      outputFolderName = paste(outputFolderName, subdirectoryValueList[i], "/", sep="")
    }
    subdirectoryValue = subdirectoryValueList[i]
    
    #Load correlation file
    correlationFileLocation = paste(outputFolderName, filePrefix, subdirectoryValue, "CorrelationFile.rds", sep= "") #get the correlation file location based on prefix 
    if(useCategoricalPermulations){
      correlationFileLocation = paste(outputFolderName, filePrefix, subdirectoryValue, "PermulationsCorrelationFile", ".rds", sep= "")
    }
    if(useCorrelationOverride){                                                     #if a correlation override was specified, replace it with that
      correlationFileLocation = paste(outputFolderName, correlationOverride, sep="")                                      
    }
    correlationData = readRDS(correlationFileLocation)                              #Import the correlation data (non-permulated)
    
    if(usePermulations){                                                            #If permualtions are being used   
      if(useCategoricalPermulations){
        if(exists(correlationData$PermP.adj)){
          correlationData$P = correlationData$permP.adj
        }else{
          correlationData$P = correlationData$permP
        }
      }else{
        if(usePermulationPValOverride){                                               #check for a location override
          permulationFileLocation = paste(outputFolderName, permulationPValOverride, sep="")                           #if so, use it 
        }else{                                                                        #if not, use the default 
          permulationFileLocation = paste(outputFolderName, filePrefix, "CombinedPrunedFastAllPermulationsPValue.rds", sep= "") #get the default location based on prefix 
        }
        permulationValues = readRDS(permulationFileLocation)                          #read the permulation file
        if(ncol(permulationValues)>1){                                              # add handling for continuous permulations having two columns
          correlationData$P = permulationValues$permpval                                         #replace the P column in the correlation data with the permulation values. This is the only column that the later function checks. 
          #correlationData$Rho = permulationValues$permstats
        }else{
          correlationData$P = permulationValues                                         #replace the P column in the correlation data with the permulation values. This is the only column that the later function checks. 
        }
        comparePermulationsTable = cbind(readRDS(correlationFileLocation), permulationValues)
      }
    }
    
    rerStats = getStat(correlationData)                                             #processes the RERs somewhat into stat values. only uses the P column, and the sign of the Rho column. 
    
    if(usePermulations == F){
      enrichmentCsvName = paste(outputFolderName, filePrefix, subdirectoryValue, "Enrichments.xlsx", sep= "") #make a filename based on the prefix and geneset
    }else{
      enrichmentCsvName = paste(outputFolderName, filePrefix, subdirectoryValue, "Permulation-Enrichments.xlsx", sep= "") #make a filename based on the prefix and geneset
    }
    if(!clusterRun){file.remove(enrichmentCsvName)}
    for(j in 1:length(gmtFileLocation)){
      #Load the gmt annotations 
      gmtAnnotations = read.gmt(gmtFileLocation[j])                                      #read the gmt file
      annotationsList = list(gmtAnnotations)                                          #reformat it into the format the next fuction expects
      enrichmentListName = substring(gmtFileLocation[j], 6, last = (nchar(gmtFileLocation[j]) - 4)) #make a geneset name based on the filename 
      names(annotationsList) = enrichmentListName                                     #name geneset list with that name 
      
      enrichmentResult = fastwilcoxGMTall(rerStats, annotationsList, outputGeneVals = T, num.g =2) #run enrichment analysis 
      
      #save the enrichment output
      if(usePermulations == F){
        enrichmentFileName = paste(outputFolderName, filePrefix, subdirectoryValue, "Enrichment-", enrichmentListName, ".rds", sep= "") #make a filename based on the prefix and geneset
      }else{
        enrichmentFileName = paste(outputFolderName, filePrefix, subdirectoryValue, "Enrichment-Permulation-", enrichmentListName, ".rds", sep= "") #make a filename based on the prefix and geneset
      }
      saveRDS(enrichmentResult, enrichmentFileName)                                   #Save the enrichment 
      gc()
      if(!i == 0){ #stops from writing two copies of the same sheet when no subdirectories. 
        if(!clusterRun){write.xlsx(enrichmentResult[[1]], file=enrichmentCsvName, sheetName=enrichmentListName, row.names=T, append = T)}
      }
    }
  }
  # --- Visualize the enrichment ----
  
  {
    #This is manual only -- run-as-script does not accept a visualize output because no way to output result. 
    #For a script version, use PvQvGoVisualize.R 
    visualize = T
    visualize = F
    clean = T
    clean = F
    
    if(visualize){
      library(stringr)
      library(insight)
      enrichmentResult2 = enrichmentResult[[1]]
      makeGOTable = function(data, collumn){
        ValueHead = head(data[order(collumn, decreasing = T),], n=20)
        ValueHead$num.genes = as.character(ValueHead$num.genes)
        ValueHead$stat = round(ValueHead$stat, digits = 5)
        ValueHead$stat = as.character(ValueHead$stat)
        ValueHead = format_table(ValueHead, pretty_names = F, digits = "scientific5")
        ValueHead
      }
      enrichHead = makeGOTable(enrichmentResult2, abs(enrichmentResult2$stat))
      enrichHead
      textplot(enrichHead[1:4], mar = c(0,0,2,0), cmar = 1.5)
      if(usePermulations){
        title(main = paste("Top pathways by permulation"))
      }else{
        title(main = paste("Top pathways by non-permulation"))
      }
      
      if(clean){
        # ---- enrichment cleaning ----
        pathwaysToRemove = grep("CANCER", rownames(enrichHead))
        rowsToKeep = (!1:nrow(enrichHead) %in% pathwaysToRemove)
        cleanedHead = enrichHead[rowsToKeep,]
        textplot(cleanedHead[1:4], mar = c(0,0,2,0), cmar = 1.5)
        if(usePermulations){
          title(main = paste("Top pathways by permulation"))
        }else{
          title(main = paste("Top pathways by non-permulation"))
        }
        CleanheadFileName = paste(outputFolderName, filePrefix, "CleanedEnrichmentHead.csv", sep= "")
        write.csv(cleanedHead, CleanheadFileName)
      }
    }
  }  
}