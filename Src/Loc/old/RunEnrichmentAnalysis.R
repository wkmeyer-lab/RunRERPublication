.libPaths("/share/ceph/wym219group/shared/libraries/R4") #add path to custom libraries to searched locations
library(RERconverge) #load RERconverge package
library(RERconverge)
library("tools")
source("Src/Reu/cmdArgImport.R")

#---- USAGE -----
#Used to generate GO Values from non-permulated RERs.  

# ARUGMENTS: 
#If an argument contains a '(' it is evaluated as code.
# r="filePrefix"                          This is the prefix attached to all files; a required argument. 
# m = gmtFileLocation.gmt                  This is the location of the main gmt file
# p = <T or F>                              This sets if the code should use permulated or unpermualted values
# f = "permulationPvalueFileLocation.rds"   This is a manual override to specify the script use a specific Permulation p-value file. 
  #If using any permulation p-value file other than "CombinedPrunedFastAll" with no run instance number, it must be specified manually.
#testing args: 
args = c('r=CVHRemake', 'p=F')
args = c('r=CVHRemake', 'm=Data/EnrichmentHsSymbolsFile.gmt')

#---- Prefix Setup -----
#------------------------
args = commandArgs(trailingOnly = TRUE)
{
  
  #- Import prefix -

  paste(args)
  message(args)
  
  #File Prefix
  if(!is.na(cmdArgImport('r'))){
    filePrefix = cmdArgImport('r')
  }else{
    stop("THIS IS AN ISSUE MESSAGE; SPECIFY FILE PREFIX")
  }
  #------------
  
  #- Make Output directory -
  
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
  
  #----------
}

{
  #----Startup----
  #----------------
  
  #--Default Values-- 
  gmtFileLocation = "Data/enrichmentGmtFile.gmt"
  usePermulations = TRUE
  usePermulationPValOverride = FALSE
  permulationPValOverride = NULL 
  
  #--Import arguments--
  #import gmt file location
  if(!is.na(cmdArgImport('m'))){
    gmtFileLocation = cmdArgImport('m')
  }else{
    paste("No gmt location arugment, using Data/enrichmentGmtFile.gmt")                          #Report using default
    message("No gmt location arugment, using Data/enrichmentGmtFile.gmt")
  }
  
  #Import permulation use
  if(!is.na(cmdArgImport('p'))){
    usePermulations = cmdArgImport('p')
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
      message("Using Manually specified permulation p-value file. NOTE this setting ignores the meta-combination argument.")
    }else{
      message("No Permulation pValue override, using standard CombinedPrunedFastAllPermulationsPValue.rds.")
    }
  }
  #--
}

#---- MAIN CODE ----
#-------------------

#Load correlation file
correlationFileLocation = paste(outputFolderName, filePrefix, "CorrelationFile.rds", sep= "")
correlationData = readRDS(correlationFileLocation)                            #Import the correlation data (non-permulated)
rerStats = getStat(correlationData)

if(usePermulations){
  if(usePermulationPValOverride){
    permulationFileLocation = permulationPValOverride
  }else{
    permulationFileLocation = paste(outputFolderName, filePrefix, "CombinedPrunedFastAllPermulationsPValue.rds", sep= "")
  }
  permulationValues = readRDS(permulationFileLocation)
  correlationData$P = permulationValues
}

rerStats = getStat(correlationData)

#Load the gmt annotations 
gmtAnnotations = read.gmt(gmtFileLocation)
annotationsList = list(gmtAnnotations)
enrichmentListName = substring(gmtFileLocation, 6, last = (nchar(gmtFileLocation) - 4))
names(annotationsList) = enrichmentListName

enrichmentResult = fastwilcoxGMTall(rerStats, annotationsList, outputGeneVals = T, num.g =10)

#save the enrichment output
enrichmentFileName = paste(outputFolderName, filePrefix, "EnrichmentFile.rds", sep= "")
saveRDS(enrichmentResult, enrichmentFileName)


# --- Visualize the enrichment ----
visualize = T
visualize = F
clean = T
clean = F
#enrichmentResult = readRDS(enrichmentFileName)
{
  #This is manual only -- run-as-script does not accept a visualize output because no way to output result. 
  #For a script version, use PvQvGoVisualize.R 
  
  
  if(visualize){
    library(stringr)
    library(insight)
    enrichmentResult2 = enrichmentResult[[1]]
    makeGOTable = function(data, collumn){
      ValueHead = head(data[order(collumn, decreasing = T),], n=40)
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

