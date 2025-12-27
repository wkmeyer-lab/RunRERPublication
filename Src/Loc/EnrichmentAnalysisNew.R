# -- Libraries 
clusterRun = F
clusterRun = T
if(clusterRun){.libPaths("/share/ceph/wym219group/shared/libraries/R4")} #add path to custom libraries to searched locations
library(RERconverge)
library("tools")
source("Src/Reu/cmdArgImport.R")
library(xlsx)


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
#'m=c("Data/tissue_specific.gmt", "Data/GSEA-c5-HsSymbols.gmt", "Data/EnrichmentHsSymbolsFile2.gmt")'
#'m=c("Data/DisGeNET.gmt", "Data/MGI_Mammalian_Phenotype_Level_4.gmt", "Data/GO_Biological_Process_2023.gmt")'
args = c('r=CategoricalDiet4Phen', "s=_Omnivore-Carnivore", 'm=c("Data/MGI_Mammalian_Phenotype_Level_4.gmt", "Data/GO_Biological_Process_2023.gmt", "Data/DisGeNET.gmt", "Data/tissue_specific.gmt", "Data/EnrichmentHsSymbolsFile2.gmt")', 'v=T', 'p=T') #This is a debug argument set. It is used to set arguments locally, when not running the code through a bash script.
args = c('r=LiverExpression',  'm=c("Data/MGI_Mammalian_Phenotype_Level_4.gmt", "Data/GO_Biological_Process_2023.gmt", "Data/DisGeNET.gmt", "Data/tissue_specific.gmt", "Data/EnrichmentHsSymbolsFile2.gmt")', 'v=T', 'p=F') #This is a debug argument set. It is used to set arguments locally, when not running the code through a bash script.
args = c('r=EcholocationUpdate',  'm=c("Data/MGI_Mammalian_Phenotype_Level_4.gmt", "Data/GO_Biological_Process_2023.gmt", "Data/DisGeNET.gmt", "Data/tissue_specific.gmt", "Data/EnrichmentHsSymbolsFile2.gmt")', 'v=T', 'p=T') #This is a debug argument set. It is used to set arguments locally, when not running the code through a bash script.
args = c('r=CategoricalDiet4Phen', 's=c("_Omnivore-Herbivore", "Carnivore-Herbivore", "_Omnivore-Insectivore", "Carnivore-Insectivore", "Herbivore-Insectivore", "_Omnivore-carnivore")', 'm=c("Data/MGI_Mammalian_Phenotype_Level_4.gmt", "Data/GO_Biological_Process_2023.gmt", "Data/DisGeNET.gmt", "Data/tissue_specific.gmt", "Data/EnrichmentHsSymbolsFile2.gmt")', 'v=T', 'p=F') #This is a debug argument set. It is used to set arguments locally, when not running the code through a bash script.
args = c('r=CategoricalDiet5Phen', 's=c("_Omnivore-Herbivore", "Carnivore-Herbivore", "_Omnivore-Insectivore", "Carnivore-Insectivore", "Herbivore-Insectivore", "_Omnivore-Piscivore", "Carnivore-Piscivore", "Herbivore-Piscivore", "Insectivore-Piscivore", "_Omnivore-carnivore")', 'm=c("Data/MGI_Mammalian_Phenotype_Level_4.gmt", "Data/GO_Biological_Process_2023.gmt", "Data/DisGeNET.gmt", "Data/tissue_specific.gmt", "Data/EnrichmentHsSymbolsFile2.gmt")', 'v=T', 'p=F') #This is a debug argument set. It is used to set arguments locally, when not running the code through a bash script.
args = c('r=CategoricalDiet3Phen', 'm=c("Data/MGI_Mammalian_Phenotype_Level_4.gmt", "Data/GO_Biological_Process_2023.gmt", "Data/DisGeNET.gmt", "Data/tissue_specific.gmt", "Data/EnrichmentHsSymbolsFile2.gmt")', 's=c("_Omnivore-Herbivore", "Carnivore-Herbivore", "_Omnivore-Carnivore", "Overall")', 'p=T')
args = c('r=LiverExpression3', 'm=c("Data/MGI_Mammalian_Phenotype_Level_4.gmt", "Data/GO_Biological_Process_2023.gmt", "Data/DisGeNET.gmt", "Data/tissue_specific.gmt", "Data/EnrichmentHsSymbolsFile2.gmt")', 'p=T', "c=LiverExpression3CorrelationDataPermulatedNamesConverted.rds")
args = c('r=EcholocationUpdate2',  'm=c("Data/MGI_Mammalian_Phenotype_Level_4.gmt", "Data/GO_Biological_Process_2023.gmt", "Data/DisGeNET.gmt", "Data/tissue_specific.gmt", "Data/EnrichmentHsSymbolsFile2.gmt")', 'v=T', 'p=T') #This is a debug argument set. It is used to set arguments locally, when not running the code through a bash script.
args = c('r=CVHRemake',  'm=c("Data/MGI_Mammalian_Phenotype_Level_4.gmt", "Data/GO_Biological_Process_2023.gmt", "Data/DisGeNET.gmt", "Data/tissue_specific.gmt", "Data/EnrichmentHsSymbolsFile2.gmt")', 'v=T', 'p=T') #This is a debug argument set. It is used to set arguments locally, when not running the code through a bash script.
args = c('r=CategoricalDiet3Phen', 'm=c("Data/MGI_Mammalian_Phenotype_Level_4.gmt", "Data/GO_Biological_Process_2023.gmt", "Data/DisGeNET.gmt", "Data/tissue_specific.gmt", "Data/EnrichmentHsSymbolsFile2.gmt")', 's=c("_Omnivore-Herbivore", "Carnivore-Herbivore", "_Omnivore-Carnivore", "Overall")', 'p=C')
args = c('r=CVHRemake',  'm=c("Data/customGeneSet.gmt")', 'v=T', 'p=T') #This is a debug argument set. It is used to set arguments locally, when not running the code through a bash script.
args = c('r=CategoricalDiet3Phen', 'm=c("Data/customGeneSet.gmt")', 's=c("_Omnivore-Herbivore", "Carnivore-Herbivore", "_Omnivore-Carnivore", "Overall")', 'p=C')
args = c('r=LiverExpression3', 'm=c("Data/customGeneSet.gmt")', 'p=T', "c=LiverExpression3CorrelationDataPermulatedNamesConverted.rds")


args = c('r=NewHiller4Phen', 'm=c("Data/MGI_Mammalian_Phenotype_Level_4.gmt", "Data/GO_Biological_Process_2023.gmt", "Data/DisGeNET.gmt", "Data/tissue_specific.gmt", "Data/EnrichmentHsSymbolsFile2.gmt")', 'p=F', 's=c("_Omnivore-Carnivore", "_Omnivore-Herbivore", "Carnivore-Herbivore", "_Omnivore-Insectivore", "Carnivore-Insectivore", "Herbivore-Insectivore", "Overall")')
args = c('r=NewHillerTestSupraPrimates', 'm=c("Data/MGI_Mammalian_Phenotype_Level_4.gmt", "Data/GO_Biological_Process_2023.gmt", "Data/DisGeNET.gmt", "Data/tissue_specific.gmt", "Data/EnrichmentHsSymbolsFile2.gmt")', 'p=F', 's=c("0-1", "Overall")')
args = c('r=NewHillerTestSupraPrimates', 'm=c("Data/MGI_Mammalian_Phenotype_Level_4.gmt", "Data/GO_Biological_Process_2023.gmt", "Data/DisGeNET.gmt", "Data/tissue_specific.gmt", "Data/EnrichmentHsSymbolsFile2.gmt")', 'p=F', 's=c("0-1", "Overall")')
args = c('r=NewHillerTestSupraPrimates', 'm=c("Data/MGI_Mammalian_Phenotype_Level_4.gmt", "Data/GO_Biological_Process_2023.gmt", "Data/DisGeNET.gmt", "Data/tissue_specific.gmt", "Data/EnrichmentHsSymbolsFile2.gmt")', 'p=F', 's=c("0-1", "Overall")')

args = c('r=MaturityLifespanPercent', 'm=c("Data/MGI_Mammalian_Phenotype_Level_4.gmt", "Data/GO_Biological_Process_2023.gmt", "Data/DisGeNET.gmt", "Data/tissue_specific.gmt", "Data/EnrichmentHsSymbolsFile2.gmt")', 'p=F')
args = c('r=MaturityLifespanPercent', 'm=c("Data/MGI_Mammalian_Phenotype_Level_4.gmt", "Data/GO_Biological_Process_2023.gmt", "Data/DisGeNET.gmt", "Data/tissue_specific.gmt", "Data/EnrichmentHsSymbolsFile2.gmt")', 'p=T', 'f=MaturityLifespanPercentPermulationPValue.rds')
args = c('r=MaturityLogRaw', 'm=c("Data/MGI_Mammalian_Phenotype_Level_4.gmt", "Data/GO_Biological_Process_2023.gmt", "Data/DisGeNET.gmt", "Data/tissue_specific.gmt", "Data/EnrichmentHsSymbolsFile2.gmt")', 'p=F')
args = c('r=PankajBodysize', 'm=c("Data/MGI_Mammalian_Phenotype_Level_4.gmt", "Data/GO_Biological_Process_2023.gmt", "Data/DisGeNET.gmt", "Data/tissue_specific.gmt", "Data/EnrichmentHsSymbolsFile2.gmt")', 'p=F')


args = c('r=MaturityLifespanPercent', 'm=c("Data/YifanGenesets.gmt")', 'p=T', 'f=MaturityLifespanPercentPermulationPValue.rds')
args = c('r=MaturityLogRaw', 'm=c("Data/YifanGenesets.gmt")', 'p=F')



args = c('r=CategoricalMobivoreTree', 'm=c("Data/MGI_Mammalian_Phenotype_Level_4.gmt", "Data/GO_Biological_Process_2023.gmt", "Data/DisGeNET.gmt", "Data/tissue_specific.gmt", "Data/EnrichmentHsSymbolsFile2.gmt")', 'p=F', 's=c("Omnivore-Vertivore", "Insectivore-Vertivore", "Herbivore-Vertivore", "Insectivore-Omnivore", "Herbivore-Omnivore", "Herbivore-Insectivore", "Overall")' )
args = c('r=CategoricalCarnivoreTree', 'm=c("Data/MGI_Mammalian_Phenotype_Level_4.gmt", "Data/GO_Biological_Process_2023.gmt", "Data/DisGeNET.gmt", "Data/tissue_specific.gmt", "Data/EnrichmentHsSymbolsFile2.gmt")', 'p=F', 's=c("Carnivore-Herbivore", "Carnivore-Omnivore", "Herbivore-Omnivore", "Overall")' )

args = c('r=CategoricalInsVertivoreTree', 'm=c("Data/MGI_Mammalian_Phenotype_Level_4.gmt", "Data/GO_Biological_Process_2023.gmt", "Data/DisGeNET.gmt", "Data/tissue_specific.gmt", "Data/EnrichmentHsSymbolsFile2.gmt")', 'p=F', 's=c("Omnivore-Vertivore", "Insectivore-Vertivore", "Herbivore-Vertivore", "Insectivore-Omnivore", "Herbivore-Omnivore", "Herbivore-Insectivore", "Overall")' )
args = c('r=CategoricalPrunedCarnivoreTree', 'm=c("Data/MGI_Mammalian_Phenotype_Level_4.gmt", "Data/GO_Biological_Process_2023.gmt", "Data/DisGeNET.gmt", "Data/tissue_specific.gmt", "Data/EnrichmentHsSymbolsFile2.gmt")', 'p=F', 's=c("Carnivore-Herbivore", "Carnivore-Omnivore", "Herbivore-Omnivore", "Overall")' )
args = c('r=CategoricalBinaryInsectivoreTree', 'm=c("Data/MGI_Mammalian_Phenotype_Level_4.gmt", "Data/GO_Biological_Process_2023.gmt", "Data/DisGeNET.gmt", "Data/tissue_specific.gmt", "Data/EnrichmentHsSymbolsFile2.gmt")', 'p=F', 's=c("Background-Insectivore", "Overall")' )
args = c('r=CategoricalBinaryVertivoreTree', 'm=c("Data/MGI_Mammalian_Phenotype_Level_4.gmt", "Data/GO_Biological_Process_2023.gmt", "Data/DisGeNET.gmt", "Data/tissue_specific.gmt", "Data/EnrichmentHsSymbolsFile2.gmt")', 'p=F', 's=c("Background-Vertivore", "Overall")' )
args = c('r=CategoricalBinaryHerbivoreTree', 'm=c("Data/MGI_Mammalian_Phenotype_Level_4.gmt", "Data/GO_Biological_Process_2023.gmt", "Data/DisGeNET.gmt", "Data/tissue_specific.gmt", "Data/EnrichmentHsSymbolsFile2.gmt")', 'p=F', 's=c("Background-Herbivore", "Overall")' )
args = c('r=CategoricalBinaryOmnivoreTree', 'm=c("Data/MGI_Mammalian_Phenotype_Level_4.gmt", "Data/GO_Biological_Process_2023.gmt", "Data/DisGeNET.gmt", "Data/tissue_specific.gmt", "Data/EnrichmentHsSymbolsFile2.gmt")', 'p=F', 's=c("Background-Omnivore", "Overall")' )


args = c('r=CategoricalInsVertivoreTree', 'm=c("Data/MGI_Mammalian_Phenotype_Level_4.gmt", "Data/GO_Biological_Process_2023.gmt", "Data/DisGeNET.gmt", "Data/tissue_specific.gmt", "Data/EnrichmentHsSymbolsFile2.gmt")', 'p=F', 's=c("Omnivore-Vertivore", "Insectivore-Vertivore", "Herbivore-Vertivore", "Insectivore-Omnivore", "Herbivore-Omnivore", "Herbivore-Insectivore", "Overall")' )

args = c('r=CategoricalInsVertivoreTree', 'm=c("Data/KeggReactome.gmt")', 'p=F', 's=c("Omnivore-Vertivore", "Insectivore-Vertivore", "Herbivore-Vertivore", "Insectivore-Omnivore", "Herbivore-Omnivore", "Herbivore-Insectivore", "Overall")' )
args = c('r=CategoricalPrunedCarnivoreTree', 'm=c("Data/KeggReactome.gmt")', 'p=F', 's=c("Carnivore-Herbivore", "Carnivore-Omnivore", "Herbivore-Omnivore", "Overall")' )

args = c('r=meanTemp', 'm=c("Data/MGI_Mammalian_Phenotype_Level_4.gmt", "Data/GO_Biological_Process_2023.gmt", "Data/DisGeNET.gmt", "Data/tissue_specific.gmt", "Data/EnrichmentHsSymbolsFile2.gmt", "Data/KeggReactome.gmt")', 'p=F')


args = c('r=CategoricalInsVertivoreTreeLiamInference', 'm=c("Data/MGI_Mammalian_Phenotype_Level_4.gmt", "Data/GO_Biological_Process_2023.gmt", "Data/DisGeNET.gmt", "Data/tissue_specific.gmt", "Data/EnrichmentHsSymbolsFile2.gmt")', 'p=F', 's=c("Omnivore-Vertivore", "Insectivore-Vertivore", "Herbivore-Vertivore", "Insectivore-Omnivore", "Herbivore-Omnivore", "Herbivore-Insectivore", "Overall")' )
args = c('r=CategoricalInsVertivoreTreeLiamInference', 'm=c("Data/KeggReactome.gmt")', 'p=F', 's=c("Omnivore-Vertivore", "Insectivore-Vertivore", "Herbivore-Vertivore", "Insectivore-Omnivore", "Herbivore-Omnivore", "Herbivore-Insectivore", "Overall")' )

args = c('r=CategoricalInsVertivoreTreeCarnivoreLiamInference', 'm=c("Data/KeggReactome.gmt")', 'p=F', 's=c("Carnivore-Omnivore", "Carnivore-Herbivore", "Herbivore-Omnivore", "Overall")' )
args = c('r=CategoricalInsVertivoreTreeCarnivoreLiamInference', 'm=c("Data/MGI_Mammalian_Phenotype_Level_4.gmt", "Data/GO_Biological_Process_2023.gmt", "Data/DisGeNET.gmt", "Data/tissue_specific.gmt", "Data/EnrichmentHsSymbolsFile2.gmt")', 'p=F', 's=c("Carnivore-Omnivore", "Carnivore-Herbivore", "Herbivore-Omnivore", "Overall")' )


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
usePermulations = TRUE
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
    
    if(length(subdirectoryValueList ==1)){outputFolderName = paste(outputFolderName, subdirectoryValueList[1], "/", sep=""); subdirectoryValue = subdirectoryValueList[1]}
    
  }else{
    message("No subdirectory specified.")
  }
}


#                   ------- Code Body -------- 
for(i in 1:length(subdirectoryValueList)){
    outputFolderName = paste("Output/",filePrefix,"/", sep = "")
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
      correlationData$P = correlationData$permP
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
  
  for(i in 1:length(gmtFileLocation)){
  #Load the gmt annotations 
  gmtAnnotations = read.gmt(gmtFileLocation[i])                                      #read the gmt file
  annotationsList = list(gmtAnnotations)                                          #reformat it into the format the next fuction expects
  enrichmentListName = substring(gmtFileLocation[i], 6, last = (nchar(gmtFileLocation[i]) - 4)) #make a geneset name based on the filename 
  names(annotationsList) = enrichmentListName                                     #name geneset list with that name 
  
  enrichmentResult = fastwilcoxGMTall(rerStats, annotationsList, outputGeneVals = T, num.g =2) #run enrichment analysis 
  
  #save the enrichment output
  enrichmentFileName = paste(outputFolderName, filePrefix, subdirectoryValue, "Enrichment-", enrichmentListName, ".rds", sep= "") #make a filename based on the prefix and geneset
  saveRDS(enrichmentResult, enrichmentFileName)                                   #Save the enrichment 
  enrichmentCsvName = enrichmentFileName = paste(outputFolderName, filePrefix, subdirectoryValue, "Enrichments.xlsx", sep= "") #make a filename based on the prefix and geneset
  gc()
  if(!i == 0){ #stops from writing two copies of the same sheet when no subdirectories. 
    write.xlsx(enrichmentResult, file=enrichmentCsvName, sheetName=enrichmentListName, row.names=T, append = T)
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
