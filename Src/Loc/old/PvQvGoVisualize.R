.libPaths("/share/ceph/wym219group/shared/libraries/R4") #add path to custom libraries to searched locations
library(RERconverge)
library("tools")
source("Src/Reu/cmdArgImport.R")
library(qvalue)
library(insight)
library(stringr)


#---- USAGE -----
#Used to generate various plots for data analysis. 
#Makes: 
  #P value table; Q value table; Histograms of p-values for positive Rho entires, histogram of p-values for negative Rho entires. 

# ARUGMENTS: 
#If an argument contains a '(' it is evaluated as code.
# 'r="filePrefix"'                          This is the prefix attached to all files; a required argument. 
# g = <T OR F>                              This sets if Gene Ontology analysis is run. Default TRUE. 
# p = <T or F>                              This sets if the code should use permulated or unpermualted values, or displays both 
# f = "permulationPvalueFileLocation.rds"   This is a manual override to specify the script use a specific Permulation p-value file. 
      #If using any permulation p-value file other than "CombinedPrunedFastAll" with no run instance number, it must be specified manually.
#testing args: 
args = c('r=CVHRemake', 'g=F')
args = c('r=Domestication', 'g=F', 'p=F')
args = c('r=CategoricalDiet', 'g=F', 'p=F')
{
  #---- Initial Setup -----
  #------------------------
  
  # --- Import prefix ----
  args = commandArgs(trailingOnly = TRUE)
  paste(args)
  message(args)
  
  #File Prefix
  if(!is.na(cmdArgImport('r'))){
    filePrefix = cmdArgImport('r')
  }else{
    stop("THIS IS AN ISSUE MESSAGE; SPECIFY FILE PREFIX")
  }
  
  #------------
  
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
  
  #----------
}

{
  #----Startup----
  #----------------
  
  #--Default Values-- 
  performGeneOntolgy = TRUE
  usePermulationPValOverride = FALSE
  permulationPValOverride = NULL 
  metacombineValue = FALSE
  usePermulations = TRUE
  enrichmentAnnotationListName = "MSigDBPathways"
  
  #--Import arguments-- 
  
  #Import if should perform gene ontology 
  if(!is.na(cmdArgImport('g'))){
    performGeneOntolgy = cmdArgImport('g')
    performGeneOntolgy = as.logical(performGeneOntolgy)
    if(is.na(performGeneOntolgy)){
      performGeneOntolgy = TRUE
      message("Include geneset enrichment value not interpretable as logical. Did you remember to capitalize? Using TRUE.")
    }
  }else{
    message("Include geneset enrichment value not specified, using TRUE.")
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

# --- Import the data --- 
correlationFileLocation = paste(outputFolderName, filePrefix, "CorrelationFile.rds", sep= "")

if(usePermulationPValOverride){
  permulationFileLocation = permulationPValOverride
}else{
  permulationFileLocation = paste(outputFolderName, filePrefix, "CombinedPrunedFastAllPermulationsPValue.rds", sep= "")
}

correlData = readRDS(correlationFileLocation)                            #Import the correlation data (non-permulated)
if(usePermulations){
  correlData$permPValue = readRDS(permulationFileLocation)                 #Add a collumn to the data with the permulation p Values
}

if(usePermulations){
  targetcolumn = "permPValue"
}else{
  targetcolumn = "p.adj"
}
correlData$permPValue
correlData[[targetcolumn]]

#Make Q values 
qValueObject = qvalue(p = correlData[[targetcolumn]])
correlData$qValue = qValueObject$qvalues

#Make separate objects for positive and negative Rho 
correlDataPositive = correlData[which(correlData$Rho > 0),]
correlDataNegative = correlData[which(correlData$Rho < 0),]
correlDataZero = correlData[which(correlData$Rho == 0),]


# --- run Gene ontology Analysis ---
if(performGeneOntolgy){
  
  
  
  
}

(correlData[[targetcolumn]])


  # --- Display the data ---
if(performGeneOntolgy){
  outputPDFLocation = paste(outputFolderName, filePrefix, "Graphs-PvalQvalGo.pdf", sep= "")
}else{
  outputPDFLocation = paste(outputFolderName, filePrefix, "Graphs-PvalQvalOnly.pdf", sep= "")
}

#function to make text plots
makeTextPlot = function(data, collumn){
  ValueHead = head(data[order(collumn),], n=40)
  ValueHead$N = as.character(ValueHead$N)
  ValueHead$Rho = round(ValueHead$Rho, digits = 5)
  ValueHead$Rho = as.character(ValueHead$Rho)
  ValueHead = format_table(ValueHead, pretty_names = F, digits = "scientific5")
  ValueHead
}


#Make a pdf file
if(performGeneOntolgy){
  pdf(file = outputPDFLocation, width = 15, height = 30)
  par(mfrow = c(4,2))
}else{
  pdf(file = outputPDFLocation, width = 15, height = 15)
  par(mfrow = c(2,2))
}


#plot the top genes
pValueHead = makeTextPlot(correlData, correlData[[targetcolumn]])
textplot(pValueHead, mar = c(0,0,2,0), cmar = 1.5)
title(main = paste("Top genes by",  targetcolumn))
qValueHead = makeTextPlot(correlData, correlData$qValue)
textplot(qValueHead, mar = c(0,0,2,0), cmar = 1.5)
title(main = paste("Top genes by", targetcolumn,  "q-Value"))

#plot p value histograms 
hist(correlDataPositive[[targetcolumn]], breaks = 40, xaxp = c(0,1,20), main = paste("Histogram of correlDataPositive",  targetcolumn))
hist(correlDataNegative[[targetcolumn]], breaks = 40, xaxp = c(0,1,20), main = paste("Histogram of correlDataNegative",  targetcolumn))


if(performGeneOntolgy){
  #Plot the GO Outputs now 
  nonpermEnrichmentFileName = paste(outputFolderName, filePrefix, "EnrichmentFile.rds", sep= "")
  if(file.exists(nonpermEnrichmentFileName)){
    enrichmentResult = readRDS(nonpermEnrichmentFileName)
    enrichmentResult2 = enrichmentResult[[1]]
    makeGOTable = function(data, collumn){
      ValueHead = head(data[order(abs(collumn)),], n=40)
      ValueHead$num.genes = as.character(ValueHead$num.genes)
      ValueHead$stat = round(ValueHead$stat, digits = 5)
      ValueHead$stat = as.character(ValueHead$stat)
      ValueHead = format_table(ValueHead, pretty_names = F, digits = "scientific5")
      #ValueHead$N = str_sub(ValueHead$N, end = -9)
      ValueHead
    }
    enrichHead = makeGOTable(enrichmentResult2, enrichmentResult2$stat)
    textplot(enrichHead, mar = c(0,0,0,0), cmar = 1.5)
    title(main = paste("Top pathways by non-permulation"))
  }else{
    message("No GO file found. Did you remember to run RunEnrichmentAnalysis.R?")
  }
}
dev.off()

