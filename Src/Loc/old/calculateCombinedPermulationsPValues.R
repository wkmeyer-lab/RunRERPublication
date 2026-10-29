.libPaths("/share/ceph/wym219group/shared/libraries/R4") #add path to custom libraries to searched locations
Sys.setenv('R_MAX_VSIZE'=20000000000)
library(RERconverge) #load RERconverge package
library(RERconverge)
library("tools")
source("Src/Reu/cmdArgImport.R")
source("Src/Reu/convertLogiToNumeric.R")
source("Src/Reu/permPValCorReport.R")

#---- USAGE -----
#used to calculate p values on combined permulations files made by CombinePermualtions.R. 

# ARUGMENTS: 
#If an argument contains a '(' it is evaluated as code.
# 'r="filePrefix"'              This is the prefix attached to all files; a required argument. 
# 'i=<number>'                  This is used to generate unique filenames for each instance of the script. Used to target different combination files. 
# 'c=F' OR 'c=T'                This is used to set if the script is being run on previously metacombined permutations. Used for parrallelization. 
# 't = <s OR f OR p>'           This sets which permulation filetype to look for. s is for slow, f is for fast, and p is for pruned-fast

# 'p = <T or F>'                This determines if the p value calculation should be run. Can be turned off to only run GO analysis
# 's = <integer>'               This sets the gene number to start at for parallelization,                 
# 'n = <integer>'               This is the number of genes to do, used to parallelization
# 'd = <T or F>'                This sets if the script should add one to the denominator or not 
# 'g = <T or F>'                This sets if the script should use same-sign-only denominators or not
# 'l = <T or F>'                This sets if the script should use clades-based correlation values, or non-clade-based correlation values 
# 'm = <Filename>'              This is an override for the correlation file target 


# 'e = <T or F>'                This determines if enrichment analysis should be run.
# 'y = <Filename>'              This is an override for the gmt file location  
#-------
#Debug setup defaults
#permulationNumberValue = 3
#Testing args:
args = c('r=demoinsectivory', 'n=3', 'e=F', 't=p')
args = c('r=allInsectivory', 'e=F', 's=1', 't=s', 'i=1')
args = c('r=carnvHerbs', 'e=F', 's=1', 't=s', 'i=1')
args = c('r=CVHRemake', 'e=T', 'p=F')
#------


#---- Prefix Setup -----
#------------------------
{
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

#----Startup----
#----------------
{
  #-Default values -
  
  runInstanceValue = NULL
  metacombineValue = FALSE
  startValue = 1
  geneNumberValue = NA #This means that by default it does all of the genes
  plusOneValue = TRUE
  sameSignValue = FALSE
  useCladesValue = FALSE
  useCorrelationOverride = FALSE
  fileType ="p"
  runPvalueCalculation = TRUE
  runErichmentAnalysis = TRUE
  gmtFileLocation = "Data/enrichmentGmtFile.gmt"
  enrichmentAnnotationListName = "MSigDBPathways"
  #-------
  
  #- Command args import -
  
  {
    #General arguments 
    {
      #instance number of the script
      if(!is.na(cmdArgImport('i'))){
        runInstanceValue = cmdArgImport('i')
      }else{
        message("This script does not have a run instance value")
      }
      
      #Import if this is being run on combined correlations 
      if(!is.na(cmdArgImport('c'))){
        metacombineValue = cmdArgImport('c')
        metacombineValue = as.logical(metacombineValue)
        if(is.na(metacombineValue)){
          metacombineValue = FALSE
          message("Metacombination value not interpretable as logical. Did you remember to capitalize? Using FALSE.")
        }
      }else{
        message("Metacombination value not specified, using FLASE. If you aren't parrallelizing, don't worry about this.")
      }
      
      #Import which filetype to use 
      if(!is.na(cmdArgImport('t'))){
        fileType = cmdArgImport('t')
        fileType = as.character(fileType)
      }else{
        message("No fileType specified, defaulting to pruned-fast (fastPermulationsData)")
      }
      
      #Import the permulation number to start at
      if(!is.na(cmdArgImport('s'))){
        startValue = cmdArgImport('s')
        startValue = as.numeric(startValue)
      }else{
        message("Start value not specified, using 1")
      }
      
      #Import the number of genes to do 
      if(!is.na(cmdArgImport('n'))){
        geneNumberValue = cmdArgImport('n')
        geneNumberValue = as.numeric(geneNumberValue)
      }else{
        message("Number ofgenes not specified, using all")
      }
      
      #Plus one denominator
      if(!is.na(cmdArgImport('d'))){
        plusOneValue = cmdArgImport('d')
        plusOneValue = as.logical(plusOneValue)
        if(is.na(plusOneValue)){
          plusOneValue = FALSE
          message("Plus one denominator value not interpretable as logical. Did you remember to capitalize? Using FALSE.")
        }
      }else{
        message("Plus one denominator value not specified, using TRUE.")
      }
      
      #same sign denominator
      if(!is.na(cmdArgImport('g'))){
        sameSignValue = cmdArgImport('g')
        sameSignValue = as.logical(sameSignValue)
        if(is.na(sameSignValue)){
          sameSignValue = FALSE
          message("Same sign denominator value not interpretable as logical. Did you remember to capitalize? Using FALSE.")
        }
      }else{
        message("Same sign denominator value not specified, using FLASE.")
      }
      
      #Use Clades
      if(!is.na(cmdArgImport('l'))){
        useCladesValue = cmdArgImport('l')
        useCladesValue = as.logical(useCladesValue)
        if(is.na(useCladesValue)){
          useCladesValue = TRUE
          message("Use Clades value not interpretable as logical. Did you remember to capitalize? Using FLASE.")
        }
      }else{
        message("Use Clades value not specified, using FALSE.")
      }
      
      #Correlation Override
      if(!is.na(cmdArgImport('m'))){
        useCorrelationOverride = TRUE
        correlationOverride = cmdArgImport('m')
      }else{
        message("No correlation override")
      }
      
      #Calculate p values
      if(!is.na(cmdArgImport('p'))){
        runPvalueCalculation = cmdArgImport('p')
        runPvalueCalculation = as.logical(runPvalueCalculation)
        if(is.na(runPvalueCalculation)){
          runPvalueCalculation = TRUE
          message("Run p-value Calculation not interpretable as logical. Did you remember to capitalize? Using TRUE.")
        }
      }else{
        message("Run p-value Calculation not specified, using TRUE.")
      }
    }
    
    #Enrichment arguments
    {
      #Run Enrichment analysis
      if(!is.na(cmdArgImport('e'))){
        runErichmentAnalysis = cmdArgImport('e')
        runErichmentAnalysis = as.logical(runErichmentAnalysis)
        if(is.na(runErichmentAnalysis)){
          runErichmentAnalysis = TRUE
          message("Run Enrichment analysis not interpretable as logical. Did you remember to capitalize? Using TRUE.")
        }
      }else{
        message("Run Enrichment analysis not specified, using TRUE.")
      }
      
      #gmt file location
      if(!is.na(cmdArgImport('y'))){
        gmtFileLocation = cmdArgImport('y')
      }else{
        paste("No gmt location argument, using Data/enrichmentGmtFile.gmt")                          #Report using default
        message("No gmt location argument, using Data/enrichmentGmtFile.gmt")
      }
    }
  }
}


#---- MAIN CODE ----
#-------------------


# --- Get input files ---

# -- Get correlation file--
if(useCorrelationOverride){                                                     #If using an override
  overrideCorellationFileName = paste(outputFolderName, filePrefix, correlationOverride, sep= "")
  if(file.exists(paste(overrideCorellationFileName, ".rds", sep=""))){
    overrideCladesCorrelation = readRDS(paste(overrideCorellationFileName, ".rds", sep=""))
  }else{
    message("Override Correlation file does not exist, p-values not calculated. (RunRERBinaryMT.R)")
    stop("Requires Override Correlation file.")
  }
  correlationSet = overrideCladesCorrelation
}else if(useCladesValue){                                                       #If using clades
  cladesCorellationFileName = paste(outputFolderName, filePrefix, "CladesCorrelationFile", sep= "")
  if(file.exists(paste(cladesCorellationFileName, ".rds", sep=""))){
    cladesCorrelation = readRDS(paste(cladesCorellationFileName, ".rds", sep=""))
  }else{
    message("Clades Correlation file does not exist, p-values not calculated. (sisterListGeneration.R)")
    stop("Requires Clades Correlation file.")
  }
  correlationSet = cladesCorrelation
}else{                                                                          #Otherwise
  noncladesCorellationFileName = paste(outputFolderName, filePrefix, "CorrelationFile", sep= "")
  if(file.exists(paste(noncladesCorellationFileName, ".rds", sep=""))){
    noncladesCorrelation = readRDS(paste(noncladesCorellationFileName, ".rds", sep=""))
  }else{
    message("Non-Clades Correlation file does not exist, p-values not calculated. (RunRERBinaryMT.R)")
    stop("Requires Non-Clades Correlation file.")
  }
  correlationSet = noncladesCorrelation
}
# ---- Get Combined Permulations filename ---
#Do the first combination:
if(fileType == 's'){
  fileTypeString = ''
}else if(fileType == 'f'){
  fileTypeString = "UnprunedFast"
}else if(fileType == 'p'){
  fileTypeString = "PrunedFast"
  message("The default permulation filetype is being used. IF CANNOT READ IN FILES, ensure the correct filetype is specified.")
}else{
  stop( "THIS IS AN ISSUE MESSAGE, IMPROPER FILETYPE ARGUMENT")
}

if(metacombineValue == F){
  combinedDataFileName = paste(outputFolderName, filePrefix, "Combined", fileTypeString,"PermulationsData", runInstanceValue, ".rds", sep="")
}else{
  combinedDataFileName = paste(outputFolderName, filePrefix, "MetaCombined", fileTypeString, "PermulationsData", runInstanceValue, ".rds", sep="")
}
message(combinedDataFileName)
timeBeforeReadIn = Sys.time()
combinedPermulationsData = readRDS(combinedDataFileName)
timeAferReadIn = Sys.time()
readingTime = timeAferReadIn - timeBeforeReadIn
message("File reading time: ", readingTime, attr(readingTime, "units"))
  
#Explict order for garbage collection 
gc() 

# -- run pValue calculation --
if(runPvalueCalculation){
  permulationPValues = permPValCorReport(correlationSet, combinedPermulationsData, startNumber = startValue, geneNumber = geneNumberValue, plusOne = plusOneValue, signedDenominator = sameSignValue)
  
  #save the permulations p values
    #Make generate a marker saying which genes were calculated 
  if(is.na(geneNumberValue)){
    geneRange = "All"
  }else{
    geneRange = paste(startValue, "-", (startValue + geneNumberValue -1), sep = "" )
  }
  
  permulationPValueFileName = paste(outputFolderName, filePrefix, "Combined", fileTypeString, geneRange, "PermulationsPValue", runInstanceValue, ".rds", sep= "")
  saveRDS(permulationPValues, file = permulationPValueFileName)
}

# -- run GO Analysis --
if(runErichmentAnalysis){
  # Set up the annotations list
  gmtAnnotations = read.gmt(gmtFileLocation)
  annotationsList = list(gmtAnnotations)
  enrichmentListName = substring(gmtFileLocation, 6, last = (nchar(gmtFileLocation) - 4))
  names(annotationsList) = enrichmentListName
  
  #Get or make the enrichment file
  nonpermEnrichmentFileName = paste(outputFolderName, filePrefix, "EnrichmentFile.rds", sep= "")
  if(!file.exists(nonpermEnrichmentFileName)){
    message("No Enchriment File found. Generating enrichment file.")
    rerStats = getStat(correlationSet)
    enrichmentResult = fastwilcoxGMTall(rerStats, annotationsList, outputGeneVals = F)
    #save the enrichment 
    enrichmentFileName = paste(outputFolderName, filePrefix, "EnrichmentFile.rds", sep= "")
    saveRDS(enrichmentResult, enrichmentFileName)
  }else{
    enrichmentResultReadIn = readRDS(nonpermEnrichmentFileName)
    generateEnrichment = F
  }
  
  
  if(file.exists(nonpermEnrichmentFileName)){
    enrichmentResultReadIn = readRDS(nonpermEnrichmentFileName)
      if(names(enrichmentResultReadIn) == enrichmentListName) {# check to make sure that the enrichment file is using the same gmt file
        enrichmentResult = enrichmentResultReadIn
      }else{
        message("Existing enrichment file using different gmt file. Replacing enrichment file.")
        generateEnrichment = T
      }
  }else{
    message("No Enchriment File found. Generating enrichment file.")
    generateEnrichment = T
  }
    
  if(generateEnrichment){
    rerStats = getStat(correlationSet)
    enrichmentResult = fastwilcoxGMTall(rerStats, annotationsList, outputGeneVals = F, , num.g =4)
    #save the enrichment 
    enrichmentFileName = paste(outputFolderName, filePrefix, "EnrichmentFile.rds", sep= "")
    saveRDS(enrichmentResult, enrichmentFileName)
  }
  
  timeAtAddEnrichmentStart = Sys.time()
  permulationsWithEnrichments = getEnrichPerms(combinedPermulationsData, enrichmentResult, annotationsList)
  timeAtAddEnrichmentEnd = Sys.time()
  timeForEnrichments = timeAtAddEnrichmentEnd - timeAtAddEnrichmentStart
  message("Time to add enrichments to permulations: ", timeForEnrichments, attr(timeForEnrichments, "units"))
  
  enrichmentPermValues = permpvalenrich(enrichmentResult, permulationsWithEnrichments)
  timeEnrichmentCalculationEnd = Sys.time()
  timeForEnrichmentCalculation = timeEnrichmentCalculationEnd - timeAtAddEnrichmentEnd 
  message("Time to calculate enrichments: ", timeForEnrichmentCalculation, attr(timeForEnrichmentCalculation, "units"))
  
  #save the output 
  enrichmentPermulationFileName = paste(outputFolderName, filePrefix, "EnrichmentPermulationsValue.rds", sep= "")
  saveRDS(enrichmentPermValues, file = enrichmentPermulationFileName)
  
}
