# -- Libraries 
.libPaths("/share/ceph/wym219group/shared/libraries/R4") #add path to custom libraries to searched locations
library(RERconverge)
library(tools)
source("Src/Reu/cmdArgImport.R")


# -- Usage:
# This text describes the purpose of the script 

# -- Command arguments list
# r = filePrefix                      This is a prefix used to organize and separate files by analysis run. Always required. 
# v = <T or F>                        This prefix is used to force the regeneration of the script's output, even if the files already exist. Not required, not always used.
# m = mainTreeFilename.txt or .rds    This sets the location of the maintrees file
# t = <ER or SYM or ARD>              This sets the model type used to estimate ancestral branches; MUST BE THE SAME AS ONE USED IN MakeCategoricalPhenotypeTree.R
# p = <auto or stationary or flat>    This sets the probability of each state at the root of the tree.   
# n = numberOfPermulations            This is the number of permulations to run in the script 
# i = runInstanceValue                This is used to generate unique filenames for each instance of the script. Used in parrallelization. 
# l = relaxationValue <int, 0-1>      This is a value between 0 and 1, the percentage off of exact a permulation is allowed to be. Increases runspeed, but reduces permulation accuracy. Recommend 0.1 if more than 3 categories


#----------------
args = c('r=CategoricalDiet3Phen', 'm=data/RemadeTreesAllZoonomiaSpecies.rds', 'v=F', 't=ER', 'n=20', 'i=Dev', 'l=0.1') #This is a debug argument set. It is used to set arguments locally, when not running the code through a bash script.

args = c('r=NewHiller4Phen', 'm=data/newHillerMainTrees.rds', 'v=F', 't=ER', 'n=4', 'i=Dev', 'l=0.1') #This is a debug argument set. It is used to set arguments locally, when not running the code through a bash script.
args = c('r=NewHiller4Phen', 'm=data/newHillerMainTrees.rds', 'v=F', 't=ER', 'n=4', 'i=Dev', 'l=0.1') #This is a debug argument set. It is used to set arguments locally, when not running the code through a bash script.
args = c('r=CIvHBinaryZoonomia', 'm=data/RemadeTreesAllZoonomiaSpecies.rds', 'v=F', 't=ER', 'n=20', 'i=Dev')

# --- Standard start-up code ---
args = commandArgs(trailingOnly = TRUE)
{  # Bracket used for collapsing purposes
  #File Prefix
  if(!is.na(cmdArgImport('r'))){                                                #This cmdArgImport script is a way to import arguments from the command line. 
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
modelType = "ER"
rootProbability = "auto"
permulationAmount = 100
runInstanceValue = NULL
useRelaxation = FALSE
relaxationValue = NULL

{ # Bracket used for collapsing purposes
  #MainTrees Location
  if(!is.na(cmdArgImport('m'))){
    mainTreesLocation = cmdArgImport('m')
  }else{
    message("No maintrees arg, using default")
  }
  #Model Type
  if(!is.na(cmdArgImport('t'))){
    modelType = cmdArgImport('t')
  }else{
    message("No model specified, using ER.")
  }
  #Root Probability
  if(!is.na(cmdArgImport('p'))){
    rootProbability = cmdArgImport('p')
  }else{
    message("No probability specified, using auto.")
  }
  
  #Number of permulations
  if(!is.na(cmdArgImport('n'))){
    permulationAmount = cmdArgImport('n')
    permulationAmount = as.numeric(permulationAmount)
  }else{
    paste("Number of permulations not specified, using 100")
  }
  #Run instance value
  if(!is.na(cmdArgImport('i'))){
    runInstanceValue = cmdArgImport('i')
  }else{
    paste("This script does not have a run instance value")
  }
  #Relaxation
  if(!is.na(cmdArgImport('l'))){
    useRelaxation = TRUE
    relaxationValue = cmdArgImport('l')
    relaxationValue = as.numeric(relaxationValue)
    source("Src/Reu/RelaxedRejectionPermFuncs2.R")
  }else{
    paste("Not relaxing permulations.")
  }
}


#                   ------- Code Body -------- 

# -- Read in files -- 
#MainTrees
if(file_ext(mainTreesLocation) == "rds"){                                       #if the tree is an RDS file
  mainTrees = readRDS(mainTreesLocation)                                        #Read as RDS
}else{                                                                          #Otherwise
  mainTrees = readTrees(mainTreesLocation)                                      #read as text
}

#Phenotype Vector 
phenotypeVectorFilename = paste(outputFolderName, filePrefix, "CategoricalPhenotypeVector.rds",sep="") #make a filename based on the prefix
phenotypeVector = readRDS(phenotypeVectorFilename)                              #Load the phenotype vector 

#RERs 
#RERFileName = paste(outputFolderName, filePrefix, "RERFile.rds", sep= "")       #Set a filename for the RERs based on the prefix
#RERObject = readRDS(RERFileName)                                                #Load the RERs 

#Paths
#pathsFilename = paste(outputFolderName, filePrefix, "CategoricalPathsFile.rds", sep= "") #make a filename based on the prefix
#pathsObject = readRDS(pathsFilename)                                            #Load the paths

#Correlations
#correlationFileName = paste(outputFolderName, filePrefix, "CorrelationFile.rds", sep= "") #Make a correlation filename based on the prefix
#correlationsObject = readRDS(correlationFileName)                               #Load the correlations

# -- Run Permulations --

permsStartTime = Sys.time()                                                     #get the time before start of permulations
if(useRelaxation){
  permulationData = categoricalPermulations(mainTrees, phenotypeVector, rm = modelType, rp = rootProbability, ntrees = permulationAmount, percent_relax = relaxationValue)
}else{
  permulationData = categoricalPermulations(mainTrees, phenotypeVector, rm = modelType, rp = rootProbability, ntrees = permulationAmount)
}
permsEndTime = Sys.time()                                                       #get time at end of permulations
message(paste("Time to run permulations: ", (permsEndTime - permsStartTime), attr(permsEndTime - permsStartTime, "units"))) #Print the time taken to calculate permulations
message(paste("Time per permulation: ", ((permsEndTime - permsStartTime)/permulationAmount), attr((permsEndTime - permsStartTime)/permulationAmount, "units")))

# -- Save Permulations -- 

saveStartTime = Sys.time()
permulationsDataFileName = paste(outputFolderName, filePrefix, "PermulationsData", runInstanceValue, ".rds", sep= "")
saveRDS(permulationData, permulationsDataFileName)
saveEndTime= Sys.time()
message(paste("Time to save permulations: ", (saveEndTime - saveStartTime), attr(saveEndTime - saveStartTime, "units")))
