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
# n = numberOfPermulations            This is the number of permulations to run in the script 
# i = runInstanceValue                This is used to generate unique filenames for each instance of the script. Used in parrallelization. 
# o = outgroupSpecies                 This is the outgroup species with the tip name as it appears on mainTrees$masterTree
# p = <T or F>                        This sets if the script should calculate permulation p vals. ONLY USE THIS IF NOT PARRALLELIZING, as it only considers permulations from this script.
#----------------
args = c('r=MaturityLifespanPercent', 'm=data/newHillerMainTrees.rds', 'v=F', 'n=20', 'i=Dev', 'o=ornAna2', 'p=T') #This is a debug argument set. It is used to set arguments locally, when not running the code through a bash script.


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
permulationAmount = 100
runInstanceValue = NULL
outgroupVal = "Debug"
calculatePVal = F

{ # Bracket used for collapsing purposes
  #MainTrees Location
  if(!is.na(cmdArgImport('m'))){
    mainTreesLocation = cmdArgImport('m')
  }else{
    message("No maintrees arg, using default")
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
  #Outgroup Species
  if(!is.na(cmdArgImport('o'))){
    outgroupVal = cmdArgImport('o')
  }else{
    stop("An outgroup is required, specify outgroup species")
  }
  #Calculate P 
  if(!is.na(cmdArgImport('p'))){
    calculatePVal = cmdArgImport('p')
  }else{
    paste("p caluclation not specified, not interally calculating p values.")
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
phenotypeVectorFilename = paste(outputFolderName, filePrefix, "ContinuousPhenotypeVector.rds",sep="") #make a filename based on the prefix
phenotypeVector = readRDS(phenotypeVectorFilename)                              #Load the phenotype vector 

#RERs 
RERFileName = paste(outputFolderName, filePrefix, "RERFile.rds", sep= "")       #Set a filename for the RERs based on the prefix
RERObject = readRDS(RERFileName)                                                #Load the RERs 

#Paths
#pathsFilename = paste(outputFolderName, filePrefix, "CategoricalPathsFile.rds", sep= "") #make a filename based on the prefix
#pathsObject = readRDS(pathsFilename)                                            #Load the paths

#Correlations
correlationFileName = paste(outputFolderName, filePrefix, "CorrelationFile.rds", sep= "") #Make a correlation filename based on the prefix
correlationsObject = readRDS(correlationFileName)                               #Load the correlations

# -- make rooted Mastertree --

rootedMaster = root.phylo(mainTrees$masterTree, outgroup = outgroupVal, resolve.root = T)
is.binary(rootedMaster)

# -- Run Permulations --

permsStartTime = Sys.time()                                                     #get the time before start of permulations
permulationData = getPermsContinuous(permulationAmount,phenotypeVector,RERObject, trees = mainTrees, mastertree = rootedMaster, calculateenrich = F)
permsEndTime = Sys.time()                                                       #get time at end of permulations
message(paste("Time to run permulations: ", (permsEndTime - permsStartTime), attr(permsEndTime - permsStartTime, "units"))) #Print the time taken to calculate permulations
message(paste("Time per permulation: ", ((permsEndTime - permsStartTime)/permulationAmount), attr((permsEndTime - permsStartTime)/permulationAmount, "units")))

# -- Save Permulations -- 

saveStartTime = Sys.time()
permulationsDataFileName = paste(outputFolderName, filePrefix, "PermulationsData", runInstanceValue, ".rds", sep= "")
saveRDS(permulationData, permulationsDataFileName)
saveEndTime= Sys.time()
message(paste("Time to save permulations: ", (saveEndTime - saveStartTime), attr(saveEndTime - saveStartTime, "units")))

# -- Calculate p values (optional) --

if(calculatePVal){
  permulatedPVals = permpvalcor(correlationsObject, permulationData)
  
  
  permulationsOutputFilename = paste(outputFolderName, filePrefix, "PermulationPValue", sep= "")
  write.csv(permulatedPVals, file= paste(permulationsOutputFilename, ".csv", sep=""), row.names = T, quote = F) #Save correlations as csv
  saveRDS(permulatedPVals, paste(permulationsOutputFilename, ".rds", sep=""))                          #and as an rds 
  
}
