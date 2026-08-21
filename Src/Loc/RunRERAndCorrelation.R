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
# r = filePrefix                                                               REQUIRED: This is a prefix used to organize and separate files by analysis run. Always required. 
# m = mainTreeFilename.txt or .rds                                             REQUIRED: This sets the location of the maintrees file
# v = <T or F>                                                                 This prefix is used to force the regeneration of the script's output, even if the files already exist. Not required, not always used.
# s = < ["b" or "binary"] or ["c" or "continuous"] or ["g" or "categorical"]>  This prefix is used to set the type of phenotype being supplied; it will detect automatically if not specified.
# c = < "diff" or "mean" or "last" >                                           This is used for continuous traits, to determine if the metic should be the difference between the nodes (diff), the mean(mean of the two nodes), or last(the downstream value). Note that Mean and Last are not phylogenetically independent, and do not have downstream processing. 
# l = <min.sp value>                                                           This sets the min.sp value to be used in the correlation. 
# f = speciesFilterText                                                        OPTIONAL OVERRIDE: This can be used to manually specify a species filer; leave blank for automatic
# p = phenotypeTreeFilename.txt or .rds                                        OPTIONAL OVERRIDE: This can be used to manually override the phenotype tree being used. For continuous analyses, this is the location of the trait vector.


#----------------

args = c("r=ComplexDietCentralAnalysis", 'm=zoonomiaAllMammalsTrees.rds', "s=g", "v=F", "l=170")


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
    message("No phenotype style specified, attempting to detect automatically.")
    
    #Check for categorical 
    phenotypeTreeFilename = paste(outputFolderName, filePrefix, "Categorical", "Tree.rds", sep="")
    if(file.exists(phenotypeTreeFilename)){
      message("Found information for Categorical phenotype, using Categorical phenotype")
      phenotypeStyle = "Categorical"
    }else{
      phenotypeTreeFilename = paste(outputFolderName, filePrefix, "Continuous", "PhenotypeVector.rds", sep="")
      if(file.exists(phenotypeTreeFilename)){
        message("Found information for continuous phenotype, using continuous phenotype")
        phenotypeStyle = "Continuous"
      }else{
        phenotypeTreeFilename = paste(outputFolderName, filePrefix, "Binary", "Tree.rds", sep="")
        if(file.exists(phenotypeTreeFilename)){
          message("Found information for binary phenotype, using binary phenotype")
          phenotypeStyle = "binary"
        }else{
          message("No phenotyep specified and no correct tree found.")
          stop()
        }
      } 
    }
    
    

    
    
    message("PhenotypeStyle not specified, using continuous")
    
  }
  
  #phenotype tree location
  if(!is.na(cmdArgImport('p'))){
    phenotypeTreeLocation = cmdArgImport('p')
  }else{                                                                        #See if a pre-made tree for this prefix and style exists 
    if(phenotypeStyle == "Continuous"){
      phenotypeTreeFilename = paste(outputFolderName, filePrefix, phenotypeStyle, "PhenotypeVector.rds", sep="")
    }else{
      phenotypeTreeFilename = paste(outputFolderName, filePrefix, phenotypeStyle, "Tree.rds", sep="")
    }
    
    if(file.exists(paste(phenotypeTreeFilename))){                              #if so, use it                
      phenotypeTreeLocation = phenotypeTreeFilename                     
      message("Pre-made Phenotype tree found, using pre-made tree.")
    }else{
      #paste("THIS IS AN ISSUE MESSAGE; SPECIFY PHENOTYPE TREE")
      stop("THIS IS AN ISSUE MESSAGE; NO PHENOTYPE TREE/VECTOR FOUND. ENSURE PHENOTYPE TREE/VECTOR FOR SELECTED STYLE AVAILABLE AND PREFIX IS CORRECT.")
    }
  }
  
  #Species Filter   
  if(!is.na(cmdArgImport('f'))){
    speciesFilter = cmdArgImport('f')
  }else{ #See if a pre-made filter for this prefix exists 
    speciesFilterFileName = paste(outputFolderName, filePrefix, "SpeciesFilter.rds",sep="") #Make the name of the location a pre-made filter would have to test for it
    if (file.exists(paste(speciesFilterFileName))){                  
      speciesFilter = readRDS(speciesFilterFileName)                       #if so, use it 
      message("Pre-made filter found, using pre-made filter.")
    }else{                                                    
      message("No speciesFilter arg, using all species")                           #if not, use no filter
    }
  }
  
  #Continuous Metric 
  if(!is.na(cmdArgImport('c'))){
    continousMetric = cmdArgImport('c')
    if(!any(continousMetric == validMetrics)){
      stop("Contious metric is not a valid option. Use 'diff', 'mean', or 'last'")
    }
  }else{                                                                        #If a continuous phenotype, report using Diff
    if(phenotypeStyle == "Continuous"){
      message("No continuous metric specified, using diff")
    }
  }
  #min.sp Value
  if(!is.na(cmdArgImport('l'))){
    minSpValue = cmdArgImport('l')
  }else{
    message("min.sp not specified, using 10")
  }
}

#                   ------- Code Body -------- 

# -- Read in the trees --
#MainTrees
if(file_ext(mainTreesLocation) == "rds"){
  if(!exists("mainTrees")){mainTrees = readRDS(mainTreesLocation)}
}else{
  if(!exists("mainTrees")){mainTrees = readTrees(mainTreesLocation)} 
}
#Phenotype tree
if(file_ext(phenotypeTreeLocation) == "rds"){                                   #if the tree is an RDS file
  phenotypeTree = readRDS(phenotypeTreeLocation)                                #Read as RDS
}else{                                                                          #Otherwise
  phenotypeTree = readTrees(phenotypeTreeLocation)                              #read as text
}

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

if(phenotypeStyle == "Binary"){                                                 #If binary
  correlation = correlateWithBinaryPhenotype(RERObject, pathsObject, min.sp =minSpValue)#Correlate with binary phenotype
  
  #Generate a phenotype vector 
  fgEdgeObjects = phenotypeTree$edge[which(phenotypeTree$edge.length>=1) ,]                                        #Make an object of the edges in the foreground. This is used as opposed to just referencing the tree directly to allow for "walking" in the final loop of the code
  foregroundNodes = which(1:length(phenotypeTree$tip.label) %in% as.vector(fgEdgeObjects))
  foregroundSpecies = phenotypeTree$tip.label[foregroundNodes]
  phenotypeVector = c(0,0);length(phenotypeVector) = length(phenotypeTree$tip.label);phenotypeVector[] = 0 
  names(phenotypeVector) = phenotypeTree$tip.label
  phenotypeVector[(names(phenotypeVector) %in% foregroundSpecies)] = 1
  phenotypeVector = phenotypeVector[names(phenotypeVector) %in% speciesFilter]
  phenotypeVectorFilename = paste(outputFolderName, filePrefix, "phenotypeVector.rds", sep="")
  saveRDS(phenotypeVector, file = phenotypeVectorFilename)
  
}else if(phenotypeStyle == "Continuous"){                                        #if continuous
  correlation = correlateWithContinuousPhenotype(RERObject, pathsObject, min.sp = minSpValue)
  
} else if(phenotypeStyle == "Categorical"){                                     #if categorical
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

