# -- Libraries 
.libPaths("/share/ceph/wym219group/shared/libraries/R4") #add path to custom libraries to searched locations
library(RERconverge)
source("Src/Reu/cmdArgImport.R")

# -- Usage:
# This text describes the purpose of the script 

# -- Command arguments list
# r = filePrefix    This is a prefix used to organize and separate files by analysis run. Always required. 
# v = <T or F>      This prefix is used to force the regeneration of the script's output, even if the files already exist. Not required, not always used.


#----------------
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
    paste("Force update not specified, not forcing update")
  }
}

# --- Argument Imports ---
# Defaults
argument1Holder = TRUE

{ # Bracket used for collapsing purposes
  #First argument
  if(!is.na(cmdArgImport('X'))){
    argument1Holder = cmdArgImport('X')
    argumentHolder = as.logical(argument1Holder)
    if(is.na(argument1Holder)){
      argument1Holder = TRUE
      message("Argument value not interpretable as logical. Did you remember to capitalize? Using TRUE.")
    }
  }else{
    message("Argument value not specified, using TRUE.")
  }
}


#                   ------- Code Body -------- 