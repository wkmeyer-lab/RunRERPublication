#Library setup 
.libPaths("/share/ceph/wym219group/shared/libraries/R4") #add path to custom libraries to searched locations
library(RERconverge) #load RERconverge package
library(RERconverge)
library("tools")
#library(paleotree)
source("Src/Reu/cmdArgImport.R")
source("Src/Reu/convertLogiToNumeric.R")


# --- Debug settup---
#batTree = readRDS("Output/allInsectivory/allInsectivoryBinaryForegroundTree.rds")


#debugPlotTree = function(debugTree){dev.off(); dev.new(); dev.new(); testplot2 = plotTreeHighlightBranches(debugTree,hlspecies=which(debugTree$edge.length== 1),hlcols="blue", main="Input tree"); edgelabels(cex = 0.7, frame="none", font=2, adj=c(0,-0.2), col="blue"); nodelabels(cex = 0.7, frame="none", font=2, adj=c(-0.2,0.3), col="dark green"); tiplabels(cex = 0.8, frame="none", font=2, adj=c(0.2,0), col="dark red")}
#debugPlotTree(inputTree)
#debugPlotTree(unprunedInputTree)
#inputTree$edge.length[which(inputTree$edge.length >1)] = 1





# ----- USAGE -----
#This is used to generate sistersList and foregroundStrings for permulations in RER converge. 
#This expects:
  # A binary foreground tree, generated using the makeTreeBinary script in this project. 
  # An argument of a file prefix, for the folder to output to and find the binary tree. 
#This produces: 
  #A sistersList list, exported as an RDS. 
  #A foreground string, exported as an RDS. 

#ARGUMENTS: 
# 'r="filePrefix"'                 This is the prefix attached to all files; a required argument.
# 'm=mainTreeFilename.txt or .rds' This is the location of the maintree file. Accepts .txt or .rds. Only used for creation of clades correlation files. 
#  t = <T or F>                    This value determines if the phenotype vector should be trimmed to the species filter
#  v = <T or F>                    This value forces the regeneration of output files which otherwise would not be (CladesPathsFile.rds; CladesCorrelationFile.rds; CladesCorrelationFile.csv). 


# ------ Command Line Imports:

#testing args: 
args = "r=demoInsectivory"
args = "r=allInsectivory"
args = c("r=carnvHerbsOldRemake", "m=Data/RemadeTreesAllZoonomiaSpecies.rds", "v=F")

# --- Import prefix --- 

args = commandArgs(trailingOnly = TRUE)

#File Prefix
if(!is.na(cmdArgImport('r'))){
  filePrefix = cmdArgImport('r')
}else{
  stop("THIS IS AN ISSUE MESSAGE; SPECIFY FILE PREFIX")
}


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

# ----- Import force update argument
forceUpdate = FALSE

#Import if update being forced with argument 
if(!is.na(cmdArgImport('v'))){
  forceUpdate = cmdArgImport('v')
  forceUpdate = as.logical(forceUpdate)
}else{
  paste("Force update not specified, not forcing update")
}

# ---------------

# ----- Command Args Import -----
speciesFilter = NULL
trimPhenotypeVector = TRUE

{
#speciesFilter
speciesFilterFileName = paste(outputFolderName, filePrefix, "SpeciesFilter.rds",sep="") #Make the name of the location a pre-made filter would have to test for it

if(!is.na(cmdArgImport('f'))){
  speciesFilter = cmdArgImport('f')
}else if (file.exists(paste(speciesFilterFileName))){                  #See if a pre-made filter for this prefix exists 
  speciesFilter = readRDS(speciesFilterFileName)                       #if so, use it 
  paste("Pre-made filter found, using pre-made filter.")
}else{                                                    
  paste("No speciesFilter arg, using NULL")                           #if not, use no filter
}

#Maintrees, with check to see if needed 
{
cladesPathsFileName = paste(outputFolderName, filePrefix, "CladesPathsFile.rds", sep= "")
cladesCorellationFileName = paste(outputFolderName, filePrefix, "CladesCorrelationFile", sep= "")
if(!file.exists(paste(cladesPathsFileName)) | !file.exists(paste(cladesCorellationFileName, ".rds", sep="")) | forceUpdate){
  if(!file.exists(paste(cladesPathsFileName)) | forceUpdate){
    mainTreesLocation = "/share/ceph/wym219group/shared/projects/MammalDiet/Zoonomia/RemadeTreesAllZoonomiaSpecies.rds"
    if(!is.na(cmdArgImport('m'))){
      mainTreesLocation = cmdArgImport('m')
    }else{
      paste("No maintrees arg, using default")                          #Report using default
      message("No maintrees arg, using default")
    }
    mainTrees = readRDS(mainTreesLocation)
  }
}
}
}

#Import if trim phenotype vector
if(!is.na(cmdArgImport('f'))){
  trimPhenotypeVector = cmdArgImport('f')
  trimPhenotypeVector = as.logical(trimPhenotypeVector)
}else{
  paste("Trimming phenotype vector not specified, using TRUE")
}





# ---- Main code ---- 
inputTreeFilename = paste(outputFolderName, filePrefix, "BinaryForegroundTree.rds", sep="")
inputTree= readRDS(inputTreeFilename) 
if(trimPhenotypeVector){
  unprunedInputTree = inputTree
  inputTree = pruneTree(inputTree, speciesFilter)
}

# -- Foreground Edges --
fgEdges = which(inputTree$edge.length>=1)                                       #Make a list of the edges in the foreground
fgEdgeObjects = inputTree$edge[fgEdges,]                                        #Make an object of the edges in the foreground. This is used as opposed to just referencing the tree directly to allow for "walking" in the final loop of the code



# -- Startup of variables for loops --
cladNumber = 1                                                                  #This is the number of clade being worked on; eg when cladNumber = 1, the output clade would be "clade1"
compeltedStartNodes = NULL                                                      #Create a completed startnodes variable
WrapperEdges = NULL                                                             #Create a variable to store the edges which are wrappers
run =1                                                                          #set the run number to one
for(i in 1:40){cladName = paste("clade",i,sep=""); rm(list=cladName)}; rm(cladName) #Remove any old clade objects. This is useful when debugging the script
cladList = NULL                                                                 #Create a list of clades, stored as an int of the start node, and a name with the clade name associated with that node
currentCladName = NULL                                                          #Create a current clade name variable, and also clear it if it existed
 
addToCladList = function(startNode, inList = cladList){
  cladEntry = startNode                                                         #Store the start node as CladEntry
  currentCladName = paste("clade", cladNumber, sep = "")                        #Update the clade name
  names(cladEntry) = currentCladName                                            #Set the name of the clade entry to the current clade 
  message("Clade entry:", cladEntry)                                            #Message the node 
  message("Clade number: ", cladNumber)                                         #Message the clade associated with that node
  outList = append(inList, cladEntry)                                           #Add this entry to the clade list
  return(outList)                                                               #Return the updated clade list
}

# ---- Loop for initial categorization and clades ----
for(i in fgEdges){
  
  # - Message - 
  message("-----")
  message("First step run # ", run)
  message("Edge name ", i)
  run=run+1
  currentCladName = paste("clade", cladNumber, sep = "")                        #Set the current clade name
  
  # - start node -
  startNode = inputTree$edge[i,1]                                               #Get the node the branch starts at 
  message("Start node: ", startNode)                                            #Report the start node 
  if(startNode %in% compeltedStartNodes){message("repeat"); next}               #If the start node has already been evaluated, skip it
  compeltedStartNodes = append(compeltedStartNodes, startNode)                  #If not, add this start node to the list of evaluated start nodes 

  
  # - end nodes -
  endNodes = fgEdgeObjects[which(fgEdgeObjects[,1] == startNode),2]             #Get the nodes which the start node connects to
  message("End nodes:", endNodes[1], " ", endNodes[2])                          #report the end nodes
  
  #If for some reason the start node has more than two child nodes, warn the user.
  if(!(length(endNodes)==1 | length(endNodes)==2)){message(paste("Something has gone wrong; foreground branch ", i, "has incorrect number of child nodes."))}
  
  # --- Classification of clade type ----
  cladType = NULL
  
  if(length(endNodes) == 1 & all(endNodes<length(inputTree$tip.label))){        #If there is only one child node, and it is a node with a tip label
    cladType = "solo"                                                           #Categorize it as a solo clade
    message(cladType)                                                           #Report clade type
  }
  
  else if(all(endNodes <= length(inputTree$tip.label))){                        #If there are two child nodes, and they both have a tip label
    cladType = "starter"                                                        #categorize it as a starter clade
    message(cladType)                                                           #Report clade type
  }
  
  else if(length(endNodes) == 1 & all(endNodes>length(inputTree$tip.label))){   #If there is only one child node, and it is not labeled (internal node)
    cladType= "internal"                                                        #categorize it as an internal clade
    #These are produced via the terminal/all/ancestral setup, and thus are handled elsewhere
    message(cladType)                                                           #Report clade type
    next                                                                        #Go to next clade 
  }
  
  else if(any(endNodes <= length(inputTree$tip.label)) & any(endNodes > length(inputTree$tip.label))){  #If there are two child nodes, one is labeled and one is not
    cladType = "firstWrapper"                                                                           #categorize it as a first wrapper clade (one clade, one species)
    WrapperEdges = append(WrapperEdges, i)                                      #Add to the list of edges for the next loop to handle
    message(cladType)                                                           #Report clade type
    next                                                                        #Go to next clade
  }
  
  else if(all(endNodes > length(inputTree$tip.label))){                         #If there are two child nodes, and they are both internal                     
    cladType = "metaWrapper"                                                    #categorize as a Metawrapper clade (a clade containing two clades)
    WrapperEdges = append(WrapperEdges, i)                                      #Add to the list of edges for the next loop to handle
    message(cladType)                                                           #Report clade type
    next                                                                        #Go to next clade
  }
  

  
  # -- Handle first-level clades --
  
  if(cladType == "solo"){                                                       #if this is a solo clade 
    speciesOne = inputTree$tip.label[endNodes[1]]                               #Set the species name as the end node's tip label
      #---- old single-species clade Ver ----
      #assign(currentCladName, c(speciesOne))                                      #assign that species name to the current clade name
      #cladList = addToCladList(startNode)
      #cladNumber = cladNumber+1                                                   #increase the clade number
      #------
    soloEntry = startNode                                                       #Store the start node as CladEntry
    names(soloEntry) = speciesOne                                               #Associate this node with he species name in the cladelist 
    message("Solo node:", soloEntry)                                            #Message the node 
    message("Species associated: ", speciesOne)                                 #Message the clade associated with that node
    cladList = append(cladList, soloEntry)                                      #Add this clade to the clade list
  }
  else if(cladType == "starter"){                                               #If this is a start clade
    speciesOne = inputTree$tip.label[endNodes[1]]                               #Set the first species name to the tip label of the first node
    speciesTwo = inputTree$tip.label[endNodes[2]]                               #Set the second species name to the tip label of the second node
    assign(currentCladName, c(speciesOne, speciesTwo))                          #Assign those species names to the current clade name
    cladList = addToCladList(startNode)                                         #Add this clade to the clade list
    cladNumber = cladNumber+1                                                   #increase the clade number
  }
}



# ----- Loop for clades which wrap other clades ----
remainingWrapperEdges = WrapperEdges                                            #Set the edges still to do to be all of the wrapper edges
stuckCycles = 0                                                                 #Set the number of stuck cycles to zero 
run = 1                                                                         #set the run number back to one 
while(length(remainingWrapperEdges >0) & stuckCycles < 20){                     #While there are remaining edges, and the loop hasn't been stuck for twenty cycles:
  for(i in remainingWrapperEdges){
    # - message - 
    message("-----")
    message("Second step Run # ", run)
    message("Edge name: ", i)
    startRemainingEdges = length(remainingWrapperEdges)                         #Set the number of start remaining edges before the loop for a later comparison
    message(paste("Length of remaining Edge Wrappers = ", startRemainingEdges))
    message(paste("number of stuck cycles = ", stuckCycles))
    run=run+1
    currentCladName = paste("clade", cladNumber, sep = "")
    
    # - start and end nodes
    startNode = inputTree$edge[i,1]                                             #get the node the edge starts at 
    message("startnode: ", startNode)
    endNodes = fgEdgeObjects[which(fgEdgeObjects[,1] == startNode),2]           #get the nodes that the start node connects to, from a list made earlier
    message("Endnodes: ", endNodes[1], " ", endNodes[2])
    
    #If the node had more than two child nodes, notify the user
    if(!(length(endNodes)==1 | length(endNodes)==2)){message(paste("Something has gone wrong; foreground branch ", i, "has incorrect number of child nodes."))}
    
    firstNode = endNodes[1]
    firstNodeValid = (firstNode <= length(inputTree$tip.label) | firstNode %in% cladList ) #Set the first node as valid if either it had a tip label, or it has been assigned a clade
    message("first node valid: ", firstNodeValid)                               #Report if the node is valid
      
    secondNode = endNodes[2]
    secondNodeValid = (secondNode <= length(inputTree$tip.label) | secondNode %in% cladList )#Set the second node as valid if either it had a tip label, or it has been assigned a clade
    message("secondNodeValid: ", secondNodeValid)                               #Report if the node is valid
    
    #- Evaluate the node - 
    
    if((firstNodeValid & secondNodeValid)){                                     #If Both end nodes have been processed, evaluate this edge. If not, it will be skipped, and tried again on a subsequent pass, hopefully after the end nodes have been processed.
      remainingWrapperEdges = remainingWrapperEdges[!remainingWrapperEdges %in% i] #Remove this edge from the to-do list
      
      #- get the first node's name - 
      if(firstNode <= length(inputTree$tip.label)){                             #If the node has a tip label
        firstNodeName = inputTree$tip.label[firstNode]                          #The name is the tip label
      }else{                                                                    #if not
        firstNodeName = names(which(cladList == firstNode))                     #The name is the node's assigned clade name
      }
      message("firstNodeName: ", firstNodeName)                                 #Report the first node name
      
      #- get the second node's name -
      if(secondNode <= length(inputTree$tip.label)){                            #If the node has a tip label
        secondNodeName = inputTree$tip.label[secondNode]                        #The name is the tip label
      }else{                                                                    #if not
        secondNodeName = names(which(cladList == secondNode))                   #The name is the node's assigned clade name
      }
      message("secondNodeName: ", secondNodeName)                               #Report the second node name
      
      #Make the clade object
      assign(currentCladName, c(firstNodeName, secondNodeName))                 #assign those names to the current clade name
      
      #Add the clade to the list
      cladList = addToCladList(startNode)                                       #Add the clade to the cladelist
      cladNumber = cladNumber+1                                                 #update the clade number
      
      
    }else{                                                                      #from before, for if one of the child nodes hasn't been processed yet
      message("Node not valid, skipping branch")                                #notify the user
    }
    
    #check if it was stuck during this loop
    endRemainingEdges = length(remainingWrapperEdges)                           #get the length of the remaining edges now that the loop has run
    if(startRemainingEdges == endRemainingEdges){                               #see if it has not changed this loop
      stuckCycles = stuckCycles+1                                               #If it hasn't, increase the number of stuck cycles
    }else{                                                                      #If it has changed
      stuckCycles = 0                                                           #Set stuck cycles to zero 
    }
  }
}

# ---- Loop to check for background nodes inside an otherwise foreground clade; which cause issues in the above loop ----
run = 1                                                                         #Set run to one
stuckCycles = 0                                                                 #Set stuck cycles to zero 
while(length(remainingWrapperEdges >0) & stuckCycles < 20){
  for(i in remainingWrapperEdges){
    # - Message - 
    message("-----")
    message("Resolve Internals Run # ", run)
    message("Edge name: ", i)
    message(paste("Remaining Edges = ", remainingWrapperEdges))
    message(paste("number of stuck cycles = ", stuckCycles))
    run=run+1
    startRemainingEdges = length(remainingWrapperEdges)
    currentCladName = paste("clade", cladNumber, sep = "")
    
    # - start and end nodes -
    startNode = inputTree$edge[i,1]                                             #get the node the edge starts at 
    message("startnode: ", startNode)
    endNodes = fgEdgeObjects[which(fgEdgeObjects[,1] == startNode),2]           #get the nodes that the start node connects to, from a list made earlier
    message("Endnodes: ", endNodes[1], " ", endNodes[2])
    
    #If the node had more than two child nodes, notify the user
    if(!(length(endNodes)==1 | length(endNodes)==2)){message(paste("Something has gone wrong; foreground branch ", i, "has incorrect number of child nodes."))}
    
    
    # ---- check if any of the grandchildren are background, if so, move the end node to the foreground grandchild -----
    #What this does is tell the code to make the associate clade a sister to the clade one step downstream of the node with a foreground child, instead of trying to work on the node with a foreground child and not understanding why there is only one end node
    
    # - check if the first end node's children are foreground
    firstNode = endNodes[1]                                                     #set the first node
    firstNodeChildren = inputTree$edge[which(inputTree$edge[,1] == firstNode),2]#set the children as the end nodes of using the "first node" as the start node. These are basically the grandchildren of the original start node.
    message("")
    message("First Node: ", firstNode)
    message("First Node Children: ", firstNodeChildren[1], " ", firstNodeChildren[2])
    
    if(any(!(firstNodeChildren %in% fgEdgeObjects))){                           #If either of the grandchildren of the first end node are not in the foreground
      if((!(firstNodeChildren[1] %in% fgEdgeObjects))){                         #If it's the first one        
        message(firstNodeChildren[1], " is not foreground")                     #tell the user
        message("Replacing ", firstNode, " with child ", firstNodeChildren[2])   #Tell the user how you are fixing the problem
        firstNode = firstNodeChildren[2]                                        #set the "first node" as the valid grandchild instead, effectively skipping the node which has a non-foreground child
        fgEdgeObjects[which(fgEdgeObjects[,1] == startNode),2][1] = firstNodeChildren[2] #Do the same to the children list which is stable outside this instance of the loop. By doing this, it can progressively walk the node down the chain one step each time it is run through, to handle multiple background child nodes.
      }
      if((!(firstNodeChildren[2] %in% fgEdgeObjects))){                         #If it's the second
        message(firstNodeChildren[2], " is not foreground")                     #tell the user
        message("Replacing ", firstNode, " with child ", firstNodeChildren[1])   #Tell the user how you are fixing the problem
        firstNode = firstNodeChildren[1]                                        #set the "first node" as the valid grandchild instead.
        fgEdgeObjects[which(fgEdgeObjects[,1] == startNode),2][1] = firstNodeChildren[1] #Do the same to the children list which is stable outside this instance of the loop.
      }
    }
    
    # - check if the second end node's children are foreground -
    
    secondNode = endNodes[2]
    secondNodeChildren = inputTree$edge[which(inputTree$edge[,1] == secondNode),2]
    message("")
    message("Second Node: ", secondNode)
    message("Second Node Children: ", secondNodeChildren[1], " ", secondNodeChildren[2])
    
    if(any(!(secondNodeChildren %in% fgEdgeObjects))){                          #If either of the grandchildren of the second end node are not in the foreground
      if((!(secondNodeChildren[1] %in% fgEdgeObjects))){                        #If it's the first one
        message(secondNodeChildren[1], " is not foreground")                    #tell the user
        message("Replacing ", secondNode, " with child ", secondNodeChildren[2]) #Tell the user how you are fixing the problem
        secondNode = secondNodeChildren[2]                                      #set the "second node" as the valid grandchild instead.
        fgEdgeObjects[which(fgEdgeObjects[,1] == startNode),2][2] = secondNodeChildren #Do the same to the children list which is stable outside this instance of the loop.
        }
      if((!(secondNodeChildren[2] %in% fgEdgeObjects))){                        #if It's the second one
        message(secondNodeChildren[2], " is not foreground")                    #tell the user
        message("Replacing ", secondNode, " with child ", secondNodeChildren[1]) #Tell the user how you are fixing the problem
        secondNode = secondNodeChildren[1]                                      #set the "second node" as the valid grandchild instead.
        fgEdgeObjects[which(fgEdgeObjects[,1] == startNode),2][2] = secondNodeChildren[1] #Do the same to the children list which is stable outside this instance of the loop.
        }
    }
    
    message("-")
    message("First Node: ", firstNode)
    message("Second Node: ", secondNode)
    # ---- Above loop to evaluate wrapper nodes ---
    firstNodeValid = (firstNode <= length(inputTree$tip.label) | firstNode %in% cladList ) #Set the first node as valid if either it had a tip label, or it has been assigned a clade
    message("first node valid: ", firstNodeValid)                               #Report if the node is valid
    
    secondNodeValid = (secondNode <= length(inputTree$tip.label) | secondNode %in% cladList )#Set the second node as valid if either it had a tip label, or it has been assigned a clade
    message("secondNodeValid: ", secondNodeValid)                               #Report if the node is valid
    
    #- Evaluate the node - 
    
    if((firstNodeValid & secondNodeValid)){                                     #If Both end nodes have been processed, evaluate this edge. If not, it will be skipped, and tried again on a subsequent pass, hopefully after the end nodes have been processed.
      remainingWrapperEdges = remainingWrapperEdges[!remainingWrapperEdges %in% i] #Remove this edge from the to-do list
      
      #- get the first node's name - 
      if(firstNode <= length(inputTree$tip.label)){                             #If the node has a tip label
        firstNodeName = inputTree$tip.label[firstNode]                          #The name is the tip label
      }else{                                                                    #if not
        firstNodeName = names(which(cladList == firstNode))                     #The name is the node's assigned clade name
      }
      message("firstNodeName: ", firstNodeName)                                 #Report the first node name
      
      #- get the second node's name -
      if(secondNode <= length(inputTree$tip.label)){                            #If the node has a tip label
        secondNodeName = inputTree$tip.label[secondNode]                        #The name is the tip label
      }else{                                                                    #if not
        secondNodeName = names(which(cladList == secondNode))                   #The name is the node's assigned clade name
      }
      message("secondNodeName: ", secondNodeName)                               #Report the second node name
      
      #Make the clade object
      assign(currentCladName, c(firstNodeName, secondNodeName))                 #assign those names to the current clade name
      
      #Add the clade to the list
      cladList = addToCladList(startNode)                                       #Add the clade to the cladelist
      cladNumber = cladNumber+1                                                 #update the clade number
      
      
    }else if ((!secondNodeValid | !firstNodeValid)){                            #Add a test for sitautions where an internal connection between foreground nodes is not itself foreground. Because this is such an odd edgecase, notify the user and do not attempt to resolve. 
      if(!firstNode %in% cladList){                                             #If the endNode is not in the clade list (this almost always happens because it on has one foreground child edge)
        if(all(firstNodeChildren %in% fgEdgeObjects)){                          #And both children are foreground
          message("The foreground node ", firstNode, "had two child nodes ", firstNodeChildren[1], " and ", firstNodeChildren[2], "but one of the connections between them is not foreground. Manual inspection reqiured.")  #Notify user. 
        }
      }else if(!secondNode %in% cladList){                                      #If the endNode is not in the clade list (this almost always happens because it on has one foreground child edge)
        if(all(secondNodeChildren %in% fgEdgeObjects)){                         #And both children are foreground
          message("The foreground node ", secondNode, "had two child nodes ", secondNodeChildren[1], " and ", secondNodeChildren[2], " which are both foreground. However, one of the connections between them is not foreground. Manual inspection reqiured.") #Notify user. 
        }
      }
      else{                                                                      #from before, for if one of the child nodes hasn't been processed yet
      message("Node not valid, skipping branch")                                #notify the user
      }
    }
    
    #check if it was stuck during this loop
    endRemainingEdges = length(remainingWrapperEdges)                           #get the length of the remaining edges now that the loop has run
    if(startRemainingEdges == endRemainingEdges){                               #see if it has not changed this loop
      stuckCycles = stuckCycles+1                                               #If it hasn't, increase the number of stuck cycles
    }else{                                                                      #If it has changed
      stuckCycles = 0                                                           #Set stuck cycles to zero 
    }
  }
}


# ---- generate a foreground string ----
foregroundNodes = which(1:1000 %in% as.vector(fgEdgeObjects))
foregroundStartNodes = foregroundNodes[foregroundNodes <= length(inputTree$tip.label)]
foregroundSpecies = inputTree$tip.label[foregroundStartNodes]

#TESTING CODE
#clade4[2] = "vs_mypDav1"
#clade6[2] = "prePar1"
#foregroundSpecies[12] = "vs_mypDav1"
#foregroundSpecies[17] = "prePar1"

foregroundSpeciesFilename = paste(outputFolderName, filePrefix, "ForegroundSpecies.rds", sep="")
saveRDS(foregroundSpecies, file = foregroundSpeciesFilename)




#

# ---- generate the sistersList ----
sisListFilename = paste(outputFolderName, filePrefix, "SistersList.rds", sep="")
cladObjectSet = ls(pattern = "clade")
sistersListExport =  mget(cladObjectSet)

# -- reorder the list to be in the "clade1, clade2" format --
orderedCladNames = "clade1"
sisterslistReOrder = sistersListExport[1]
for(i in 2:length(sistersListExport)){
  newCladNameObject = paste("clade", i, sep = '')
  orderedCladNames = append(orderedCladNames, newCladNameObject)
  sisterslistReOrder = append(sisterslistReOrder, sistersListExport[grep(orderedCladNames[i], names(sistersListExport))[1]])
}
sistersListExport = sisterslistReOrder
# --- save the list --
saveRDS(sistersListExport, file = sisListFilename)


#---- generate a number of internal nodes 
#internalNodes = foregroundNodes[foregroundNodes >= (length(inputTree$tip.label)+1)] #Get the number of internal nodes
#internalNodeNumber = length(internalNodes)
#message("internal node number: ", length(internalNodes))                        #message the number of internal nodes
#Save the internal node number 
#internalNodeFilename = paste(outputFolderName, filePrefix, "internalNodeNumber.rds", sep="")
#saveRDS(internalNodeNumber, file = internalNodeFilename)



#---- generate a phenotypeVector (named int of all tips with0/1 indicating foreground)-----
phenotypeVector = c(0,0)
length(phenotypeVector) = length(inputTree$tip.label)
phenotypeVector[] = 0 
names(phenotypeVector) = inputTree$tip.label
phenotypeVector[(names(phenotypeVector) %in% foregroundSpecies)] = 1
if(trimPhenotypeVector){
  phenotypeVector = phenotypeVector[names(phenotypeVector) %in% speciesFilter]
}
#Save the phenotypeVector
phenotypeVectorFilename = paste(outputFolderName, filePrefix, "phenotypeVector.rds", sep="")
saveRDS(phenotypeVector, file = phenotypeVectorFilename)




# --- Clades Correlation ---
#This correlation uses the Clades version of the path, and thus cannot be imported from the normal RER script.
#This is used in pValue calculations, so it is generated here if it does not already exist.

cladesPathsFileName = paste(outputFolderName, filePrefix, "CladesPathsFile.rds", sep= "")
cladesCorellationFileName = paste(outputFolderName, filePrefix, "CladesCorrelationFile", sep= "")

if(!file.exists(paste(cladesPathsFileName)) | !file.exists(paste(cladesCorellationFileName, ".rds", sep="")) | forceUpdate){

  # --Clades Paths --
  if(!file.exists(paste(cladesPathsFileName)) | forceUpdate){
    #get the main tree
    
    #make a foreground tree
    foregroundCladeTree = foreground2TreeClades(foregroundSpecies, sistersListExport, mainTrees, plotTree = T, )
    
    pathCladesObject = tree2PathsClades(foregroundCladeTree, mainTrees)
    saveRDS(pathCladesObject, file = cladesPathsFileName)
  }else{
    pathCladesObject = readRDS(cladesPathsFileName)
  }
  
  # -- Clades Correlations
  if(!file.exists(paste(cladesCorellationFileName, ".rds", sep="")) | forceUpdate){
    
    #Import RERs
    RERFileName = paste(outputFolderName, filePrefix, "RERFile.rds", sep= "")
    if(file.exists(paste(RERFileName))){
      RERObject = readRDS(RERFileName)
    }else{
      message("RERObject not found, generate an RER Object. (RunRERBinaryMT.R)")
    }
    
    cladesCorrelation = correlateWithBinaryPhenotype(RERObject, pathCladesObject, min.sp =35)
    write.csv(cladesCorrelation, file= paste(cladesCorellationFileName, ".csv", sep =""), row.names = T, quote = F)
    saveRDS(cladesCorrelation, file= paste(cladesCorellationFileName, ".rds", sep=""))
  }
}

#

# -- Generate a readable tree pdf -- 

#dev.off(); dev.new(); dev.new(); testplot2 = plotTreeHighlightBranches(inputTree,hlspecies=which(inputTree$edge.length== 3),hlcols="blue", main="Marine mammals trait tree"); edgelabels(cex = 0.7, frame="none", font=2, adj=c(0,-0.2), col="blue"); nodelabels(cex = 0.7, frame="none", font=2, adj=c(-0.2,0.3), col="dark green"); tiplabels(cex = 0.8, frame="none", font=2, adj=c(0.2,0), col="dark red")
#testplot2 = plotTreeHighlightBranches(inputTree,hlspecies=which(inputTree$edge.length== 3),hlcols="blue", main="Marine mammals trait tree");


#testTreeDisplayable = inputTree
#replace(testTreeDisplayable$edge.length, testTreeDisplayable$edge.length==0, 0.5)
#replace(testTreeDisplayable$edge.length, testTreeDisplayable$edge.length==1, 4)

# plotTreeHighlightBranches(testTreeDisplayable, hlspecies=which(inputTree$edge.length==1), hlcols="blue",main="PhenotypeVector input tree")
# 
# 
# 
# 
# 
# plotTreeHighlightBranches(testTreeDisplayable, hlspecies= fgEdges, hlcols="blue",main="PhenotypeVector input tree")
# plotTreeHighlightBranches(testTreeDisplayable, hlspecies= WrapperEdges, hlcols="blue",main="PhenotypeVector input tree")
# 
# debugPlotTree(inputTree)
# 
# foregroundSpecies
# sistersListExport
# testFGTree = foreground2TreeClades(foregroundSpecies, sistersListExport, mainTrees)
# marineplot1 = plotTreeHighlightBranches(testFGTree,
#                                         hlspecies=which(testFGTree$edge.length==1),
#                                         hlcols="blue", main="Marine mammals trait tree")
# if(trimPhenotypeVector){
#   unprunedTestFGTree = testFGTree
#   testFGTree = pruneTree(testFGTree, speciesFilter)
# }
# testFGDisplay = plotTreeHighlightBranches(testFGTree, hlspecies=which(testFGTree$edge.length==1),hlcols="blue")
# edgelabels(cex = 0.7, frame="none", font=2, adj=c(0,-0.2), col="blue")
# 
# 
# debugPlotTree(testFGTree)
# 
# slexport1 = sistersListExport
# all.equal(slexport1, sistersListExport)

#
#manualSistersList = list(clade1, clade2, clade3, clade4, clade5, clade6, clade7, clade8, clade9, clade10, clade11, clade12,clade13,clade14,clade15,clade16,clade17,clade18,clade19,clade20)


#foregroundString = foregroundSpecies
#sistersList = manualSistersList


#sistersListExport[!sistersListExport %in% manualSistersList]
#setdiff(sistersListExport, manualSistersList)

#sisterslistReOrder = sistersListExport[1]
#sisterslistReOrder = append(sisterslistReOrder, sistersListExport[12:19])
#sisterslistReOrder = append(sisterslistReOrder, sistersListExport[2:11])
#orderedCladNames = c("clade1", "clade2", "clade3", "clade4","clade5","clade6","clade7","clade8","clade9","clade10","clade11","clade12","clade13","clade14","clade15","clade16","clade17","clade18","clade19","clade20","clade21","clade22","clade23","clade24","clade25","clade26","clade27","clade28","clade29","clade30","clade31","clade32","clade33","clade34","clade35","clade36","clade37","clade38","clade39","clade40","clade41","clade42","clade43","clade44","clade45","clade46","clade47","clade48","clade49","clade50")
#orderedCladNames = "clade1"
#sisterslistReOrder = sistersListExport[1]
#for(i in 2:length(sisterListExport)){
#  newCladNameObject = paste("clade", i, sep = '')
#  orderedCladNames = append(orderedCladNames, newCladNameObject)
#  sisterslistReOrder = append(sisterslistReOrder, sistersListExport[grep(orderedCladNames[i], names(sistersListExport))])
#}

